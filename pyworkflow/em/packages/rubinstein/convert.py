from os.path import join, basename, exists
import numpy
from collections import OrderedDict
from pyworkflow.utils.path import (createLink, cleanPath, copyFile,
                                   replaceBaseExt, getExt, removeExt)
import pyworkflow.em as em
from pyworkflow.em.data import Micrograph
import pyworkflow.em.metadata as md
from pyworkflow.object import *

ACQUISITION_DICT = OrderedDict([
       ("_amplitudeContrast", md.RLN_CTF_Q0),
       ("_sphericalAberration", md.RLN_CTF_CS),
       ("_voltage", md.RLN_CTF_VOLTAGE),
        ("_magnification", md.RLN_CTF_MAGNIFICATION)
       ])
COOR_DICT = OrderedDict([
    ("_x", md.RLN_IMAGE_COORD_X),
    ("_y", md.RLN_IMAGE_COORD_Y)
])

CTF_DICT = OrderedDict([
       ("_defocusU", md.RLN_CTF_DEFOCUSU),
       ("_defocusV", md.RLN_CTF_DEFOCUSV),
       ("_defocusAngle", md.RLN_CTF_DEFOCUS_ANGLE)
       ])

ALIGNMENT_DICT = OrderedDict([
       ("_rlnOriginX", md.RLN_ORIENT_ORIGIN_X),
       ("_rlnOriginY", md.RLN_ORIENT_ORIGIN_Y),
       ("_rlnOriginZ", md.RLN_ORIENT_ORIGIN_Z),
       ("_rlnAngleRot", md.RLN_ORIENT_ROT),
       ("_rlnAngleTilt", md.RLN_ORIENT_TILT),
       ("_rlnAnglePsi", md.RLN_ORIENT_PSI),
       ])

def writeSetOfParticles(imgSet, starFile,
                        outputDir, **kwargs):
    """ This function will write a SetOfImages as meta
    Params:
        imgSet: the SetOfImages instance.
        starFile: the filename where to write the meta
        filesMapping: this dict will help when there is need to replace images names
    """
    filesDict = convertBinaryFiles(imgSet, outputDir)
    kwargs['filesDict'] = filesDict
    partMd = md.MetaData()
    setOfImagesToMd(imgSet, partMd, particleToRow, **kwargs)
    blockName = kwargs.get('blockName', 'Particles')
    partMd.write('%s@%s' % (blockName, starFile))

def convertBinaryFiles(imgSet, outputDir, extension='mrcs'):
    """ Convert binary images files to a format read by alignparts.
    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, no conversion was done.
    """
    filesDict = {}
    ih = em.ImageHandler()

    def getUniqueFileName(fn, extension):
        """ Get an unique file for either link or convert files.
        It is possible that the base name overlap if they come
        from different runs. (like particles.mrcs after relion preprocess)
        """
        newFn = join(outputDir, replaceBaseExt(fn, extension))
        newRoot = removeExt(newFn)

        values = filesDict.values()
        counter = 1

        while newFn in values:
            counter += 1
            newFn = '%s_%05d.%s' % (newRoot, counter, extension)

        return newFn

    def createBinaryLink(fn):
        """ Just create a link named .mrcs to Relion understand
        that it is a binary stack file and not a volume.
        """
        newFn = getUniqueFileName(fn, extension)
        createLink(fn, newFn)
        return newFn

    def convertStack(fn):
        """ Convert from a format that is not read by Relion
        to an spider stack.
        """
        newFn = getUniqueFileName(fn, 'stk')
        ih.convertStack(fn, newFn)
        return newFn

    ext = getExt(imgSet.getFirstItem().getFileName())[1:]  # remove dot in extension

    if ext == extension:
        mapFunc = createBinaryLink
        print "convertBinaryFiles: creating soft links."
    elif ext == 'mrc' and extension == 'mrcs':
        mapFunc = createBinaryLink
        print "convertBinaryFiles: creating soft links (mrcs -> mrc)."
    elif ext.endswith('hdf'):  # assume eman .hdf format
        mapFunc = convertStack
        print "convertBinaryFiles: converting stacks. (%s -> %s)" % (extension, ext)
    else:
        mapFunc = None

    if mapFunc is not None:
        for fn in imgSet.getFiles():
            newFn = mapFunc(fn)  # convert or link
            filesDict[fn] = newFn  # map new filename
            print "   %s -> %s" % (newFn, fn)

    return filesDict

def setOfImagesToMd(imgSet, imgMd, imgToFunc, **kwargs):
    """ This function will fill metadata from a SetOfMicrographs
    Params:
        imgSet: the set of images to be converted to metadata
        md: metadata to be filled
        rowFunc: this function can be used to setup the row before
            adding to meta
    """
    if 'alignType' not in kwargs:
        kwargs['alignType'] = imgSet.getAlignment()

    for img in imgSet:
        objId = imgMd.addObject()
        imgRow = md.Row()
        imgToFunc(img, imgRow, **kwargs)
        imgRow.writeToMd(imgMd, objId)

def particleToRow(part, partRow, **kwargs):
    """ Set labels values from Particle to md row. """
    coord = part.getCoordinate()
    if coord is not None:
        coordinateToRow(coord, partRow, copyId=False)
    if part.hasMicId():
        partRow.setValue(md.RLN_MICROGRAPH_ID, long(part.getMicId()))
    imageToRow(part, partRow, md.RLN_IMAGE_NAME, **kwargs)

def coordinateToRow(coord, coordRow, copyId=True):
    """ Set labels values from Coordinate coord to md row. """
    if copyId:
        setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT)
    if coord.getMicId():
        coordRow.setValue(md.RLN_MICROGRAPH_NAME, str(coord.getMicName()))

def objectToRow(obj, row, attrDict):
    """ This function will convert an EMObject into a XmippMdRow.
    Params:
        obj: the EMObject instance (input)
        row: the XmippMdRow instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    if obj.isEnabled():
        enabled = True
    else:
        enabled = False
    row.setValue(md.RLN_IMAGE_ENABLED, enabled)

    for attr, label in attrDict.iteritems():
        if hasattr(obj, attr):
            valueType = md.label2Python(label)
            row.setValue(label, valueType(getattr(obj, attr).get()))

def imageToRow(img, imgRow, imgLabel=md.RLN_IMAGE_NAME, **kwargs):
    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, imgRow)
    setRowId(imgRow, img)  # Set the id in the metadata as MDL_ITEM_ID
    index, fn = img.getLocation()
    # check if the is a file mapping
    filesDict = kwargs.get('filesDict', {})
    filename = filesDict.get(fn, fn)
    imgRow.setValue(imgLabel, filenamewithindex(index, filename))

    if kwargs.get('writeCtf', True) and img.hasCTF():
        ctfModelToRow(img.getCTF(), imgRow)

    if kwargs.get('writeAcquisition', True) and img.hasAcquisition():
        acquisitionToRow(img.getAcquisition(), imgRow)

    objectToRow(img, imgRow, {})

    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, imgRow)

def setRowId(mdRow, obj, label=md.RLN_IMAGE_ID):
    mdRow.setValue(label, long(obj.getObjId()))

def filenamewithindex(index, filename):
    if index != em.NO_INDEX:
        return "%06d@%s" % (index, filename)
    return filename

def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """
    objectToRow(ctfModel, ctfRow, CTF_DICT)

def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)
