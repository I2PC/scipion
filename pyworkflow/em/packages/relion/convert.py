# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains converter functions that will serve to:
1. Write from base classes to Relion specific files
2. Read from Relion files to base classes
"""

import os
from os.path import join
from pyworkflow.utils import Environ
from pyworkflow.utils.path import createLink, cleanPath
from pyworkflow.em import ImageHandler
from pyworkflow.em.data import SetOfClasses2D, SetOfClasses3D, SetOfParticles
from constants import *
# Since Relion share several conventions with Xmipp, we will reuse 
# the xmipp3 tools implemented for Scipion here
import xmipp


def getEnviron():
    """ Setup the environment variables needed to launch Relion. """
    environ = Environ(os.environ)
    environ.update({
            'PATH': join(os.environ['RELION_HOME'], 'bin'),
            'LD_LIBRARY_PATH': join(os.environ['RELION_HOME'], 'lib') + ":" + join(os.environ['RELION_HOME'], 'lib64'),
#             'LD_LIBRARY_PATH': join(os.environ['RELION_HOME'], 'lib64'),
            }, position=Environ.BEGIN)
    return environ


def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    from pyworkflow.em.packages.xmipp3.convert import locationToXmipp
    return locationToXmipp(index, filename)


def relionToLocation(filename):
    from pyworkflow.em.packages.xmipp3.convert import xmippToLocation
    return xmippToLocation(filename)


_xmippLabelsDict = {} # Dictionary to store mappings replaced

            
def restoreXmippLabels():
    global _xmippLabelsDict
    for k, v in _xmippLabelsDict.iteritems():
        xmipp.addLabelAlias(k, v, True)
    _xmippLabelsDict = {}


class ParticleAdaptor():
    """ Class used to convert a set of particles for Relion.
    It will write an stack in Spider format and also
    modify the output star file to point to the new stack.
    """
    def __init__(self, imgSet, filesMapping={}, originalSet=None, 
                 preprocessImageRow=None,
                 postprocessImageRow=None):
        self._rowCount = 1
        self._ih = ImageHandler()
        self._imgSet = imgSet
        self._filesMapping = filesMapping
        self._originalSet = originalSet
        self._preprocessImageRow = preprocessImageRow
        self._postprocessImageRow = postprocessImageRow
        
    def preprocessImageRow(self, img, imgRow):
        """ If binary images where converted or just 
        created links, replace the orignal names.
        """
        # preprocessImageRow is a way relion protocol have to
        # do some preprocessing on image object or row 
        # before the normal setup
        if self._preprocessImageRow:
            self._preprocessImageRow(img, imgRow)
            
        # Re-write the row with the new location
        if self._filesMapping:
            # If the filesMapping is not empty
            # we assume that the stack files have been converted
            # to be properly read by relion and also that
            # each particle preserve its index and only change the filename
            newFile = self._filesMapping[img.getFileName()]
            newLoc = (img.getIndex(), newFile)
            img.setLocation(newLoc)

    def postprocessImageRow(self, img, imgRow):
        """ If binary images where converted or just 
        created links, replace the orignal names.
        """            
        if img.hasMicId():
            imgRow.setValue('rlnMicrographName', 'fake_micrograph_%06d.mrc' % img.getMicId())
            imgRow.setValue(xmipp.MDL_MICROGRAPH_ID, long(img.getMicId()))
            
        if imgRow.hasLabel(xmipp.MDL_PARTICLE_ID):
            # Write stuff for particle polishing
            movieId = imgRow.getValue(xmipp.MDL_MICROGRAPH_ID)
            movieName = 'fake_micrograph_%06d_movie.mrcs' % movieId 
            frameId = imgRow.getValue(xmipp.MDL_FRAME_ID) + 1
            micName = '%06d@%s' % (frameId, movieName)
            particleId = imgRow.getValue(xmipp.MDL_PARTICLE_ID)
            particle = self._originalSet[particleId]
            if particle is None:
                particleName = 'None'
                raise Exception("Particle with id %d not found!!!" % particleId)
            else:
                particleName = locationToRelion(*particle.getLocation())
                
            imgRow.setValue('rlnMicrographName', micName)
            imgRow.setValue('rlnParticleName', particleName)
            
        for label, _ in imgRow:
            if not label in XMIPP_RELION_LABELS:
                imgRow.removeLabel(label)
        self._rowCount += 1
        
        if self._postprocessImageRow:
            self._postprocessImageRow(img, imgRow)
    

def addRelionLabels(replace=False, extended=False):
    """ Add relion labels as aliases for Xmipp metadata. """
    global _xmippLabelsDict
    _xmippLabelsDict = {}
    for k, v in XMIPP_RELION_LABELS.iteritems():
        _xmippLabelsDict[k] = xmipp.label2Str(k) # store original label string
        xmipp.addLabelAlias(k, v, replace)
    return
    if extended:
        for k, v in XMIPP_RELION_LABELS_EXTRA.iteritems():    
            _xmippLabelsDict[k] = xmipp.label2Str(k) # store original label string
            xmipp.addLabelAlias(k, v, replace)

            
def addRelionLabelsToEnviron(env):
    """ create an string that can be used for XMIPP_EXTRA_ALIASES
    for adding the labels of Relion.
    """
    pairs = []
    for k, v in XMIPP_RELION_LABELS.iteritems():
        pairs.append('%s=%s' % (xmipp.label2Str(k), v))
    return 
    for k, v in XMIPP_RELION_LABELS_EXTRA.iteritems():
        pairs.append('%s=%s' % (xmipp.label2Str(k), v))        
    varStr = ';'.join(pairs)
    env['XMIPP_EXTRA_ALIASES'] = varStr
    
    
def readSetOfParticles(filename, partSet, **kwargs):
    import pyworkflow.em.packages.xmipp3 as xmipp3
    addRelionLabels(replace=True)
    xmipp3.readSetOfParticles(filename, partSet, **kwargs)
    restoreXmippLabels()


def writeSetOfParticles(imgSet, starFile,
                        outputDir,
                        originalSet=None, **kwargs):
    """ This function will write a SetOfImages as Relion metadata.
    Params:
        imgSet: the SetOfImages instance.
        starFile: the filename where to write the metadata.
        filesMapping: this dict will help when there is need to replace images names
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    addRelionLabels(replace=True)
    filesMapping = convertBinaryFiles(imgSet, outputDir)
    pa = ParticleAdaptor(imgSet, filesMapping, originalSet,
                         preprocessImageRow=kwargs.get('preprocessImageRow', None),
                         postprocessImageRow=kwargs.get('postprocessImageRow', None))
    # Intercept the pre and post processImageRow
    kwargs['preprocessImageRow'] = pa.preprocessImageRow
    kwargs['postprocessImageRow'] = pa.postprocessImageRow
    xmipp3.writeSetOfParticles(imgSet, starFile, **kwargs)
    restoreXmippLabels()


def createClassesFromImages(inputImages, inputStar, classesFn, ClassType, 
                            classLabel, classFnTemplate, iter, preprocessRow=None):
    """ From an intermediate dataXXX.star file produced by relion, create
    the set of classes in which those images are classified.
    Params:
        inputImages: the SetOfImages that were initially classified by relion. 
        inputStar: the filename with relion star format.
        classesFn: filename where to write the classes.
        ClassType: the output type of the classes set ( usually SetOfClass2D or 3D )
        classLabel: label that is the class reference.
        classFnTemplate: the template to get the classes averages filenames
        iter: the iteration number, just used in Class template
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    addRelionLabels(replace=True, extended=True)
    # We asume here that the volumes (classes3d) are in the same folder than imgsFn
    # rootDir here is defined to be used expanding locals()
    xmipp3.createClassesFromImages(inputImages, inputStar, classesFn, ClassType, 
                                   classLabel, classFnTemplate, iter, preprocessRow)
    restoreXmippLabels()    
    
#def createClassesFromImages(inputImages, inputStar, classesFn, ):
    

def readSetOfClasses2D(classes2DSet, filename, **args):
    clsSet = SetOfClasses2D(filename=filename)
    classes2DSet.appendFromClasses(clsSet)


def readSetOfClasses3D(classes3DSet, filename, **args):
    clsSet = SetOfClasses3D(filename=filename)
    classes3DSet.appendFromClasses(clsSet)
    

def writeSqliteIterData(imgStar, imgSqlite):
    """ Given a Relion images star file (from some iteration)
    create the corresponding SetOfParticles (sqlite file)
    for this iteration. This file can be visualized sorted
    by the LogLikelihood.
    """
    addRelionLabels(replace=True, extended=True)
    cleanPath(imgSqlite)
    imgSet = SetOfParticles(filename=imgSqlite)
    readSetOfParticles(imgStar, imgSet)
    imgSet.write()
    restoreXmippLabels()
    
    
def writeIterAngularDist(self, inDataStar, outAngularDist, numberOfClasses, prefixes):
    """ Write the angular distribution. Should be overriden in subclasses. """
    addRelionLabels(replace=True, extended=True)
    
    md = xmipp.MetaData(inDataStar)
    
    refsList = range(1, numberOfClasses+1) 
    print "refsList: ", refsList
    for ref3d in refsList:
        for prefix in prefixes:
            mdGroup = xmipp.MetaData()
            mdGroup.importObjects(md, xmipp.MDValueEQ(xmipp.MDL_REF, ref3d))
            mdDist = xmipp.MetaData()
            mdDist.aggregateMdGroupBy(mdGroup, xmipp.AGGR_COUNT, 
                                      [xmipp.MDL_ANGLE_ROT, xmipp.MDL_ANGLE_TILT], 
                                      xmipp.MDL_ANGLE_ROT, xmipp.MDL_WEIGHT)
            mdDist.setValueCol(xmipp.MDL_ANGLE_PSI, 0.0)
            blockName = '%sclass%06d_angularDist@' % (prefix, ref3d)
            print "Writing angular distribution to: ", blockName + data_angularDist
            mdDist.write(blockName + data_angularDist, xmipp.MD_APPEND) 
            
    restoreXmippLabels()
    
    
def splitInCTFGroups(self, imgStar, groups):
    """ Add a new colunm in the image star to separate the particles into ctf groups """
    addRelionLabels(replace=True)
    
    md = xmipp.MetaData(imgStar)
    defocus = [md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId) for objId in md]
    
    minDef = min(defocus)
    maxDef = max(defocus)
    
    increment = (maxDef - minDef) / groups
    counter = 1
    part = 1
    md.sort(xmipp.MDL_CTF_DEFOCUSU)
    
    for objId in md:
        partDef = md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId)
        currDef = minDef + increment
        
        if partDef <= currDef:
            part = part + 1
            md.setValue(xmipp.MDL_SERIE, "ctfgroup_%05d" % counter, objId)
        else:
            if part < 100:
                increment = md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId) - minDef
                part = part + 1
                md.setValue(xmipp.MDL_SERIE, "ctfgroup_%05d" % counter, objId)
            else:
                part = 1
                minDef = md.getValue(xmipp.MDL_CTF_DEFOCUSU, objId)
                counter = counter + 1
                increment = (maxDef - minDef) / (groups - counter)
                md.setValue(xmipp.MDL_SERIE, "ctfgroup_%05d" % counter, objId)
    md.write(imgStar)
    restoreXmippLabels()
    
        
def getMdFirstRow(starfile):
    addRelionLabels(replace=True, extended=True)
    from pyworkflow.em.packages.xmipp3.utils import getMdFirstRow
    row = getMdFirstRow(starfile)
    restoreXmippLabels()
    return row


def findImagesPath(starFile):
    """ Find the path of the images relative to some star file. """
    absPath = os.path.dirname(os.path.abspath(starFile))
    row = getMdFirstRow(starFile)
    _, imgFile = relionToLocation(row.getValue(xmipp.MDL_IMAGE))
    
    while absPath is not None and absPath != '/':
        if os.path.exists(os.path.join(absPath, imgFile)):
            return os.path.relpath(absPath)
        absPath = os.path.dirname(absPath)
        
    return None


def prependToFileName(imgRow, prefixPath):
    """ Prepend some root name to imageRow filename. """
    index, imgPath = relionToLocation(imgRow.getValue(xmipp.MDL_IMAGE))
    newLoc = locationToRelion(index, os.path.join(prefixPath, imgPath))
    imgRow.setValue(xmipp.MDL_IMAGE, newLoc)


def relativeFromFileName(imgRow, prefixPath):
    """ Remove some prefix from filename in row. """
    index, imgPath = relionToLocation(imgRow.getValue(xmipp.MDL_IMAGE))
    imgPath = os.path.relpath(imgPath, prefixPath)
    newLoc = locationToRelion(index, imgPath)
    imgRow.setValue(xmipp.MDL_IMAGE, newLoc)


def setupCTF(imgRow, sampling):
    """ Do some validations and set some values
    for Relion import.
    """
    imgRow.setValue(xmipp.MDL_SAMPLINGRATE, sampling)
    #TODO: check if we want to move this behaviour to setup CTFModel by default
    hasDefocusU = imgRow.containsLabel(xmipp.MDL_CTF_DEFOCUSU)
    hasDefocusV = imgRow.containsLabel(xmipp.MDL_CTF_DEFOCUSV)
    hasDefocusAngle = imgRow.containsLabel(xmipp.MDL_CTF_DEFOCUS_ANGLE)
    
    if hasDefocusU or hasDefocusV:
        if not hasDefocusU:
            imgRow.setValue(xmipp.MDL_CTF_DEFOCUSU, imgRow.getValue(xmipp.MDL_CTF_DEFOCUSV))
        if not hasDefocusV:
            imgRow.setValue(xmipp.MDL_CTF_DEFOCUSV, imgRow.getValue(xmipp.MDL_CTF_DEFOCUSU))
        if not hasDefocusAngle:
            imgRow.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, 0.)
            

def convertBinaryFiles(imgSet, outputDir):
    """ Convert binary images files to a format read by Relion.
    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesMapping = {}
    # This approach can be extended when
    # converting from a binary file format that
    # is not read from Relion
    ext = imgSet.getFirstItem().getFileName()
    if ext.endswith('.mrc'):
        uniqueFiles = list(imgSet.getFiles())
        for f in uniqueFiles:
            newFile = os.path.join(outputDir, os.path.basename(f)+'s')
            createLink(f, newFile)
            filesMapping[f] = newFile
    elif ext.endswith('.hdf'):#assume eman
        uniqueFiles = list(imgSet.getFiles())
        from pyworkflow.em.convert import convertStack
        for f in uniqueFiles:
            newFile = os.path.join(outputDir, os.path.basename(f.rsplit( ".", 1 )[ 0 ])+'.spi')
            convertStack(f, newFile)
            filesMapping[f] = newFile

    return filesMapping
