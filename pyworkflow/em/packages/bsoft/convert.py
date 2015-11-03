'''
Created on Feb 21, 2014

@author: airen
'''
import os
from xmipp import *
from pyworkflow.utils.path import join, dirname, replaceBaseExt
from pyworkflow.em import *
from pyworkflow.em.packages.xmipp3.convert import rowToCoordinate, rowFromMd


def getEnviron():
    """ Setup the environment variables needed to launch Bsoft. """
    environ = Environ(os.environ)
    BSOFT_HOME = os.environ['BSOFT_HOME']
    
    environ.update({
            'BSOFT': BSOFT_HOME,
            'BPARAM': join(BSOFT_HOME, 'params'),
            'PATH': join(BSOFT_HOME, 'bin'),
            'LD_LIBRARY_PATH': join(BSOFT_HOME, 'lib')
            }, position=Environ.BEGIN)
    return environ

# Map from Xmipp labels to Relion labels names
XMIPP_BSOFT_LABELS = {
                        MDL_MICROGRAPH:        'micrograph.file_name'
                       ,MDL_MICROGRAPH_ID:        'micrograph.id' 
                       ,MDL_IMAGE:        'particle.filename'
                       ,MDL_PARTICLE_ID:        'particle.id'           
                       ,MDL_XCOOR:         'particle.x'
                       ,MDL_YCOOR:        'particle.y'
                       ,MDL_ZCOOR:        'particle.z'
                       ,MDL_SHIFT_X:    'particle.origin_x'
                       ,MDL_SHIFT_Y:    'particle.origin_y'
                       ,MDL_SHIFT_Z:    'particle.origin_z'
                       ,MDL_ENABLED:    'particle.select'
                       ,MDL_PICKING_PARTICLE_SIZE:   'particle.origin_x'
                       ,MDL_MAGNIFICATION:  'particle.magnification'
                       
                       }



def addBsoftLabelAliases():
    for k, v in XMIPP_BSOFT_LABELS.iteritems():
        addLabelAlias(k, v, True)
        
_xmippLabelsDict = {} # Dictionary to store mappings replaced

def restoreXmippLabels():
    global _xmippLabelsDict
    for k, v in _xmippLabelsDict.iteritems():
        xmipp.addLabelAlias(k, v, True)
    _xmippLabelsDict = {}

        
def readSetOfCoordinates(outputDir, micSet, coordSet):
    """ Read from Bsoft .star files.
    Params:
        outputDir: the directory where the .star files are.
           
        micSet: the SetOfMicrographs to associate the .star, which 
            name should be the same of the micrographs.
        coordSet: the SetOfCoordinates that will be populated.
    """

    addBsoftLabelAliases()
    boxSize = None
    for mic in micSet:
        outputFile = join(outputDir, replaceBaseExt(mic.getFileName(), 'star'))
        #scipionPosFile = join(outputDir, "scipion_" + replaceBaseExt(mic.getFileName(), 'pos'))
        if exists(outputFile):
            posMd = xmipp.MetaData(outputFile)
        else:
            posMd = xmipp.MetaData()
        
        for objId in posMd:
            coord = rowToCoordinate(rowFromMd(posMd, objId))
            boxSize = 2 * posMd.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, objId)
            coord.setMicrograph(mic)
            coord.setX(coord.getX())
            coord.setY(coord.getY())
            
            coordSet.append(coord)      
            # Add an unique ID that will be propagated to particles
            posMd.setValue(xmipp.MDL_ITEM_ID, long(coord.getObjId()), objId)
#         if not posMd.isEmpty():
#             posMd.write("particles@%s"  % scipionPosFile)
            
            
             #reading origin.x value and converting to particle size, can change, we take last value
            coordSet.setBoxSize(boxSize)
            
            
class ParticleAdaptor():
    """ Class used to convert a set of particles for Bsoft.
    It will write an stack in Spider format and also
    modify the output star file to point to the new stack.
    """
    def __init__(self, imgSet, stackFile=None):
        self._rowCount = 1
        self._ih = ImageHandler()
        self._imgSet = imgSet
        self._stackFile = stackFile
        
        import pyworkflow.em.packages.xmipp3 as xmipp3
        self._particleToRow = xmipp3.particleToRow
        
    def setupRow(self, img, imgRow):
        """ Convert image and modify the row. """
        if self._stackFile is not None:
            newLoc = (self._rowCount, self._stackFile)
            #TODO: Check whether the input image format is valid for Bsoft
            self._ih.convert(img.getLocation(), newLoc)
            img.setLocation(newLoc)
            # Re-write the row with the new location
        self._particleToRow(img, imgRow) #TODO: CHECK why not the following, writeAlignment=False)
            
        for label, _ in imgRow:
            if not label in XMIPP_BSOFT_LABELS:
                imgRow.removeLabel(label)
        self._rowCount += 1  
   

def writeSetOfParticles(imgSet, starFile, stackFile):
    """ This function will write a SetOfImages as Bsoft metadata.
    Params:
        imgSet: the SetOfImages instance.
        starFile: the filename where to write the metadata.
    """
    import pyworkflow.em.packages.xmipp3 as xmipp3
    print "In writeSetOfParticles Bsoft"
    addBsoftLabelAliases()
#     pa = ParticleAdaptor(imgSet, stackFile)
#     xmipp3.writeSetOfParticles(imgSet, starFile, rowFunc=pa.setupRow, writeAlignment=False)
    md = xmipp.MetaData()
    md.setColumnFormat(False)
    blockName = 'stack'
    imgRow = XmippMdRow()
    imgRow.setValue(MDL_MICROGRAPH_ID, long(1))
    imgRow.setValue(MDL_IMAGE, str(stackFile))
    imgRow.writeToMd(md, md.addObject())  
    imgSet._bsoftStar = String(starFile)
    restoreXmippLabels()
    
def createBsoftInputParticles(imgSet, starFile, stackFile): 
    """ Ensure that in 'starFile' it is a valid STAR files with particles.
    If the imgSet comes from Bsoft, just create a link.
    If not, then write the proper file.
    """
    imgsStar = getattr(imgSet, '_bsoftStar', None)
    if imgsStar is None:
        writeSetOfParticles(imgSet, starFile, stackFile)
    else:
        imgsFn = imgsStar.get()
        createLink(imgsFn, imgsStar.get())