# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This modules contains basic hierarchy
for specific Xmipp3 EM data objects
"""

from pyworkflow.em import *
from xmipp3 import XmippMdRow, XmippSet, findRowById

from pyworkflow.utils.path import replaceBaseExt, exists, dirname, join
import xmipp
    
      
class XmippImageLocation(ImageLocation):
    
    def __str__(self):
        if self._index:
            return "%d@%s" % (self._index, self._filename)
        else:
            return self._filename
    
    def parse(self, locStr):
        self._index, self._filename = xmipp.FileName(locStr).decompose()        


class XmippImage(XmippMdRow, Image):
    """Xmipp implementation for Image"""
    _label = xmipp.MDL_IMAGE
    _labelCTF = xmipp.MDL_CTF_MODEL
    _labelId = xmipp.MDL_ITEM_ID
    
    def __init__(self, filename=None, **args):
        XmippMdRow.__init__(self)
        
        self.setId(0) # By default value 0 means no Id
        Image.__init__(self, filename, **args)
    
    def getLabelValue(self):
        return 
    
    def getLocation(self):
        return XmippImageLocation(locStr=self.getValue(self._label))
    
    def setLocation(self, loc):
        locStr = str(XmippImageLocation(loc=loc))
        self.setValue(self._label, locStr)
        
    def setFileName(self, filename):
        if filename is not None:
            self.setLocation(XmippImageLocation(filename=filename))
        
    def getFileName(self):
        return self.getLocation().getFileName()
    
    def hasCTF(self):
        return self.hasLabel(self._labelCTF)
    
    def getCTF(self):
        """ Return the CTF model """
        return XmippCTFModel(self.getValue(self._labelCTF))
    
    def setCTF(self, newCTF):
        # TODO: implement the set of CTF, probably need a conversion step
        raise Exception("Not implemented setCTF for XmippImage")
    
    def getId(self):
        return self.getValue(self._labelId)
        
    def setId(self, imgId):
        """ This id identifies the element inside a set """
        self.setValue(self._labelId, long(imgId))
                
    @staticmethod
    def _convert(itemClass, img, ctfFn):
        """Convert from Micrograph to XmippMicrograph"""
        if isinstance(img, itemClass):
            return img
        imgXmipp = itemClass(img.getFileName())
        if img.hasCTF():
            ctf = img.getCTF()
            ctfXmipp = XmippCTFModel.convert(ctf, ctfFn)
            imgXmipp.setValue(xmipp.MDL_CTF_DEFOCUSU, ctf.defocusU.get())
            imgXmipp.setValue(xmipp.MDL_CTF_DEFOCUSV, ctf.defocusV.get())

            imgXmipp.setValue(xmipp.MDL_IMAGE1, ctf.getPsdFile())
            imgXmipp.setValue(xmipp.MDL_CTF_MODEL, ctfFn)
            
        if img.hasId():
            imgXmipp.setId(img.getId())
        else:
            raise Exception("XmippImage conversion failed: Image id not found.")
        # TODO: copyInfo??
        # from img to imgXmipp??  
        return imgXmipp
    
    @staticmethod
    def convert(img, ctfFn):
        return XmippImage._convert(XmippImage, img, ctfFn)
        
        
        
        
class XmippSetOfImages(SetOfImages):
    """Represents a set of Images for Xmipp"""
    _label = xmipp.MDL_IMAGE
    
    def __init__(self, filename=None, **args):
        SetOfImages.__init__(self, filename, **args)
        self._set = XmippSet(XmippImage)
        self._setPairs = XmippSet(XmippTiltedPair)
           
    def getSize(self):
        """Return the number of images"""
        self.loadIfEmpty()
        return self._set.getSize()
    
    def load(self):
        """ Load extra data from files. """
        if self.getFileName() is None:
            raise Exception("Set filename before calling load()")
        self._set.read(self._getListBlock())
        if self.hasTiltPairs():
            self._setPairs.read(self._getTiltedBlock())
        
    def loadIfEmpty(self):
        """ Load data only if the main set is empty. """
        if self._set.isEmpty():
            self.load()
        
    def sort(self):
        """Sort the set according to MDL_IMAGE"""
        self._set.sort(self._label)
        
    def write(self):
        """ Write metadatas to file """
        if self.getFileName() is None:
            raise Exception("Set filename before calling load()")
        self._set.write(self._getListBlock(), xmipp.MD_OVERWRITE)
        if self.hasTiltPairs():
            self._setPairs.write(self._getTiltedBlock(), xmipp.MD_APPEND)   
         
    def append(self, img):
        """Add a new micrograph to the set"""          
        if not isinstance(img, XmippImage):
            path = dirname(self.getFileName())
            name = replaceBaseExt(img.getFileName(), 'ctfparam')
            ctfFn = join(path, 'extra', name) # TODO: check later if we want the 'extra' logic here
            xmippImg = self._set.getItemClass().convert(img, ctfFn)
        else:
            xmippImg = img
            if not img.hasId():
                self._idCount += 1
                xmippImg.setId(self._idCount)
        
        self._set.append(xmippImg)
        
    def appendFromMd(self, md):
        """ Add all entries found in the md. """
        self._set._md.unionAll(md)   
        
    def appendPair(self, idU, idT):
        """ Add a new tilted pair to the set"""
        self._setPairs.append(XmippTiltedPair(idU, idT))
        
    def __iter__(self):
        """ Iterate over the set of images. """
        self.loadIfEmpty()
        for img in self._set:
            img.setSamplingRate(self.samplingRate.get())
            yield img
        
    def iterTiltPairs(self):
        """Iterate over the tilt pairs if is the case"""
        #self.__loadFiles()        
        self.loadIfEmpty()
        for tp in self._setPairs:
            #retrieve index for untilted and tilted micrographs
            yield (tp.getUntilted(), tp.getTilted())
            
    def __getitem__(self, imgId):
        """ Return an element from the set.
        Params:
         imgId: the id of the image
        """
        self.loadIfEmpty()
#        md  = xmipp.MetaData()
#        md.importObjects(self._set._md, xmipp.MDValueEQ(self._label, imgId))
        # TODO: Make a Md query to get the element with this filename
        for item in self._set:
            if item.getId() == imgId:
                return item
            
        raise Exception("getitem for XmippSetOfImages not implemented")
    
    def _getListBlock(self):
        return 'Images@' + self.getFileName()
    
    def _getTiltedBlock(self):
        return 'TiltedPairs@' + self.getFileName()
                    
    @staticmethod
    def _convert(setClass, setOfImgs, filename):
        if isinstance(setOfImgs, setClass):
            return setOfImgs
        
        xmippImgs = setClass(filename)
        xmippImgs.copyInfo(setOfImgs)
        
        for item in setOfImgs:
            xmippImgs.append(item)
        
        # If there are tilt pairs add tilt pairs metadata
        if setOfImgs.hasTiltPairs():
            for iU, iT in setOfImgs.iterTiltPairs():
                xmippImgs.appendPair(iU, iT)
    
        xmippImgs.write()
        
        return xmippImgs
    
    @staticmethod
    def convert(setOfImgs, filename):
        return XmippSetOfImages._convert(XmippSetOfImages, setOfImgs, filename)
    
                
class XmippMicrograph(XmippImage, Micrograph):
    """Xmipp implementation for Micrograph"""
    _label = xmipp.MDL_MICROGRAPH
    
    def __init__(self, filename=None, **args):
        XmippImage.__init__(self, filename, **args)
        Micrograph.__init__(self, filename, **args)
            
    @staticmethod
    def convert(img, ctfFn):
        return XmippImage._convert(XmippMicrograph, img, ctfFn)


class XmippVolume(XmippImage, Volume):
    """Xmipp implementation for Volume"""
    
    def __init__(self, filename=None, **args):
        XmippImage.__init__(self, filename, **args)
        Volume.__init__(self, filename, **args)    
    
    
class XmippSetOfMicrographs(XmippSetOfImages, SetOfMicrographs):
    """Represents a set of Micrographs for Xmipp"""
    _label = xmipp.MDL_MICROGRAPH
    def __init__(self, filename=None, **args):
        SetOfMicrographs.__init__(self, filename, **args)
        XmippSetOfImages.__init__(self, filename, **args)
        
        self._set = XmippSet(XmippMicrograph)

    @staticmethod
    def convert(setOfImgs, filename):
        return XmippSetOfImages._convert(XmippSetOfMicrographs, setOfImgs, filename)
        
    def _getListBlock(self):
        return 'Micrographs@' + self.getFileName()
    

class XmippSetOfParticles(XmippSetOfImages, SetOfParticles):
    """Represents a set of Micrographs for Xmipp"""
    def __init__(self, filename=None, **args):
        SetOfParticles.__init__(self, filename, **args)
        XmippSetOfImages.__init__(self, filename, **args)

    @staticmethod
    def convert(setOfImgs, filename):
        return XmippSetOfImages._convert(XmippSetOfParticles, setOfImgs, filename)
    

class XmippSetOfVolumes(XmippSetOfImages, SetOfVolumes):
    """Represents a set of Volumes for Xmipp"""
    def __init__(self, filename=None, **args):
        SetOfVolumes.__init__(self, filename, **args)
        XmippSetOfImages.__init__(self, filename, **args)
                
        
class XmippCoordinate(Coordinate, XmippMdRow):
    """This class holds the (x,y) position and other information
    associated with a Xmipp coordinate"""
    
    _label_coordId = xmipp.MDL_ITEM_ID
    _label_coordX = xmipp.MDL_XCOOR
    _label_coordY = xmipp.MDL_YCOOR
    
    def __init__(self, **args):
        Coordinate.__init__(self, **args)
        XmippMdRow.__init__(self)
        
    def setId(self, newId):
        #self.id = newId
        self.setValue(self._label_coordId, newId)
        
    def getId(self):
        #return self.id
        return self.getValue(self._label_coordId)
    
    def getPosition(self):
        """Return the position of the coordinate.
        mode: select if the position is the center of the box
          or in the top left corner."""
        return (self.getValue(self._label_coordX), 
                self.getValue(self._label_coordY))

    def setPosition(self, x, y):
        self.setValue(self._label_coordX, x)
        self.setValue(self._label_coordY, y)
    
    def getMicrograph(self):
        """Return the micrograph object to which
        this coordinate is associated"""
        return self._micrograph
    
    def setMicrograph(self, micrograph):
        """Set the micrograph to which this coordinate belongs"""
        self._micrograph = micrograph
    

class XmippCTFModel(CTFModel, XmippMdRow):
    # Mapping to general CTFModel attributes
    ctfParams = {
                 "samplingRate":xmipp.MDL_CTF_SAMPLING_RATE,
                 "voltage":xmipp.MDL_CTF_VOLTAGE,
                 "defocusU":xmipp.MDL_CTF_DEFOCUSU,
                 "defocusV":xmipp.MDL_CTF_DEFOCUSV,
                 "defocusAngle":xmipp.MDL_CTF_DEFOCUS_ANGLE,
                 "sphericalAberration":xmipp.MDL_CTF_CS
                 }

    def __init__(self, filename = None, **args):
        CTFModel.__init__(self, **args)
        XmippMdRow.__init__(self)
        self.getFileName = self.get
        self.setFileName = self.set
        self.setFileName(filename)  
          
        if filename is not None:
            # Use the object value to store the filename

            md = xmipp.MetaData(filename)
            objId = md.firstObject()
            self.getFromMd(md, objId)
            
            for key, label in  self.ctfParams.iteritems():
                mdVal = md.getValue(label, objId)
                getattr(self, key).set(mdVal)
      
    
    def write(self, fn):

        for key, label in  self.ctfParams.iteritems():
            self.setValue(label, getattr(self, key).get())
        
        md = xmipp.MetaData()
        md.setColumnFormat(False)
        objId = md.addObject()
        self.setToMd(md, objId)
        md.write(fn)
        
        self.setFileName(fn)
    
    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths

    @staticmethod
    def convert(ctfModel, filename):
        """ Method to convert from a general ctfModel to XmippCtfModel """

        if isinstance(ctfModel, XmippCTFModel):
            return ctfModel
            
        xmippCTFModel = XmippCTFModel()
        xmippCTFModel.copyInfo(ctfModel)
        
        xmippCTFModel.write(filename)
        
        return xmippCTFModel          
    
    
class XmippSetOfCoordinates(SetOfCoordinates):
    """Implementation of SetOfCoordinates for Xmipp"""
    def __init__(self, filename=None, **args):
        # Use object value to store filename
        # Here filename is the path to a metadata where pos filePaths are linked to micrographs ids
        SetOfCoordinates.__init__(self, value=filename, **args)
        
    def getFileName(self):
        return self.get()       
    
    def getSize(self):
        """ Return the number of coordinates on the set """
        size = 0
        for posFn in self.iterCoordinatesFile():
            mdPos = xmipp.MetaData('particles@%s' % posFn)
            size += mdPos.size()
        return size      
        
    def iterMicrographCoordinates(self, micrograph):
        """ Iterates over the set of coordinates belonging to that micrograph. """
        pathPos = self.getMicrographCoordFile(micrograph.getId())
        
        if pathPos is not None and exists(pathPos):
            mdPos = xmipp.MetaData('particles@' + pathPos)
                            
            for objId in mdPos:
                #x = mdPos.getValue(xmipp.MDL_XCOOR, objId)
                #y = mdPos.getValue(xmipp.MDL_YCOOR, objId)
                #coorId = mdPos.getValue(xmipp.MDL_ITEM_ID, objId)
                coordinate = XmippCoordinate()
                coordinate.getFromMd(mdPos, objId)
                #coordinate.setPosition(x, y)
                coordinate.setMicrograph(micrograph)
                coordinate.setBoxSize(self.boxSize.get())
                #coordinate.setId(coorId)        
                yield coordinate
                
    def iterCoordinates(self):
        """Iterates over the whole set of coordinates.
        If the SetOfMicrographs has tilted pairs, the coordinates
        should have the information related to its paired coordinate."""

        for mic in self.getMicrographs():
            for coord in self.iterMicrographCoordinates(mic):
                yield coord
                
    def iterCoordinatesFile(self):
        """ Iterates over the micrographs_coordinates file
        returning each position file.
        """
        micPosMd = xmipp.MetaData(self.getFileName())
        for objId in micPosMd:
            yield micPosMd.getValue(xmipp.MDL_MICROGRAPH_PARTICLES, objId)        
        
    
    def getFiles(self):
        filePaths = set()

        for posFn in self.iterCoordinatesFile():
            if exists(posFn):         
                filePaths.add(posFn)
        return filePaths
    
    def hasTiltPairs(self):
        return self.getMicrographs().hasTiltPairs()

    def iterTiltPairs(self):
        """Iterate over the tilt pairs if is the case"""     
        mics = self.getMicrographs()
        for micUId, micTId in mics.iterTiltPairs():
            micU = mics[micUId]
            micT = mics[micTId]
            for coordU, coordT in zip(self.iterMicrographCoordinates(micU),
                                      self.iterMicrographCoordinates(micT)):
                yield (coordU.getId(), coordT.getId())
                     
                
    def getMicrographCoordFile(self, micId):
        """ This function will return the pos file corresponding to a micrograph item id"""
        micPosMd = xmipp.MetaData(self.getFileName())
        micRow = findRowById(micPosMd, micId)
        if micRow is not None:    
            return micRow.getValue(xmipp.MDL_MICROGRAPH_PARTICLES)
        return None
    
    @staticmethod
    def convert(setOfCoords, filename):
        if isinstance(setOfCoords, XmippSetOfCoordinates):
            return setOfCoords
        
        print "CONVERTING: ", type(setOfCoords), " TO XmippSetOfCoordinates"
        print "   ", type(setOfCoords) is XmippSetOfCoordinates
        
        xmippCoords = XmippSetOfCoordinates(filename)
        xmippCoords.setMicrographs(setOfCoords.getMicrographs())
        xmippCoords.setBoxSize(setOfCoords.getBoxSize())
        
        extraPath = dirname(filename)
        posMd = xmipp.MetaData()
        # write a position file per micrograph that has coordinates
        for mic in setOfCoords.getMicrographs():
            posFile = join(extraPath, replaceBaseExt(mic.getFileName(), 'pos'))
            mdPosFile = xmipp.MetaData()
            hasCoords = False
            for coord in setOfCoords.iterMicrographCoordinates(mic):
                coordXmipp = XmippCoordinate()
                coordXmipp.copyInfo(coord)
                objId = mdPosFile.addObject()
                coordXmipp.setToMd(mdPosFile, objId)
                hasCoords = True            
            if hasCoords:                           
                mdPosFile.write('particles@%s' % posFile)
                posId = posMd.addObject()
                posMd.setValue(xmipp.MDL_ITEM_ID, long(mic.getId()), posId)
                posMd.setValue(xmipp.MDL_MICROGRAPH_PARTICLES, str(posFile), posId)
            
        # write the micrograph_coordinates file
        posMd.write(filename) 
        
        return xmippCoords    

            
class XmippTiltedPair(XmippMdRow):
    """ Tilted Pairs relations in Xmipp are stored in a MetaData row. """
    def __init__(self, uFn=None, tFn=None, **args):
        XmippMdRow.__init__(self, **args)
        if uFn is not None:
            self.setUntilted(uFn)
        if tFn is not None:
            self.setTilted(tFn)            
        
    def setTilted(self, value):
        self.setValue(xmipp.MDL_MICROGRAPH_TILTED, value)
        
    def getTilted(self):
        return self.getValue(xmipp.MDL_MICROGRAPH_TILTED)
        
    def setUntilted(self, value):
        self.setValue(xmipp.MDL_MICROGRAPH, value)
        
    def getUntilted(self):
        return self.getValue(xmipp.MDL_MICROGRAPH)
        
#    @staticmethod
#    def convert((uFn, tFn)):
#        """ Create a XmippTiltedPair from a general TiltedPair instance. """
#        xmippPair = XmippTiltedPair()
#        xmippPair.setValue(xmipp.MDL_MICROGRAPH, uFn)
#        xmippPair.setValue(xmipp.MDL_MICROGRAPH_TILTED, tFn)
#            
#        return xmippPair
    
    
class XmippImageClassAssignment(ImageClassAssignment, XmippMdRow):
    """ Image-class assignments in Xmipp are stored in a MetaData row. """
    def __init__(self, **args):
        XmippMdRow.__init__(self, **args)
        ImageClassAssignment.__init__(self, **args)
        
    def setImage(self, image):
        """ Set associated image. """
        self.setValue(xmipp.MDL_IMAGE, image.getFileName())
        
    def getImage(self):
        """ Get associated image. """
        xmippImg = XmippImage(self.getValue(xmipp.MDL_IMAGE))
        xmippImg.setValue(xmipp.MDL_CTF_MODEL, self.getValue(xmipp.MDL_CTF_MODEL))
        return xmippImg
    
    
class XmippClass2D(Class2D):
    """ Xmipp classes are stored in blocks of the MetaData. """
    def __init__(self, classNumber, filename, representative=None, **args):
        Class2D.__init__(self, **args)
        self._number = classNumber
        self._filename = filename
        self._representative = representative
        
    def __iter__(self):
        md = xmipp.MetaData('class%06d_images@%s' % 
                            (self._number, self._filename))
        for objId in md:
            imgCA = XmippImageClassAssignment()
            imgCA.getFromMd(md, objId)
            yield imgCA
    
    def getClassRepresentative(self):
        return self._representative
        
        
class XmippClassification2D(Classification2D):
    """ Store results from a 2D classification. """
    def __init__(self, filename=None, classesBlock='classes', **args):
        Classification2D.__init__(self, **args)
        self.getFileName = self.get
        self.setFileName = self.set
        self.setFileName(filename)
        self.classesBlock = String(classesBlock)
        
    def __iter__(self):
        fn = self.getFileName()
        block = self.classesBlock.get()
        md = xmipp.MetaData('%(block)s@%(fn)s' % locals())
        for objId in md:
            ref = md.getValue(xmipp.MDL_REF, objId)
            img = XmippImage(md.getValue(xmipp.MDL_IMAGE, objId))
            yield XmippClass2D(ref, fn, img)
            
    def getClassesMdFileName(self):
        """ Return the filename with block pointing
        to the classes. 
        """
        return "%s@%s" % (self.classesBlock.get(),
                          self.getFileName())
  
