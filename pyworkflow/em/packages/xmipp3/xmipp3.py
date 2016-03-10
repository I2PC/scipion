# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
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

import os
import sys
import platform
from os.path import join
from collections import OrderedDict
from constants import *

import xmipp
import pyworkflow.dataset as ds
from pyworkflow.object import ObjectWrap

from pyworkflow.utils import Environ
from pyworkflow.dataset import COL_RENDER_CHECKBOX, COL_RENDER_TEXT, COL_RENDER_IMAGE,\
    COL_RENDER_VOLUME

from xmipp import MetaData, MetaDataInfo, MDL_IMAGE, MDL_IMAGE1, MDL_IMAGE_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_REF, \
        MDL_SHIFT_X, MDL_SHIFT_Y, MDL_FLIP, MD_APPEND, MDL_MAXCC, MDL_ENABLED, MDL_CTF_MODEL, MDL_SAMPLINGRATE, DT_DOUBLE, \
        MDL_ANGLE_ROT, MDL_SHIFT_Z, Euler_angles2matrix, Image, FileName, getBlocksInMetaDataFile, label2Str

LABEL_TYPES = { 
               xmipp.LABEL_SIZET: long,
               xmipp.LABEL_DOUBLE: float,
               xmipp.LABEL_INT: int,
               xmipp.LABEL_BOOL: bool              
               }

def getEnviron(xmippFirst=True):
    """ Create the needed environment for Xmipp programs. """
    environ = Environ(os.environ)
    pos = Environ.BEGIN if xmippFirst else Environ.END
    environ.update({
            'PATH': join(os.environ['XMIPP_HOME'], 'bin'),
            'LD_LIBRARY_PATH': join(os.environ['XMIPP_HOME'], 'lib')
            }, position=pos)
    return environ
    
def getLabelPythonType(label):
    """ From xmipp label to python variable type """
    labelType = xmipp.labelType(label)
    return LABEL_TYPES.get(labelType, str)

def getXmippPath(*paths):
    '''Return the path the the Xmipp installation folder
    if a subfolder is provided, will be concatenated to the path'''
    if os.environ.has_key('XMIPP_HOME'):
        return os.path.join(os.environ['XMIPP_HOME'], *paths)  
    else:
        raise Exception('XMIPP_HOME environment variable not set')
    
def getMatlabEnviron(*toolPaths):
    """ Return an Environment prepared for launching Matlab
    scripts using the Xmipp binding.
    """
    env = getEnviron()
    env.set('PATH', os.environ['MATLAB_BINDIR'], Environ.BEGIN)
    env.set('LD_LIBRARY_PATH', os.environ['MATLAB_LIBDIR'], Environ.BEGIN)
    for toolpath in toolPaths:
        env.set('MATLABPATH', toolpath, Environ.BEGIN)
    env.set('MATLABPATH', os.path.join(os.environ['XMIPP_HOME'], 'libraries', 'bindings', 'matlab'),
            Environ.BEGIN)
    
    return env
    
    
class XmippProtocol():
    """ This class groups some common functionalities that
    share some Xmipp protocols, like converting steps.
    """
           
    def _insertConvertStep(self, inputName, xmippClass, resultFn):
        """ Insert the convertInputToXmipp if the inputName attribute
        is not an instance of xmippClass.
        It will return the result filename, if the 
        conversion is needed, this will be input resultFn.
        If not, it will be inputAttr.getFileName()
        """
        inputAttr = getattr(self, inputName)
        if not isinstance(inputAttr, xmippClass):
            self._insertFunctionStep('convertInputToXmipp', inputName, xmippClass, resultFn)
            return resultFn
        return inputAttr.getFileName()
         
    def convertInputToXmipp(self, inputName, xmippClass, resultFn):
        """ This step can be used whenever a conversion is needed.
        It will receive the inputName and get this attribute from self,
        invoke the convert function and check the result files if
        conversion was done (otherwise the input was already in Xmipp format).
        """
        inputAttr = getattr(self, inputName)
        inputXmipp = xmippClass.convert(inputAttr, resultFn)
         
        if inputXmipp is not inputAttr:
            self._insertChild(inputName + 'Xmipp', inputXmipp)
            return [resultFn] # validate resultFn was produced if converted
         
    def getConvertedInput(self, inputName):
        """ Retrieve the converted input, it can be the case that
        it is the same as input, when not conversion was done. 
        """
        return getattr(self, inputName + 'Xmipp', getattr(self, inputName))
        

class XmippMdRow():
    """ Support Xmipp class to store label and value pairs 
    corresponding to a Metadata row. 
    """
    def __init__(self):
        self._labelDict = OrderedDict() # Dictionary containing labels and values
        self._objId = None # Set this id when reading from a metadata
        
    def getObjId(self):
        return self._objId
    
    def hasLabel(self, label):
        return self.containsLabel(label)
    
    def containsLabel(self, label):
        # Allow getValue using the label string
        if isinstance(label, basestring):
            label = xmipp.str2Label(label)
        return label in self._labelDict
    
    def removeLabel(self, label):
        if self.hasLabel(label):
            del self._labelDict[label]
    
    def setValue(self, label, value):
        """args: this list should contains tuples with 
        MetaData Label and the desired value"""
        # Allow setValue using the label string
        if isinstance(label, basestring):
            label = xmipp.str2Label(label)
        self._labelDict[label] = value
            
    def getValue(self, label, default=None):
        """ Return the value of the row for a given label. """
        # Allow getValue using the label string
        if isinstance(label, basestring):
            label = xmipp.str2Label(label)
        return self._labelDict.get(label, default)
    
    def getValueAsObject(self, label, default=None):
        """ Same as getValue, but making an Object wrapping. """
        return ObjectWrap(self.getValue(label, default))
    
    def readFromMd(self, md, objId):
        """ Get all row values from a given id of a metadata. """
        self._labelDict.clear()
        self._objId = objId
        
        for label in md.getActiveLabels():
            self._labelDict[label] = md.getValue(label, objId)
            
    def writeToMd(self, md, objId):
        """ Set back row values to a metadata row. """
        for label, value in self._labelDict.iteritems():
            # TODO: Check how to handle correctly unicode type
            # in Xmipp and Scipion
            t = type(value)
            
            if t is unicode:
                value = str(value)
                
            if t is int and xmipp.labelType(label) == xmipp.LABEL_SIZET:
                value = long(value)
                
            try:
                md.setValue(label, value, objId)
            except Exception, ex:
                print >> sys.stderr, "XmippMdRow.writeToMd: Error writting value to metadata."
                print >> sys.stderr, "                     label: %s, value: %s, type(value): %s" % (label2Str(label), value, type(value))
                raise ex
            
    def readFromFile(self, fn):
        md = xmipp.MetaData(fn)
        self.readFromMd(md, md.firstObject())
        
    def copyFromRow(self, other):
        for label, value in other._labelDict.iteritems():
            self.setValue(label, value)
            
    def __str__(self):
        s = '{'
        for k, v in self._labelDict.iteritems():
            s += '  %s = %s\n' % (xmipp.label2Str(k), v)
        return s + '}'
    
    def __iter__(self):
        return self._labelDict.iteritems()
            
    def printDict(self):
        """ Fancy printing of the row, mainly for debugging. """
        print str(self)
    
    
class RowMetaData():
    """ This class is a wrapper for MetaData in row mode.
    Where only one object is used.
    """
    def __init__(self, filename=None):
        self._md = xmipp.MetaData()
        self._md.setColumnFormat(False)
        self._id = self._md.addObject()
        
        if filename:
            self.read(filename)
        
    def setValue(self, label, value):
        self._md.setValue(label, value, self._id)
        
    def getValue(self, label):
        return self._md.getValue(label, self._id)
        
    def write(self, filename, mode=xmipp.MD_APPEND):
        self._md.write(filename, mode)
        
    def read(self, filename):
        self._md.read(filename)
        self._md.setColumnFormat(False)
        self._id = self._md.firstObject()
        
    def __str__(self):
        return str(self._md)
    
        
def findRow(md, label, value):
    """ Query the metadata for a row with label=value.
    Params:
        md: metadata to query.
        label: label to check value
        value: value for equal condition
    Returns:
        XmippMdRow object of the row found.
        None if no row is found with label=value
    """
    mdQuery = xmipp.MetaData() # store result
    mdQuery.importObjects(md, xmipp.MDValueEQ(label, value))
    n = mdQuery.size()
    
    if n == 0:
        row = None
    elif n == 1:
        row = XmippMdRow()
        row.readFromMd(mdQuery, mdQuery.firstObject())
    else:
        raise Exception("findRow: more than one row found matching the query %s = %s" % (xmipp.label2Str(label), value))
    
    return row

def findRowById(md, value):
    """ Same as findRow, but using MDL_ITEM_ID for label. """
    return findRow(md, xmipp.MDL_ITEM_ID, long(value))
  
  
class XmippSet():
    """ Support class to store sets in Xmipp base on a MetaData. """
    def __init__(self, itemClass):
        """ Create new set, base on a Metadata.
        itemClass: Class that represent the items.
        A method .getFileName should be available to store the md.
        Items contained in XmippSet are supposed to inherit from XmippMdRow.
        """
        self._itemClass = itemClass 
        #self._fileName = fileName      
        self._md = xmipp.MetaData()
        
        
    def __iter__(self):
        """Iterate over the set of images in the MetaData"""
        #self._md.read(self._fileName)
        
        for objId in self._md:  
            item = self._itemClass()
            item.readFromMd(self._md, objId)  
            #m = Image(md.getValue(xmipp.MDL_IMAGE, objId))
            #if self.hasCTF():
            #    m.ctfModel = XmippCTFModel(md.getValue(xmipp.MDL_CTF_MODEL, objId)) 
            yield item
        
#        
#    def setFileName(self, filename):
#        self._fileName = filename
        
    def setMd(self, md):
        self._md = md
        
    def write(self, filename, mode):
        self._md.write(filename, mode)
        
    def append(self, item):
        """Add a new item to the set, the item can be of a base class (EM)
        and will be try to convert it to the respective _itemClass."""
        objId = self._md.addObject()
        # Convert to xmipp micrograph if necessary
        if isinstance(item, self._itemClass):
            itemXmipp = item
        else:
            itemXmipp = self._itemClass.convert(item)
        itemXmipp.writeToMd(self._md, objId)
        
    def sort(self, label):
        self._md.sort(label)
        
    def read(self, filename):
        self._md.read(filename)
        
    def isEmpty(self):
        return self._md.isEmpty()
    
    def getItemClass(self):
        """ Return the Class of the items in the Set. """
        return self._itemClass
    
    def getSize(self):
        return self._md.size()
        
    def convert(self, xmippSetClass, filename):
        """ Convert from a generic set to a xmippSetClass.
        In particular a filename is requiered to store the result MetaData.
        It is also asummed that this class have a .copyInfo method.
        """
        if isinstance(self, xmippSetClass):
            return self
        
        setOut = xmippSetClass(filename)
        setOut.copyInfo(self)
        
        for item in self:
            setOut.append(item)
        setOut.write()
        
        return setOut   
    
    
class XmippDataSet(ds.DataSet):
    """ Provide a DataSet implementation based on Xmipp xmd file.
    The tables of the dataset will be the blocks in the metadata.
    Each block is a table on the dataset and is read as an Xmipp metadata. 
    """
    def __init__(self, filename):
        self._filename = filename
        blocks = xmipp.getBlocksInMetaDataFile(filename)
        
        if not blocks: # If there are no block, add an empty one
            blocks = ['']
        ds.DataSet.__init__(self, blocks)
        
    def _loadTable(self, tableName):
        if tableName:
            mdFn = tableName + "@" + self._filename
        else:
            mdFn = self._filename
        md = xmipp.MetaData(mdFn)
        return self._convertMdToTable(md)
        
    def _getLabelRenderType(self, label, md):
        """ Return the way to render each label. """
        if xmipp.labelIsImage(label):
            value = md.getValue(label, md.firstObject())
            if value.endswith('.vol'):
                return COL_RENDER_VOLUME
            return COL_RENDER_IMAGE
        
        labelType = xmipp.labelType(label)
        if labelType == xmipp.LABEL_BOOL:
            return COL_RENDER_CHECKBOX
        
        return COL_RENDER_TEXT
        
        
    def _convertLabelToColumn(self, label, md):
        """ From an Xmipp label, create the corresponding column. """
        return ds.Column(xmipp.label2Str(label), 
                         getLabelPythonType(label),
                         renderType=self._getLabelRenderType(label, md))
        
    def _convertMdToTable(self, md):
        """ Convert a metatada into a table. """
           
        labels = md.getActiveLabels()
        hasTransformation = self._hasTransformation(labels)  
        columns = [self._convertLabelToColumn(l, md) for l in labels]        
        #NAPA de LUXE (xmipp deberia saber a que campo va asignado el transformation matrix)             
        if hasTransformation:
            columns.append(ds.Column("image_transformationMatrix", str))
            
        labelsStr = [col.getName() for col in columns]        
        
        table = ds.Table(*columns)
        
        for objId in md:
            values = [md.getValue(l, objId) for l in labels]
            if hasTransformation:
                values.append(self._getTransformation(md, objId))
            d = dict(zip(labelsStr, values))
            table.addRow(objId, **d)
            
        return table
    
    def _convertTableToMd(self, table):
        colLabels = [(col.getName(), xmipp.str2Label(col.getName())) 
                     for col in table.iterColumns()]
        md = xmipp.MetaData()
        
        for row in table.iterRows():
            objId = md.addObject()
            for col, label in colLabels:
                if col != 'id':
                    value = getattr(row, col)
                    md.setValue(label, value, objId)
                
        return md
        
    def writeTable(self, tableName, table):
        """ Write changes made to a table. """
        md = self._convertTableToMd(table)
        md.write("%s@%s" % (tableName, self._filename), xmipp.MD_APPEND)

        
    def _hasTransformation(self, labels):
        for l in [xmipp.MDL_SHIFT_X, xmipp.MDL_SHIFT_Y, xmipp.MDL_SHIFT_Z, xmipp.MDL_ANGLE_ROT, xmipp.MDL_ANGLE_TILT, xmipp.MDL_ANGLE_PSI]:
            if l in labels:
                return True
        return False
        
    def _getTransformation(self, md, objId):
        rot  = md.getValue(xmipp.MDL_ANGLE_ROT, objId) or 0.
        tilt = md.getValue(xmipp.MDL_ANGLE_TILT, objId) or 0.
        psi  = md.getValue(xmipp.MDL_ANGLE_PSI, objId) or 0.

        tMatrix = xmipp.Euler_angles2matrix(rot, tilt, psi)
        x = md.getValue(xmipp.MDL_SHIFT_X, objId) or 0.
        y = md.getValue(xmipp.MDL_SHIFT_Y, objId) or 0.
        z = md.getValue(xmipp.MDL_SHIFT_Z, objId) or 0.

        matrix = [tMatrix[0][0], tMatrix[0][1], tMatrix[0][2], x,
                  tMatrix[1][0], tMatrix[1][1], tMatrix[1][2], y,
                  tMatrix[2][0], tMatrix[2][1], tMatrix[2][2], z]

        return matrix
        
        
class ProjMatcher():
    """ Base class for protocols that use a projection """
    
    def projMatchStep(self, volume, angularSampling, symmetryGroup, images, fnAngles, Xdim):
        from pyworkflow.utils.path import cleanPath
        # Generate gallery of projections        
        fnGallery = self._getExtraPath('gallery.stk')
        if volume.endswith('.mrc'):
            volume+=":mrc"
        
        self.runJob("xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                   % (volume, fnGallery, angularSampling, symmetryGroup, images))
    
        # Assign angles
        self.runJob("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 1000 --search5d_shift %s --search5d_step  %s --append"\
                   % (images, fnAngles, fnGallery, str(Xdim/2), str(int(Xdim/10)), str(int(Xdim/25))))
        
        cleanPath(self._getExtraPath('gallery_sampling.xmd'))
        cleanPath(self._getExtraPath('gallery_angles.doc'))
        cleanPath(self._getExtraPath('gallery.doc'))
    
        # Write angles in the original file and sort
        MD=MetaData(fnAngles)
        for id in MD:
            galleryReference = MD.getValue(xmipp.MDL_REF,id)
            MD.setValue(xmipp.MDL_IMAGE_REF, "%05d@%s" % (galleryReference+1,fnGallery), id)
        MD.write(fnAngles)
        
    def produceAlignedImagesStep(self, volumeIsCTFCorrected, fn, images):
        
        from numpy import array, dot
        fnOut = 'classes_aligned@' + fn
        MDin = MetaData(images)
        MDout = MetaData()
        n = 1
        hasCTF = MDin.containsLabel(xmipp.MDL_CTF_MODEL)
        for i in MDin:
            fnImg = MDin.getValue(xmipp.MDL_IMAGE,i)
            fnImgRef = MDin.getValue(xmipp.MDL_IMAGE_REF,i)
            maxCC = MDin.getValue(xmipp.MDL_MAXCC,i)
            rot =  MDin.getValue(xmipp.MDL_ANGLE_ROT,i)
            tilt = MDin.getValue(xmipp.MDL_ANGLE_TILT,i)
            psi =-1.*MDin.getValue(xmipp.MDL_ANGLE_PSI,i)
            flip = MDin.getValue(xmipp.MDL_FLIP,i)
            if flip:
                psi = -psi
            eulerMatrix = Euler_angles2matrix(0.,0.,psi)
            x = MDin.getValue(xmipp.MDL_SHIFT_X,i)
            y = MDin.getValue(xmipp.MDL_SHIFT_Y,i)
            shift = array([x, y, 0])
            shiftOut = dot(eulerMatrix, shift)
            [x,y,z]= shiftOut
            if flip:
                x = -x
            id = MDout.addObject()
            MDout.setValue(xmipp.MDL_IMAGE, fnImg, id)
            MDout.setValue(xmipp.MDL_IMAGE_REF, fnImgRef, id)
            MDout.setValue(xmipp.MDL_IMAGE1, "%05d@%s"%(n, self._getExtraPath("diff.stk")), id)
            if hasCTF:
                fnCTF = MDin.getValue(xmipp.MDL_CTF_MODEL,i)
                MDout.setValue(xmipp.MDL_CTF_MODEL,fnCTF,id)
            MDout.setValue(xmipp.MDL_MAXCC, maxCC, id)
            MDout.setValue(xmipp.MDL_ANGLE_ROT, rot, id)
            MDout.setValue(xmipp.MDL_ANGLE_TILT, tilt, id)
            MDout.setValue(xmipp.MDL_ANGLE_PSI, psi, id)
            MDout.setValue(xmipp.MDL_SHIFT_X, x,id)
            MDout.setValue(xmipp.MDL_SHIFT_Y, y,id)
            MDout.setValue(xmipp.MDL_FLIP,flip,id)
            MDout.setValue(xmipp.MDL_ENABLED,1,id)
            n+=1
        MDout.write(fnOut,xmipp.MD_APPEND)
        
        # Actually create the differences
        img = Image()
        imgRef = Image()
        if hasCTF and volumeIsCTFCorrected:
            Ts = MDin.getValue(xmipp.MDL_SAMPLINGRATE, MDin.firstObject())
    
        for i in MDout:
            img.readApplyGeo(MDout,i)
            imgRef.read(MDout.getValue(xmipp.MDL_IMAGE_REF,i))
            if hasCTF and volumeIsCTFCorrected:
                fnCTF = MDout.getValue(xmipp.MDL_CTF_MODEL,i)
                imgRef.applyCTF(fnCTF,Ts)
                img.convert2DataType(DT_DOUBLE)
            imgDiff = img-imgRef
            imgDiff.write(MDout.getValue(xmipp.MDL_IMAGE1,i))

class HelicalFinder():
    """ Base class for protocols that find helical symmetry """

    def getSymmetry(self,dihedral):
        if dihedral:
            return "helicalDihedral"
        else:
            return "helical"

    def runCoarseSearch(self,fnVol,dihedral,heightFraction,z0,zF,zStep,rot0,rotF,rotStep,Nthr,fnOut,cylinderInnerRadius,cylinderOuterRadius,height,Ts):
        args="-i %s --sym %s --heightFraction %f -z %f %f %f --rotHelical %f %f %f --thr %d -o %s --sampling %f"%(fnVol,self.getSymmetry(dihedral),heightFraction,
                                                                                z0,zF,zStep,rot0,rotF,rotStep,Nthr,fnOut,Ts)
        if cylinderOuterRadius>0 and cylinderInnerRadius<0:
            args+=" --mask cylinder %d %d"%(-cylinderOuterRadius,-height)
        elif cylinderOuterRadius>0 and cylinderInnerRadius>0:
            args+=" --mask tube %d %d %d"%(-cylinderInnerRadius,-cylinderOuterRadius,-height)
        self.runJob('xmipp_volume_find_symmetry',args)

    def runFineSearch(self, fnVol, dihedral, fnCoarse, fnFine, heightFraction, z0, zF, rot0, rotF, cylinderInnerRadius, cylinderOuterRadius, height, Ts):
        md=MetaData(fnCoarse)
        objId=md.firstObject()
        rotInit=md.getValue(MDL_ANGLE_ROT,objId)
        zInit=md.getValue(MDL_SHIFT_Z,objId)
        args="-i %s --sym %s --heightFraction %f --localHelical %f %f -o %s -z %f %f 1 --rotHelical %f %f 1 --sampling %f"%(fnVol,self.getSymmetry(dihedral),heightFraction,
                                                                                           zInit,rotInit,fnFine,z0,zF,rot0,rotF,Ts)
        if cylinderOuterRadius>0 and cylinderInnerRadius<0:
            args+=" --mask cylinder %d %d"%(-cylinderOuterRadius,-height)
        elif cylinderOuterRadius>0 and cylinderInnerRadius>0:
            args+=" --mask tube %d %d %d"%(-cylinderInnerRadius,-cylinderOuterRadius,-height)
        self.runJob('xmipp_volume_find_symmetry',args)

    def runSymmetrize(self, fnVol, dihedral, fnParams, fnOut, heightFraction, cylinderInnerRadius, cylinderOuterRadius, height, Ts):
        md=MetaData(fnParams)
        objId=md.firstObject()
        rot0=md.getValue(MDL_ANGLE_ROT,objId)
        z0=md.getValue(MDL_SHIFT_Z,objId)
        args="-i %s --sym %s --helixParams %f %f --heightFraction %f -o %s --sampling %f"%(fnVol,self.getSymmetry(dihedral),z0,rot0,heightFraction,fnOut,Ts)
        self.runJob('xmipp_transform_symmetrize',args)
        doMask=False
        if cylinderOuterRadius>0 and cylinderInnerRadius<0:
            args="-i %s --mask cylinder %d %d"%(fnVol,-cylinderOuterRadius,-height)
            doMask=True
        elif cylinderOuterRadius>0 and cylinderInnerRadius>0:
            args="-i %s --mask tube %d %d %d"%(fnVol,-cylinderInnerRadius,-cylinderOuterRadius,-height)
            doMask=True
        if doMask:
            self.runJob('xmipp_transform_mask',args)
            

# ---------------- Legacy code from 'protlib_xmipp.py' ----------------------

def xmippExists(path):
    return xmipp.FileName(path).exists()
            
            
class XmippScript():
    ''' This class will serve as wrapper around the XmippProgram class
    to have same facilities from Python scripts'''
    def __init__(self, runWithoutArgs=False):
        self._prog = xmipp.Program(runWithoutArgs)
        
    def defineParams(self):
        ''' This function should be overwrited by subclasses for 
        define its own parameters'''
        pass
    
    def readParams(self):
        ''' This function should be overwrited by subclasses for 
        and take desired params from command line'''
        pass
    
    def checkParam(self, param):
        return self._prog.checkParam(param)
    
    def getParam(self, param, index=0):
        return self._prog.getParam(param, index)
    
    def getIntParam(self, param, index=0):
        return int(self._prog.getParam(param, index))
    
    def getDoubleParam(self, param, index=0):
        return float(self._prog.getParam(param, index))
    
    def getListParam(self, param):
        return self._prog.getListParam(param)
    
    def addUsageLine(self, line, verbatim=False):
        self._prog.addUsageLine(line, verbatim)

    def addExampleLine(self, line, verbatim=True):
        self._prog.addExampleLine(line, verbatim)
        
    def addParamsLine(self, line):
        self._prog.addParamsLine(line)
    
    def run(self):
        ''' This function should be overwrited by subclasses and
        it the main body of the script'''   
        pass
     
    def tryRun(self):
        ''' This function should be overwrited by subclasses and
        it the main body of the script'''
        try:
            self.defineParams()
            doRun = self._prog.read(sys.argv)
            if doRun:
                self.readParams()
                self.run()
        except Exception:
            import traceback
            traceback.print_exc(file=sys.stderr)
            
        
class ScriptIJBase(XmippScript):
    def __init__(self, name):
        XmippScript.__init__(self)
        self.name = name
    
    def defineOtherParams(self):
        pass
    
    def readOtherParams(self):
        pass
    
    def defineParams(self):
        self.addParamsLine('  --input <...>            : Input files to show');
        self.addParamsLine('         alias -i;');
        self.addParamsLine('  [--memory <mem="2g">]    : Memory amount for JVM');
        self.addParamsLine('         alias -m;');
        self.defineOtherParams()
    
    def readInputFiles(self):
        self.inputFiles = self.getListParam('-i')
        
    def readParams(self):
        self.readInputFiles()
        self.memory = self.getParam('--memory')
        files = []
        missingFiles = []
        for f in self.inputFiles:
            if xmippExists(f):
                files.append('"%s"' % f) # Escape with " for filenames containing spaces
            else:
                missingFiles.append(f)
        self.inputFiles = files
        if len(missingFiles):
            print "Missing files: \n %s" % '  \n'.join(missingFiles)
        
        self.args = "-i %s" % ' '.join(self.inputFiles)
        self.readOtherParams()
 
 
class ScriptPluginIJ(ScriptIJBase):
    def __init__(self, macro):
        ScriptIJBase.__init__(self, macro)
                  
    def run(self):
        runImageJPlugin(self.memory, self.name, self.args)
     
     
class ScriptAppIJ(ScriptIJBase):
    def __init__(self, name):
        ScriptIJBase.__init__(self, name)
                  
    def run(self):
        if len(self.inputFiles) > 0:
            runJavaIJapp(self.memory, self.name, self.args, batchMode=False)
        else:
            print "No input files. Exiting..."
            
    
def getImageJPluginCmd(memory, macro, args, batchMode=False):
    if len(memory) == 0:
        memory = "1g"
        print "No memory size provided. Using default: " + memory
    imagej_home = getXmippPath('external', 'imagej')
    plugins_dir = os.path.join(imagej_home, "plugins")
    macro = os.path.join(imagej_home, "macros", macro)
    imagej_jar = os.path.join(imagej_home, "ij.jar")
    cmd = """ java -Xmx%s -Dplugins.dir=%s -jar %s -macro %s "%s" """ % (memory, plugins_dir, imagej_jar, macro, args)
    if batchMode:
        cmd += " &"
    return cmd


def runImageJPlugin(memory, macro, args, batchMode=False):
    os.system(getImageJPluginCmd(memory, macro, args, batchMode))


def getArchitecture():
    arch = platform.architecture()[0]
    for a in ['32', '64']:
        if a in arch:
            return a
    return 'NO_ARCH' 
    
def getJavaIJappCmd(memory, appName, args, batchMode=False):
    '''Launch an Java application based on ImageJ '''
    if len(memory) == 0:
        memory = "2g"
        print "No memory size provided. Using default: " + memory
    imagej_home = getXmippPath("external", "imagej")
    lib = getXmippPath("lib")
    javaLib = getXmippPath('java', 'lib')
    plugins_dir = os.path.join(imagej_home, "plugins")
    arch = getArchitecture()
    cmd = "java -Xmx%(memory)s -d%(arch)s -Djava.library.path=%(lib)s -Dplugins.dir=%(plugins_dir)s -cp %(imagej_home)s/*:%(javaLib)s/* %(appName)s %(args)s" % locals()
    if batchMode:
        cmd += " &"
    return cmd
    
def runJavaIJapp(memory, appName, args, batchMode=True):
    cmd = getJavaIJappCmd(memory, appName, args, batchMode)
    print cmd
    os.system(cmd)
    
def runJavaJar(memory, jarName, args, batchMode=True):
    jarPath = getXmippPath(jarName)
    runJavaIJapp(memory, '-jar %s' % jarPath, args, batchMode)



class ScriptShowJ(ScriptAppIJ):
    def __init__(self, viewer='xmipp.viewer.Viewer'):
        ScriptAppIJ.__init__(self, viewer)
        
    def defineOtherParams(self):
        self.addParamsLine('  [--mode <mode_value=image>]           : List of params ')
        self.addParamsLine('     where <mode_value> image gallery metadata rotspectra')
        self.addParamsLine('         alias -d;')
        self.addParamsLine('  [--poll]                            : Keeps checking for changes on input files  (for image mode only!)')
        self.addParamsLine('         alias -p;')
        self.addParamsLine('  [--render <...>]    : Specifies image columns to render (for metadata mode only)')
        self.addParamsLine('                          : by default the first one that can be visualized is rendered')
        self.addParamsLine('  [--visible <...>]    : Specifies visible labels')
        self.addParamsLine('  [--order <...>]    : Specifies labels order')
        self.addParamsLine('  [--labels <...>]    : Specifies labels to display')
        self.addParamsLine('  [--sortby <...>]    : Specifies label to sort by. asc or desc mode can be added')
        
        self.addParamsLine('         alias -e;')
        self.addParamsLine('  [--rows <rows>]                            : number of rows in table')
        self.addParamsLine('         alias -r;')
        self.addParamsLine('  [--columns <columns>]                            : number of columns in table')
        self.addParamsLine('         alias -c;')
        self.addParamsLine('  [--zoom <zoom>]                            : zoom for images')
        self.addParamsLine('         alias -z;')
        self.addParamsLine('  [--view <axis="z">]                        : Viewer position (for volumes only)')
        self.addParamsLine('     where <axis> z y x z_pos y_pos x_pos')
        self.addParamsLine('  [--dont_apply_geo]                        : Does not read geometrical information(for metadata only)')
        self.addParamsLine('  [--dont_wrap]                             : Does not wrap (for metadata only)')
        self.addParamsLine('  [--debug] : debug')
        self.addParamsLine('  [--mask_toolbar] : Open mask toolbar (only valid for images)')
        self.addParamsLine('  [--label_alias <alias_string>]  : Activate some metadata label alias, for example')
        self.addParamsLine('                                  : anglePsi=aPsi;shiftX=sX;shiftY:sY')
        self.addParamsLine('  [--label_relion]                : Activates the mapping to Relion labels')
        self.addParamsLine('  [--label_bsoft]                 : Activates the mapping to Bsoft labels')
        
    def readOtherParams(self):
        #FIXME: params seems to be they cannot be passed directly to java
        params = ['--mode', '--rows', '--columns', '--zoom', '--view', '--sortby']
        for p in params:
            if self.checkParam(p):
                self.args += " %s %s" % (p, self.getParam(p))
        params = [ '--render', '--visible', '--order', '--labels']
        pvalues = ''
        for p in params:
            if self.checkParam(p):
                for pvalue in self.getListParam(p):
                    pvalues = '%s %s'%(pvalues, pvalue)
                    
                self.args += " %s %s" % (p, pvalues)
        params = ['--poll', '--debug', '--dont_apply_geo', '--dont_wrap', '--mask_toolbar']
        for p in params:
            if self.checkParam(p):
                self.args += " %s" % p
                
        # Set environment var for extra label alias
        if self.checkParam('--label_alias'):
            os.environ['XMIPP_EXTRA_ALIASES'] = self.getParam('--label_alias')
        elif self.checkParam('--label_bsoft'):
            from protlib_import import bsoftLabelString
            os.environ['XMIPP_EXTRA_ALIASES'] = bsoftLabelString()
        elif self.checkParam('--label_relion') or self.getParam('-i').endswith('.star'):
            from protlib_import import relionLabelString
            os.environ['XMIPP_EXTRA_ALIASES'] = relionLabelString()
            

def createMetaDataFromPattern(pattern, isStack=False, label="image"):
    ''' Create a metadata from files matching pattern'''
    import glob
    files = glob.glob(pattern)
    files.sort()

    label = xmipp.str2Label(label) #Check for label value
    
    mD = xmipp.MetaData()
    inFile = xmipp.FileName()
    
    nSize = 1
    for file in files:
        fileAux=file
        if isStack:
            if file.endswith(".mrc"):
                fileAux=file+":mrcs"
            x, x, x, nSize = xmipp.getImageSize(fileAux)
        if nSize != 1:
            counter = 1
            for jj in range(nSize):
                inFile.compose(counter, fileAux)
                objId = mD.addObject()
                mD.setValue(label, inFile, objId)
                mD.setValue(xmipp.MDL_ENABLED, 1, objId)
                counter += 1
        else:
            objId = mD.addObject()
            mD.setValue(label, fileAux, objId)
            mD.setValue(xmipp.MDL_ENABLED, 1, objId)
    return mD            
