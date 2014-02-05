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
"""
This sub-package will contains Xmipp3.0 specific protocols
"""

import os
from collections import OrderedDict

import xmipp
from constants import *
import pyworkflow.dataset as ds

from xmipp import MetaData, MetaDataInfo, MDL_IMAGE, MDL_IMAGE1, MDL_IMAGE_REF, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_REF, \
        MDL_SHIFT_X, MDL_SHIFT_Y, MDL_FLIP, MD_APPEND, MDL_MAXCC, MDL_ENABLED, MDL_CTF_MODEL, MDL_SAMPLINGRATE, DT_DOUBLE, \
        Euler_angles2matrix, Image, FileName, getBlocksInMetaDataFile
from protlib_utils import runJob
from protlib_filesystem import deleteFile, findAcquisitionInfo
import os

LABEL_TYPES = { 
               xmipp.LABEL_SIZET: long,
               xmipp.LABEL_DOUBLE: float,
               xmipp.LABEL_INT: int,
               xmipp.LABEL_BOOL: bool              
               }

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
    
    
class XmippProtocol():
    """ This class groups some common functionalities that
    share some Xmipp protocols, like converting steps.
    """
           
    def _insertConvertStep(self, inputName, xmippClass, resultFn):
        """ Insert the convertInputToXmipp if the inputName attribute
        is not an instance of xmippClass.
        It will return the result filename, if the 
        convertion is needed, this will be input resultFn.
        If not, it will be inputAttr.getFileName()
        """
        inputAttr = getattr(self, inputName)
        print "inputAttr.getClassName()", inputAttr.getClassName()
        print "xmippClass", xmippClass
        if not isinstance(inputAttr, xmippClass):
            self._insertFunctionStep('convertInputToXmipp', inputName, xmippClass, resultFn)
            return resultFn
        return inputAttr.getFileName()
         
    def convertInputToXmipp(self, inputName, xmippClass, resultFn):
        """ This step can be used whenever a convertion is needed.
        It will receive the inputName and get this attribute from self,
        invoke the convert function and check the result files if
        convertion was done (otherwise the input was already in Xmipp format).
        """
        inputAttr = getattr(self, inputName)
        inputXmipp = xmippClass.convert(inputAttr, resultFn)
         
        if inputXmipp is not inputAttr:
            print "======== CONVERTIN........."
            self._insertChild(inputName + 'Xmipp', inputXmipp)
            return [resultFn] # validate resultFn was produced if converted
         
    def getConvertedInput(self, inputName):
        """ Retrieve the converted input, it can be the case that
        it is the same as input, when not convertion was done. 
        """
        return getattr(self, inputName + 'Xmipp', getattr(self, inputName))
        

class XmippMdRow():
    """ Support Xmipp class to store label and value pairs 
    corresponding to a Metadata row. It can be used as base
    for classes that maps to a MetaData row like XmippImage, XmippMicrograph..etc. 
    """
    def __init__(self):
        self._labelDict = OrderedDict() # Dictionary containing labels and values
    
    def hasLabel(self, label):
        return label in self._labelDict
    
    def setValue(self, label, value):
        """args: this list should contains tuples with 
        MetaData Label and the desired value"""
        self._labelDict[label] = value
            
    def getValue(self, label):
        return self._labelDict[label]
    
    def readFromMd(self, md, objId):
        """ Get all row values from a given id of a metadata. """
        self._labelDict.clear()
        for label in md.getActiveLabels():
            self._labelDict[label] = md.getValue(label, objId)
            
    def writeToMd(self, md, objId):
        """ Set back row values to a metadata row. """
        for label, value in self._labelDict.iteritems():
            # TODO: Check how to handle correctly unicode type
            # in Xmipp and Scipion
            if type(value) is unicode:
                value = str(value)
            md.setValue(label, value, objId)
            
    def readFromFile(self, fn):
        md = xmipp.MetaData(fn)
        self.readFromMd(md, md.firstObject())
        
    def copyFromRow(self, other):
        for label, value in other._labelDict.iteritems():
            self.setValue(label, value)
            
    def __str__(self):
        s = '{'
        for k, v in self._labelDict.iteritems():
            s += '%s: %s, ' % (xmipp.label2Str(k), v)
        return s + '}'
    
    
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
        Items contained in XmippSet are suposed to inherit from XmippMdRow.
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
    """ this class will implement a dataset base on
    a metadata file, which can contains several blocks.
    Each block is a table on the dataset and is read
    as a Xmipp metadata. 
    """
    def __init__(self, filename):
        self._filename = filename
        blocks = xmipp.getBlocksInMetaDataFile(filename)
        if len(blocks) == 0: # If there are no block, at an empty one
            blocks = ['']
        ds.DataSet.__init__(self, blocks)
        
    def _loadTable(self, tableName):
        if len(tableName):
            mdFn = tableName + "@" + self._filename
        else:
            mdFn = self._filename
        md = xmipp.MetaData(mdFn)
        return self._convertMdToTable(md)
        
    def _convertLabelToColumn(self, label):
        """ From an Xmipp label, create the corresponding column. """
        return ds.Column(xmipp.label2Str(label), 
                         getLabelPythonType(label))
        
    def _convertMdToTable(self, md):
        """ Convert a metatada into a table. """
        labels = md.getActiveLabels()
        hasTransformation = self._hasTransformation(labels)  
        columns = [self._convertLabelToColumn(l) for l in labels]        
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
        print "dondecarajo","%s@%s" % (tableName, self._filename)
        md.write("%s@%s" % (tableName, self._filename), xmipp.MD_APPEND)

        
    def _hasTransformation(self, labels):
        for l in [xmipp.MDL_SHIFT_X, xmipp.MDL_SHIFT_Y, xmipp.MDL_SHIFT_Z, xmipp.MDL_ANGLE_ROT, xmipp.MDL_ANGLE_TILT, xmipp.MDL_ANGLE_PSI]:
            if l in labels:
                return True
        return False
        
    def _getTransformation(self, md, objId):
        rot  = md.getValue(xmipp.MDL_ANGLE_ROT,objId)
        tilt = md.getValue(xmipp.MDL_ANGLE_TILT,objId)
        psi  = md.getValue(xmipp.MDL_ANGLE_PSI,objId)
        if rot is  None:
            rot = 0
        if tilt is  None:
            tilt = 0
        if psi is  None:
            psi = 0

        tMatrix=xmipp.Euler_angles2matrix(rot,tilt,psi)
        x = md.getValue(xmipp.MDL_SHIFT_X, objId)
        y = md.getValue(xmipp.MDL_SHIFT_Y, objId)
        z = md.getValue(xmipp.MDL_SHIFT_Z, objId)

        matrix = [tMatrix[0][0], tMatrix[0][1], tMatrix[0][2], x if x!=None else 0,
                tMatrix[1][0], tMatrix[1][1], tMatrix[1][2], y if y!=None else 0,
                tMatrix[2][0], tMatrix[2][1], tMatrix[2][2], z if z!=None else 0]
        
        print matrix

        return matrix
        
    def getTypeOfColumn(self, label):
        if (label == "id"):
            return "id"
        elif (label!='image_transformationMatrix' and xmipp.labelIsImage(str(label))):
            return "image"
        elif (label == "enabled"):
            return "checkbox"
        else:
            return "text" 
        
        
class ProjMatcher():
    """ Base class for protocols that use a projection """
    
    def projMatchStep(self, volume):
        from pyworkflow.utils.path import cleanPath
        
#        Xdim=MetaDataInfo(Images)[0]
    
        # Generate gallery of projections        
        fnGallery = self._getExtraPath('gallery.stk')
#        images = "classes@%s" % self.fn
        
        runJob(None,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                   % (volume, fnGallery, self.angularSampling.get(), self.symmetryGroup.get(), self.images), self.numberOfMpi.get())
    
        # Assign angles
        runJob(None,"xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift 1000 --search5d_shift %s --search5d_step  %s --append"\
                   % (self.images, self.fnAngles, fnGallery, str(self.Xdim/2), str(int(self.Xdim/10)), str(int(self.Xdim/25))), self.numberOfMpi.get())
        
        cleanPath(self._getExtraPath('gallery_sampling.xmd'))
        cleanPath(self._getExtraPath('gallery_angles.doc'))
        cleanPath(self._getExtraPath('gallery.doc'))
        
    
        # Write angles in the original file and sort
        MD=MetaData(self.fnAngles)
        for id in MD:
            galleryReference=MD.getValue(MDL_REF,id)
            MD.setValue(MDL_IMAGE_REF,"%05d@%s"%(galleryReference+1,fnGallery),id)
        MD.write(self.fnAngles)
        
    def produceAlignedImagesStep(self, volumeIsCTFCorrected):
        
        from numpy import array, dot
        fnOut = 'classes_aligned@'+self.fn
        MDin=MetaData(self.images)
        MDout=MetaData()
        n=1
        hasCTF=MDin.containsLabel(MDL_CTF_MODEL)
        for i in MDin:
            fnImg=MDin.getValue(MDL_IMAGE,i)
            fnImgRef=MDin.getValue(MDL_IMAGE_REF,i)
            maxCC=MDin.getValue(MDL_MAXCC,i)
            rot =  MDin.getValue(MDL_ANGLE_ROT,i)
            tilt = MDin.getValue(MDL_ANGLE_TILT,i)
            psi =-1.*MDin.getValue(MDL_ANGLE_PSI,i)
            flip = MDin.getValue(MDL_FLIP,i)
            if(flip):
                psi =-psi
            eulerMatrix = Euler_angles2matrix(0.,0.,psi)
            x = MDin.getValue(MDL_SHIFT_X,i)
            y = MDin.getValue(MDL_SHIFT_Y,i)
            shift = array([x, y, 0])
            shiftOut = dot(eulerMatrix, shift)
            [x,y,z]= shiftOut
            if flip:
                x = -x
            id=MDout.addObject()
            MDout.setValue(MDL_IMAGE, fnImg, id)
            MDout.setValue(MDL_IMAGE_REF, fnImgRef, id)
            MDout.setValue(MDL_IMAGE1, "%05d@%s"%(n,self._getExtraPath("diff.stk")), id)
            if hasCTF:
                fnCTF=MDin.getValue(MDL_CTF_MODEL,i)
                MDout.setValue(MDL_CTF_MODEL,fnCTF,id)
            MDout.setValue(MDL_MAXCC, maxCC, id)
            MDout.setValue(MDL_ANGLE_ROT, rot, id)
            MDout.setValue(MDL_ANGLE_TILT, tilt, id)
            MDout.setValue(MDL_ANGLE_PSI, psi, id)
            MDout.setValue(MDL_SHIFT_X, x,id)
            MDout.setValue(MDL_SHIFT_Y, y,id)
            MDout.setValue(MDL_FLIP,flip,id)
            MDout.setValue(MDL_ENABLED,1,id)
            n+=1
        MDout.write(fnOut,MD_APPEND)
        
        # Actually create the differences
        img=Image()
        imgRef=Image()
        if hasCTF and volumeIsCTFCorrected:
            fnAcquisition=findAcquisitionInfo(fnOut)
            if not fnAcquisition:
                hasCTF=False
            else:
                mdAcquisition=MetaData(fnAcquisition)
                Ts=mdAcquisition.getValue(MDL_SAMPLINGRATE,mdAcquisition.firstObject())
    
        for i in MDout:
            img.readApplyGeo(MDout,i)
            imgRef.read(MDout.getValue(MDL_IMAGE_REF,i))
            if hasCTF and volumeIsCTFCorrected:
                fnCTF=MDout.getValue(MDL_CTF_MODEL,i)
                imgRef.applyCTF(fnCTF,Ts)
                img.convert2DataType(DT_DOUBLE)
            imgDiff=img-imgRef
            imgDiff.write(MDout.getValue(MDL_IMAGE1,i))


