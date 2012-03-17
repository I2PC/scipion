#!/usr/bin/env xmipp_python
'''
/***************************************************************************
 * Authors:     Roberto Marabini Ruiz
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
'''
import CifFile
import StarFile
import sys
from transformations import *
import numpy
#from emx_struct import ParticlePickingStructEmx,\
#                       ParticlePickingStructXmd,\
#                       CtfMicrographStructEmx,\
#                       CtfMicrographStructXmd,\
#                       BlockNamesEMX,\
#                       BlockNamesXMD
                       

##########################################################################
#   General Class for Star data handling
##########################################################################

class EmxBase:
    
    blockNameListEMX=[]
    blockNameListXMD=[]
    itemNameListEMX=[]
    itemNameListXMD=[]
    """some constants"""
    xmippStartVersion = 'XMIPP_STAR_1'
    emxVersion        = 'EMX1.0'
    contactMail       = 'xmipp@cnb.csic.es'
        
    def checkVersion(self):
        """ Check first line for EMX or XMIPP magic word. Abort if
        neither of these two words are available"""
        import os
        if not os.path.exists(self.inputFileName):
            print >> sys.stderr, "File: ", self.inputFileName, "does not exists."
            exit(0)
        #set buffering size to 0 otherwise stdin is read in a buffer
        #and it is unavailable for the next open
        fin = open(self.inputFileName, "r",0)
        firstLine = fin.readline()
        #Is EMX?
        result = firstLine.find(self.emxVersion)
        if result != -1:
            fin.close()
            return (True)
        #is xmipp?
        result = firstLine.find(self.xmippStartVersion)
        if result != -1:
            fin.close()
            return (False) # isEMX
        #if not emx or xmipp abort
        print >> sys.stderr,  "Error: Metadata Files should contain the string ",\
               self.emxVersion, "or", self.xmippStartVersion, \
            " in the first line.  Exiting program.\n", \
            "First line: ", firstLine
        exit(1)
    
    def saveFileXMD(self):
        """Auxiliary funtion for saving metadata as xmipp file"""
        comment  =   "# XMIPP_STAR_1 *"
        comment += "\n##########################################################################"         
        comment +=  "\n#  Converted from " + self.emxVersion + " to " + self.xmippStartVersion
        comment +=  "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################" 
        outfile = open(self.outputFileName,"w")
        outfile.write(self.outMetadata.WriteOut(comment=comment,_email=self.contactMail))
        
    def saveFilePerBlockXMD(self):
        """Auxiliary funtion for saving metadata as xmipp file"""
        comment  =   "# XMIPP_STAR_1 *"
        comment += "\n##########################################################################"         
        comment +=  "\n#  Converted from " + self.emxVersion + " to " + self.xmippStartVersion
        counter = 1
        for blockName in self.outMetadata.keys():
            comment1  =  comment + "\n#  Inputfile: " + blockName
            comment1 += "@" + self.inputFileName
            comment1 += "\n##########################################################################" 
            filenameAux = self.outputFileName 
            if filenameAux != '/dev/stdout':
                filenameAux += "_" + str(counter).zfill(4)
            counter += 1
            outfile = open(filenameAux,"w")
            block=self.outMetadata[blockName]
            cf = CifFile.CifFile() 
            cf[blockName] = block 
            outfile.write(cf.WriteOut(comment=comment1,_email=self.contactMail))

    def saveFilePerLineXMD(self,attribute):
        """Auxiliary funtion for saving metadata as xmipp file"""
        comment  =   "# XMIPP_STAR_1 *"
        comment += "\n##########################################################################"         
        comment +=  "\n#  Converted from " + self.emxVersion + " to " + self.xmippStartVersion
        counter = 1
        for blockName in self.outMetadata.keys():
            block=self.outMetadata[blockName]
            lb = block.GetLoop(attribute)
            cf = CifFile.CifFile()
            for line in lb:
                comment1  =  comment + "\n#  Inputfile: " + blockName + "_"
                comment1 += attribute
                comment1 += "@" + self.inputFileName
                comment1 += "\n##########################################################################" 
                filenameAux = self.outputFileName 
                if filenameAux != '/dev/stdout':
                    filenameAux += "_" + str(counter).zfill(4)
                counter += 1
                outfile = open(filenameAux,"w")
                #convert line to list of list
                lineList=[]
                for i in line:
                   lineList.append([i])
                cf.clear()
                myblock = CifFile.CifBlock()
                cf[blockName] = myblock
                cf[blockName].AddCifItem((block.loopnames(),[lineList]))
                outfile.write(cf.WriteOut(comment=comment1,_email=self.contactMail))
        
    def saveFileEMX(self):
        """Auxiliary funtion for saving metadata as emx file"""
        outfile = open(self.outputFileName,"w")
        comment =  "#  Converted from " + self.xmippStartVersion + " to " + self.emxVersion
        comment += "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################"
        outfile.write(self.outMetadata.WriteOut(_add=True,comment=comment,_email=self.contactMail))

    def run(self):
        self.emx2xmipp = self.checkVersion()
        #do not change the order: first checkversion then ciffile
        #otherwise stdin will be lost
        #1 read star file in star object
        self.inMetadata     = CifFile.CifFile(self.inputFileName)
        self.outMetadata    = CifFile.CifFile()

        if self.emx2xmipp:
            #2 parse star object to EMX struct
            self.readBlocksEMX()
            #2bis parse star object to EMX struct
            self.startObject2EMX()
            #3 convert from EMX to XMIPP
            self.EMX2XMD()
            #convert xmipp to star object- this time save as xmipp star format 1
            self.XMD2startObject()
            #save file
            self.saveFileXMD()
        else:
            #2 parse star object to EMX struct
            self.readBlocksXMD()
            #2bis parse star object to EMX struct
            self.startObject2XMD()
            #3 convert from EMX to XMIPP
            self.XMD2EMX()
            #convert xmipp to star object- this time save as xmipp star format 1
            self.EMX2startObject()
            #save file
            self.saveFileEMX()



###########################################################################
##   Class Related to Alignment Conversion
###########################################################################
class ParticleAlignmentConverter(EmxBase):    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
        self.blockNameListEMX=[]
        self.blockNameListXMD=[]
        self.itemNameListEMX=[]
        self.itemNameListXMD=[]
    
    def run(self):
        self.emx2xmipp = self.checkVersion()
        #do not change the order: first checkversion then ciffile
        #otherwise stdin will be lost
        #1 read star file in star object
        self.inMetadata     = CifFile.CifFile(self.inputFileName)
        self.outMetadata    = CifFile.CifFile()

        if self.emx2xmipp:
            #2 parse star object to EMX struct
            self.readBlocksEMX()
            #2bis parse star object to EMX struct
            self.startObject2EMX()
#            #3 convert from EMX to XMIPP
#            self.EMX2XMD()
#            #convert xmipp to star object- this time save as xmipp star format 1
#            self.XMD2startObject()
#            #save file
#            self.saveFileXMD()
#        else:
            #2 parse star object to EMX struct
#            self.readBlocksXMD()
#            #2bis parse star object to EMX struct
#            self.startObject2XMD()
#            #3 convert from EMX to XMIPP
#            self.XMD2EMX()
#            #convert xmipp to star object- this time save as xmipp star format 1
#            self.EMX2startObject()
#            #save file
#            self.saveFileEMX()
                       
    def readBlocksEMX(self):
        self.myStructEMX = ParticlePickingStructEmx()
        label = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListEMX.append(BlockNamesEMX(BlockName=blockName,
                                                    size=len(_auxList)))

    def readBlocksXMD(self):
        self.myStructXMD = ParticlePickingStructXmd()
        label = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListXMD.append(BlockNamesXMD(BlockName=blockName,
                                                    size=len(_auxList)))
    def startObject2EMX(self):
        for blockName in self.inMetadata.keys():
           #for i in range(len(self.myStruct._fields_)-1):
            url = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _url=self.inMetadata[blockName].GetLoopItem(url)

            transformation_matrix_1_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_1_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_1)
            transformation_matrix_1_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_1_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_2)
            transformation_matrix_1_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_1_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_3)
            transformation_matrix_offset_x = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_offset_x=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_x)

            transformation_matrix_2_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_2_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_1)
            transformation_matrix_2_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_2_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_2)
            transformation_matrix_2_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_2_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_3)
            transformation_matrix_offset_y = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_offset_y=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_y)

            transformation_matrix_3_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_3_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_1)
            transformation_matrix_3_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_3_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_2)
            transformation_matrix_3_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_3_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_3)
            transformation_matrix_offset_z = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _transformation_matrix_offset_z=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_z)
#enable,
#FOM,
#BlockName

            url = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _url=self.inMetadata[blockName].GetLoopItem(url)

            url = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _url=self.inMetadata[blockName].GetLoopItem(url)

            url = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _url=self.inMetadata[blockName].GetLoopItem(url)

            for i in range(len(_auxCoordenateX)):
                self.itemNameListEMX.append(ParticlePickingStructEmx
                                           (
                                            coordinate_x=float(_auxCoordenateX[i]),
                                            coordinate_y=float(_auxCoordenateY[i]),
                                            BlockName = blockName
                                           )
                                        )
    def startObject2XMD(self):
        for blockName in self.inMetadata.keys():
           #for i in range(len(self.myStruct._fields_)-1):
            self.Xcoor = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
            _auxCoordenateX=self.inMetadata[blockName].GetLoopItem(self.Xcoor)
            self.Ycoor = self.myStructXMD.prefix + self.myStructXMD._fields_[1][0]
            _auxCoordenateY=self.inMetadata[blockName].GetLoopItem(self.Ycoor)
            
            for i in range(len(_auxCoordenateX)):
                self.itemNameListXMD.append(ParticlePickingStructXmd
                                           (
                                            Xcoor=int(_auxCoordenateX[i]),
                                            Ycoor=int(_auxCoordenateY[i]),
                                            BlockName = blockName
                                           )
                                        )
#        for item in self.itemNameListXMD:
#            print item.Xcoor, item.Ycoor, item.BlockName
    def EMX2XMD(self):
        
        for item in self.itemNameListEMX:
            self.itemNameListXMD.append(ParticlePickingStructXmd
                                           (
                                            Xcoor=int(round(item.coordinate_x)),
                                            Ycoor=int(round(item.coordinate_y)),
                                            #BlockName = "kk"+item.BlockName
                                            BlockName = item.BlockName
                                           )
                                          )
        #here you may cahnge the block name if needed
        for item in self.blockNameListEMX:
            self.blockNameListXMD.append(BlockNamesXMD(
#                                                       BlockName="kk"+item.BlockName,
                                                       BlockName=item.BlockName,
                                                       size=item.size)
                                         )

#        for item in self.itemNameListXMD:
#             print item.Xcoor, item.Ycoor, item.BlockName

    def XMD2EMX(self):
        
        for item in self.itemNameListXMD:
            self.itemNameListEMX.append(ParticlePickingStructEmx
                                           (
                                            coordinate_x=float(item.Xcoor),
                                            coordinate_y=float(item.Ycoor),
                                            #BlockName = "kk"+item.BlockName
                                            BlockName = item.BlockName
                                           )
                                          )
        #here you may change the block name if needed
        for item in self.blockNameListXMD:
            self.blockNameListEMX.append(BlockNamesEMX(
#                                                       BlockName="kk"+item.BlockName,
                                                       BlockName=item.BlockName,
                                                       size=item.size)
                                         )

#        for item in self.itemNameListXMD:
#             print item.Xcoor, item.Ycoor, item.BlockName

    def XMD2startObject(self):
        #create blocks in empty star file
        for block in self.blockNameListXMD:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = ParticlePickingStructXmd()
        _Xcorr_ = xmdStruct.prefix + xmdStruct._fields_[0][0]
        _Ycoor_ = xmdStruct.prefix + xmdStruct._fields_[1][0]
        #fill block
        for block in self.blockNameListXMD:
            XcoorList =[]
            YcoorList =[]
            for item in self.itemNameListXMD:
                if block.BlockName == item.BlockName:
                    XcoorList.append(str(item.Xcoor))
                    YcoorList.append(str(item.Ycoor))
            self.outMetadata[block.BlockName].AddCifItem(([[_Xcorr_,_Ycoor_]],[[XcoorList,YcoorList]]))

    def EMX2startObject(self):
        #create blocks in empty star file
        for block in self.blockNameListEMX:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = ParticlePickingStructEmx()
        _Xcorr_ = xmdStruct.prefix + xmdStruct._fields_[0][0]
        _Ycoor_ = xmdStruct.prefix + xmdStruct._fields_[1][0]
        #fill block
        for block in self.blockNameListEMX:
            XcoorList =[]
            YcoorList =[]
            for item in self.itemNameListEMX:
                if block.BlockName == item.BlockName:
                    XcoorList.append(str(item.coordinate_x))
                    YcoorList.append(str(item.coordinate_y))
            self.outMetadata[block.BlockName].AddCifItem(([[_Xcorr_,_Ycoor_]],[[XcoorList,YcoorList]]))
#

#    
#    def __init__(self, inputFn, outputFn):
#        self.inputFileName  = inputFn
#        self.outputFileName = outputFn       
#        #next two lines should go after checkVersion
#    
#    def run(self):
#        emx2xmipp = self.checkVersion()
#        self.inMetadata     = CifFile.CifFile(self.inputFileName)
#        self.outMetadata    = CifFile.CifFile()
#
#        if emx2xmipp:
#            self.convertAllBlocksEMX2XMIPP()
#            self.saveFileEMX2XMIPP()
#        else:
#            self.convertAllBlocksXMIPP2EMX()
#            self.saveFileXMIPP2EMX()          
#                       
#    def convertAllBlocksEMX2XMIPP(self):
#        """loop over the blocks and write them"""
#        for blockName in self.inMetadata.keys():
#            #create output block
#            myblock = CifFile.CifBlock()
#            self.outMetadata[blockName] = myblock
#            self.cbOut = self.outMetadata[blockName]
#            self.convertLoopEMX2XMIPP(blockName)
#            
#    def convertAllBlocksXMIPP2EMX(self):
#        """loop over the blocks and write them"""
#        for keyName in self.inMetadata.keys():
#            #create output block
#            myblock = CifFile.CifBlock()
#            self.outMetadata[keyName] = myblock
#            self.cbOut = self.outMetadata[keyName]
#            self.cbOut['_emx_micrograph.url'] = keyName
#            self.createDataHeaderXMIPP2EMX(keyName)
#            self.convertLoopXMIPP2EMX(keyName)
#            
#    def createDataHeaderXMIPP2EMX(self,micrographName):
#        """Data header is the xmipp data block name with the right label"""
#        self.cbOut['_emx_micrograph.url'] = micrographName
#
#    def convertLoopEMX2XMIPP(self,blockName):
#        #get item names and values
#        loopitems = self.inMetadata[blockName].GetLoop(self.needed_itemsEMX[0])  
#        #get block list
#        ListBlockKeys = loopitems.keys()
#        #get loop length
#        blockLenght = len (self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[1]))
#        #image is compulsory 
#        _imageUrl   = self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0])
#        #store matrices here
#        ListEulerMatrices = []
#        for loopitem in range(blockLenght):
#            ListEulerMatrices.append(numpy.identity(4))
#        #fill shifts
#        if (self.needed_itemsEMX[4]  in ListBlockKeys):
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[4])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][0][3]=aux[loopitem]
#        
#        if (self.needed_itemsEMX[8]  in ListBlockKeys):
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[8])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][1][3]=aux[loopitem]
#        
#        if (self.needed_itemsEMX[12]  in ListBlockKeys):
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[12])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][2][3]=aux[loopitem]
#        # fill 3D rotations
#        _2Dmatrix=True
#        if (self.needed_itemsEMX[3]  in ListBlockKeys):
#            _2Dmatrix=False
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[3])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][0][2]=aux[loopitem]
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[7])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][1][2]=aux[loopitem]
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[9])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][2][0]=aux[loopitem]
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[10])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][2][1]=aux[loopitem]
#            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[11])
#            for loopitem in range(blockLenght):
#                ListEulerMatrices[loopitem][2][2]=aux[loopitem]
#        #2D rotations
#        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[1])
#        for loopitem in range(blockLenght):
#            ListEulerMatrices[loopitem][0][0]=aux[loopitem]
#        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[2])
#        for loopitem in range(blockLenght):
#            ListEulerMatrices[loopitem][0][1]=aux[loopitem]
#        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[5])
#        for loopitem in range(blockLenght):
#            ListEulerMatrices[loopitem][1][0]=aux[loopitem]
#        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[6])
#        for loopitem in range(blockLenght):
#            ListEulerMatrices[loopitem][1][1]=aux[loopitem]
#
#        #handle enable
#        _enable = [1]*blockLenght
#        if (self.needed_itemsEMX[13]  in ListBlockKeys):
#            _enable= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[13])#enable
#        #handle FOM
#        _fom = [1.]*blockLenght
#        if (self.needed_itemsEMX[14]  in ListBlockKeys):
#            _fom= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[14])#enable
#        
#        #availabel in _imageUrl
#        #_image     = []
#        _angleRot  = []
#        _angleTilt = []
#        _anglePsi  = []
#        _shiftX    = []
#        _shiftY    = []
#        _shiftZ    = []
#        _flip      = [0]*blockLenght
#        _scale     = []
#
#        for loopitem in range(blockLenght):
#            #create 4x4 matrix
#            scale, shear, angles, trans, persp = decompose_matrix(ListEulerMatrices[loopitem])
#            #images should not present shear or persp
#            if (math.fabs(numpy.linalg.norm(shear))>smallNumber):
#                print >> sys.stderr,  "The input matrix ", ListEulerMatrices[loopitem],\
#                "presents shear. This is not supported"
#                exit(1)
#            if (math.fabs(numpy.linalg.norm(persp)-1.)>smallNumber):
#                print >> sys.stderr,  "The input matrix ", ListEulerMatrices[loopitem],\
#                "presents perpective. This is not supported"
#                exit(1)
#
#            _angleRot.append ("%0.4f" % (angles[0]* 180./math.pi))
#            _angleTilt.append("%0.4f" % (angles[1]* 180./math.pi))
#            _anglePsi.append ("%0.4f" % (angles[2]* 180./math.pi))
#            _shiftX.append ("%0.4f" % (trans[0]))
#            _shiftY.append ("%0.4f" % (trans[1]))
#            _shiftZ.append ("%0.4f" % (trans[2]))
#            #flip only in 2D
#            if(_2Dmatrix):
#                if (numpy.linalg.det()< 0):
#                    _flip[loopitem] = 1
#            #check if different in each direction
#            scaleAverage = (scale[0]+scale[1]+scale[2])/3.
#            if(math.fabs(scaleAverage-scale[0])>smallNumber):
#                print >> sys.stderr,  "Reading Iamge:", _imageUrl[loopitem], 
#                print >> sys.stderr,  "Different scale factor in each axis.", scale, "This is not supported", 
#                print >> sys.stderr,  "scale along x axis will be assigned"
#            _scale.append ("%0.4f" % (scaleAverage))
#
##        self.cbOut.AddCifItem(([self.needed_itemsXMIPP],[[_XList,_YList]])) 
#            
#        self.cbOut.AddCifItem(( [self.needed_itemsXMIPP], 
#                                [[
#                                _imageUrl,
#                                _angleRot,
#                                _angleTilt,
#                                _anglePsi,
#                                _shiftX,
#                                _shiftY,
#                                _shiftZ,
#                                _flip,
#                                _scale,
#                                _enable,
#                                _fom
#                                ]]
#                              )) 
#       
#    def convertLoopXMIPP2EMX(self,blockName):
#        #get item names and values
#        loopitems = self.inMetadata[blockName].GetLoop(self.needed_itemsXMIPP[0])  
#        #get block list
#        ListBlockKeys = loopitems.keys()
#        #get loop length
#        blockLenght = len (self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[0]))        
#        #image is compulsory 
#        _imageUrl   = self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[0])
#        #store intermediate results here
#        shift=[]
#        angle=[]
#        scale=[]
#
#        _x=[]
#        _y=[]
#        _z=[]
#        for loopitem in range(blockLenght):
#            _x.append(0)
#            _y.append(0)
#            _z.append(0)
#        
#        #fill shifts
#        if (self.needed_itemsXMIPP[4]  in ListBlockKeys):
#            _x= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[4])
#        if (self.needed_itemsXMIPP[5]  in ListBlockKeys):
#            _y= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[5])
#        if (self.needed_itemsXMIPP[6]  in ListBlockKeys):
#            _z= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[6])
#        for loopitem in range(blockLenght):
#            shift.append((_x[loopitem],_y[loopitem],_z[loopitem]))
#
#        _rot=[]
#        _tilt=[]
#        _psi=[]
#        
#        for loopitem in range(blockLenght):
#            _rot.append(0)
#            _tilt.append(0)
#            _psi.append(0)
#        #fill rotations
#        if (self.needed_itemsXMIPP[1]  in ListBlockKeys):
#            _rot= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[1])
#        if (self.needed_itemsXMIPP[2]  in ListBlockKeys):
#            _tilt= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[2])
#        if (self.needed_itemsXMIPP[3]  in ListBlockKeys):
#            _psi= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[3])
#        for loopitem in range(blockLenght):
#            angle.append( 
#                          (
#                           (float(_rot[loopitem] )) *math.pi/180.,\
#                           (float(_tilt[loopitem])) *math.pi/180.,\
#                           (float(_psi[loopitem] )) *math.pi/180.
#                          )  
#                        )
#                        
#                        
#
#        _scale=[]
#        
#        for loopitem in range(blockLenght):
#            _scale.append(1.)
#        #fill scale
#        if (self.needed_itemsXMIPP[8]  in ListBlockKeys):
#            _scale= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[8])
#        for loopitem in range(blockLenght):
#            scale.append((_scale[loopitem],_scale[loopitem],_scale[loopitem]))
#        #enable
#        _enable=[]
#        for loopitem in range(blockLenght):
#            _enable.append(1)
#        if (self.needed_itemsXMIPP[9]  in ListBlockKeys):
#            _enable= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[9])
#        #fom
#        _fom=[]
#        for loopitem in range(blockLenght):
#            _fom.append(1.)
#        if (self.needed_itemsXMIPP[10]  in ListBlockKeys):
#            _fom= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[10])
#        
#        #convert this in matrices
#        #store matrices here
#        ListEulerMatrices = []
#
#        for loopitem in range(blockLenght):
#            ListEulerMatrices.append(compose_matrix(scale=scale[loopitem],
#                    angles=angle[loopitem],
#                    translate=shift[loopitem]))
#
#        #Convert matrices to emx
#        _emx_particle_url = _imageUrl
#        _emx_particle_transformation_matrix_offset_x = _x
#        _emx_particle_transformation_matrix_offset_y = _y
#        _emx_particle_transformation_matrix_offset_z = _z
#        _emx_particle_transformation_matrix_1_1=[]
#        _emx_particle_transformation_matrix_1_2=[]
#        _emx_particle_transformation_matrix_1_3=[]
#        _emx_particle_transformation_matrix_2_1=[]
#        _emx_particle_transformation_matrix_2_2=[]
#        _emx_particle_transformation_matrix_2_3=[]
#        _emx_particle_transformation_matrix_3_1=[]
#        _emx_particle_transformation_matrix_3_2=[]
#        _emx_particle_transformation_matrix_3_3=[]
#        
#        for loopitem in range(blockLenght):
#            _emx_particle_transformation_matrix_1_1.append(ListEulerMatrices[loopitem][0][0])
#            _emx_particle_transformation_matrix_1_2.append(ListEulerMatrices[loopitem][0][1])
#            _emx_particle_transformation_matrix_1_3.append(ListEulerMatrices[loopitem][0][2])
#            _emx_particle_transformation_matrix_2_1.append(ListEulerMatrices[loopitem][1][0])
#            _emx_particle_transformation_matrix_2_2.append(ListEulerMatrices[loopitem][1][1])
#            _emx_particle_transformation_matrix_2_3.append(ListEulerMatrices[loopitem][1][2])
#            _emx_particle_transformation_matrix_3_1.append(ListEulerMatrices[loopitem][2][0])
#            _emx_particle_transformation_matrix_3_2.append(ListEulerMatrices[loopitem][2][1])
#            _emx_particle_transformation_matrix_3_3.append(ListEulerMatrices[loopitem][2][2])
#        
#        _emx_particle_enable  = _enable
#        _emx_particle_FOM     = _fom
#        self.cbOut.AddCifItem(( [self.needed_itemsEMX], [
#                               [
#                                _emx_particle_url,
#                                _emx_particle_transformation_matrix_1_1,
#                                _emx_particle_transformation_matrix_1_2,
#                                _emx_particle_transformation_matrix_1_3,
#                                _emx_particle_transformation_matrix_offset_x,
#                                _emx_particle_transformation_matrix_2_1,
#                                _emx_particle_transformation_matrix_2_2,
#                                _emx_particle_transformation_matrix_2_3,
#                                _emx_particle_transformation_matrix_offset_y,
#                                _emx_particle_transformation_matrix_3_1,
#                                _emx_particle_transformation_matrix_3_2,
#                                _emx_particle_transformation_matrix_3_3,
#                                _emx_particle_transformation_matrix_offset_z,
#                                _emx_particle_enable,
#                                _emx_particle_FOM
#                                ]
#                              ])) 
#
###########################################################################
##   Class Related to Particle Picking Conversion
###########################################################################
#class ParticleClassConverter(EmxBase):    
#    needed_itemsXMIPP =  (
#          "_image",
#          "_ref"#int
#          )
#    needed_itemsEMX = (
#          "_emx_particle.url"
#          )
#
#    
#    def __init__(self, inputFn, outputFn):
#        self.inputFileName  = inputFn
#        self.outputFileName = outputFn
#
#        self.imageList = []
#        self.refList   = []
#        #next two lines should go after checkVersion
#    
#    def run(self):
#        emx2xmipp = self.checkVersion()
#        self.inMetadata     = CifFile.CifFile(self.inputFileName)
#        self.outMetadata    = CifFile.CifFile()
#
#        if emx2xmipp:
#            self.convertAllBlocksEMX2XMIPP()
#            self.saveFileEMX2XMIPP()
#        else:
#            self.convertAllBlocksXMIPP2EMX()
#            self.saveFileXMIPP2EMX()          
#                       
#    def convertAllBlocksEMX2XMIPP(self):
#        """loop over the blocks and write them"""
#        ref=1
#        myblock = CifFile.CifBlock()
#        self.outMetadata['class'] = myblock
#        self.cbOut = self.outMetadata['class']#alias
#        #self.cbOut.AddCifItem(self.needed_itemsEMX)
#
#        for blockName in self.inMetadata.keys():
#            self.convertLoopEMX2XMIPP(blockName,ref)
#            ref += 1
#        self.cbOut.AddCifItem(( [self.needed_itemsXMIPP],[[self.imageList,self.refList]] )) 
#    
#    def convertAllBlocksXMIPP2EMX(self):
#        """loop over the blocks and write them"""
#        #get all different references from xmipp
#        #just one block
#        _referenceList = self.inMetadata.first_block()
#        blockNameXmipp = self.inMetadata.keys()[0]
#        _referenceList = _referenceList.GetLoopItem(self.needed_itemsXMIPP[1])
#        #get a list with unique classes
#        #_referenceList = self.inMetadata[blockNameXmipp].GetLoopItem(self.needed_itemsXMIPP[0])
#        _referenceList = list(set(_referenceList))
#        for blockNameEmx in _referenceList:
#            #create output block
#            myblock = CifFile.CifBlock()
#            self.outMetadata[blockNameEmx] = myblock
#            self.cbOut = self.outMetadata[blockNameEmx]#for emx
#            self.createDataHeaderXMIPP2EMX(blockNameEmx)
#            self.convertLoopXMIPP2EMX(blockNameXmipp,blockNameEmx)#for xmipp
#
#    def createDataHeaderXMIPP2EMX(self,class_id):
#        """Data header is the xmipp data block name with the right label"""
#        self.cbOut['_emx_class.id'] = class_id
#            
#    def convertLoopEMX2XMIPP(self,blockName,ref):
#         _imageList = self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX)
#         _ref=[ref]*len(_imageList)
#         self.imageList += _imageList
#         self.refList   += _ref
#         #add to cif should be done at the end
#         
#    def convertLoopXMIPP2EMX(self,blockNameXmipp,blockNameEmx):
#        #get item names and values
#        loopimages  = self.inMetadata[blockNameXmipp].GetLoopItem(self.needed_itemsXMIPP[0])  
#        loopclasses = self.inMetadata[blockNameXmipp].GetLoopItem(self.needed_itemsXMIPP[1])  
#        #get block list
#        _imageList = []
#        for imageItem in range (len(loopclasses)):
#            if(loopclasses[imageItem]==blockNameEmx):
#                _imageList.append(loopimages[imageItem])
#        
##        print [(self.needed_itemsEMX)], [
##                               [_imageList
##                                ]
##                              ]
#
#        self.cbOut.AddCifItem((self.needed_itemsEMX, _imageList )) 
#
###########################################################################
##   Class Related to CTF conversion: defocus per micrograph
###########################################################################
class CtfMicrographConverter(EmxBase):    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
        self.blockNameListEMX=[]
        self.blockNameListXMD=[]
        self.itemNameListEMX=[]
        self.itemNameListXMD=[]
    
    def run(self):
        self.emx2xmipp = self.checkVersion()
        #do not change the order: first checkversion then ciffile
        #otherwise stdin will be lost
        #1 read star file in star object
        self.inMetadata     = CifFile.CifFile(self.inputFileName)
        self.outMetadata    = CifFile.CifFile()

        if self.emx2xmipp:
            #2 parse star object to EMX struct
            self.readBlocksEMX()
            #2bis parse star object to EMX struct
            self.startObject2EMX()
            #3 convert from EMX to XMIPP
            self.EMX2XMD()
            #convert xmipp to star object- this time save as xmipp star format 1
            self.XMD2startObject()
            #save file
            self.saveFileXMD()
            #each ctf in a file
            myStructXmd = CtfMicrographStructXmd()
            label = myStructXmd.prefix + myStructXmd._fields_[0][0]
            self.saveFilePerLineXMD(label)

        else:
            #2 parse star object to EMX struct
            self.readBlocksXMD()
            #2bis parse star object to EMX struct
            self.startObject2XMD()
            #3 convert from EMX to XMIPP
            self.XMD2EMX()
            #convert xmipp to star object- this time save as xmipp star format 1
            self.EMX2startObject()
            #save file
            self.saveFileEMX()
                       
    def readBlocksEMX(self):
        self.myStructEMX = CtfMicrographStructEmx()
        label = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListEMX.append(BlockNamesEMX(BlockName=blockName,
                                                    size=len(_auxList)))

    def readBlocksXMD(self):
        """do not think there is more than one block"""
        self.myStructXMD = CtfMicrographStructXmd()
        label = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListXMD.append(BlockNamesXMD(BlockName=blockName,
                                                    size=len(_auxList)))
            
    def startObject2EMX(self):
        for blockName in self.inMetadata.keys():
            self.url = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _url=self.inMetadata[blockName].GetLoopItem(self.url)
            self.magnification = self.myStructEMX.prefix + self.myStructEMX._fields_[1][0]
            _magnification=self.inMetadata[blockName].GetLoopItem(self.magnification)
            self.scanner_pixel_size = self.myStructEMX.prefix + self.myStructEMX._fields_[2][0]
            _scanner_pixel_size=self.inMetadata[blockName].GetLoopItem(self.scanner_pixel_size)
            self.defocusU = self.myStructEMX.prefix + self.myStructEMX._fields_[3][0]
            _defocusU=self.inMetadata[blockName].GetLoopItem(self.defocusU)
            self.defocusV = self.myStructEMX.prefix + self.myStructEMX._fields_[4][0]
            _defocusV=self.inMetadata[blockName].GetLoopItem(self.defocusV)
            self.astigmatism_angle = self.myStructEMX.prefix + self.myStructEMX._fields_[5][0]
            _astigmatism_angle=self.inMetadata[blockName].GetLoopItem(self.astigmatism_angle)
            self.voltage = self.myStructEMX.prefix + self.myStructEMX._fields_[6][0]
            _voltage=self.inMetadata[blockName].GetLoopItem(self.voltage)
            self.Cs = self.myStructEMX.prefix + self.myStructEMX._fields_[7][0]
            _Cs=self.inMetadata[blockName].GetLoopItem(self.Cs)
            self.amplitude_contrast = self.myStructEMX.prefix + self.myStructEMX._fields_[8][0]
            _amplitude_contrast=self.inMetadata[blockName].GetLoopItem(self.amplitude_contrast)
            for i in range(len(_Cs)):
                self.itemNameListEMX.append(CtfMicrographStructEmx
                                           (
                                            url=_url[i],
                                            magnification=float(_magnification[i]),
                                            scanner_pixel_size=float(_scanner_pixel_size[i]),
                                            defocusU=float(_defocusU[i]),
                                            defocusV=float(_defocusV[i]),
                                            astigmatism_angle=float(_astigmatism_angle[i]),
                                            voltage=float(_voltage[i]),
                                            Cs=float(_Cs[i]),
                                            amplitude_contrast=float(_amplitude_contrast[i]),
                                            BlockName = blockName
                                           )
                                        )

    def startObject2XMD(self):
        for blockName in self.inMetadata.keys():

            image= self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
            _image=self.inMetadata[blockName].GetLoopItem(image)

            CTF_Sampling_rate= self.myStructXMD.prefix + self.myStructXMD._fields_[1][0]
            _CTF_Sampling_rate=self.inMetadata[blockName].GetLoopItem(CTF_Sampling_rate)

            CTF_Defocus_U= self.myStructXMD.prefix + self.myStructXMD._fields_[2][0]
            _CTF_Defocus_U=self.inMetadata[blockName].GetLoopItem(CTF_Defocus_U)

            CTF_Defocus_V= self.myStructXMD.prefix + self.myStructXMD._fields_[3][0]
            _CTF_Defocus_V=self.inMetadata[blockName].GetLoopItem(CTF_Defocus_V)

            CTF_Defocus_angle= self.myStructXMD.prefix + self.myStructXMD._fields_[4][0]
            _CTF_Defocus_angle=self.inMetadata[blockName].GetLoopItem(CTF_Defocus_angle)

            CTF_Voltage= self.myStructXMD.prefix + self.myStructXMD._fields_[5][0]
            _CTF_Voltage=self.inMetadata[blockName].GetLoopItem(CTF_Voltage)

            CTF_Spherical_aberration= self.myStructXMD.prefix + self.myStructXMD._fields_[6][0]
            _CTF_Spherical_aberration=self.inMetadata[blockName].GetLoopItem(CTF_Spherical_aberration)

            CTF_Q0= self.myStructXMD.prefix + self.myStructXMD._fields_[7][0]
            _CTF_Q0=self.inMetadata[blockName].GetLoopItem(CTF_Q0)

            
            for i in range(len(_CTF_Q0)):
                self.itemNameListXMD.append(CtfMicrographStructXmd
                                           (
                                            image=_image[i],
                                            CTF_Sampling_rate=float(_CTF_Sampling_rate[i]),
                                            CTF_Defocus_U=float(_CTF_Defocus_U[i]),
                                            CTF_Defocus_V=float(_CTF_Defocus_V[i]),
                                            CTF_Defocus_angle=float(_CTF_Defocus_angle[i]),
                                            CTF_Voltage=float(_CTF_Voltage[i]),
                                            CTF_Spherical_aberration=float(_CTF_Spherical_aberration[i]),
                                            CTF_Q0=float(_CTF_Q0[i]),
                                            BlockName = blockName
                                           )
                                        )
#        for item in self.itemNameListXMD:
#            print item.image, item.CTF_Spherical_aberration, item.BlockName
    def EMX2XMD(self):
        
        for item in self.itemNameListEMX:
            self.itemNameListXMD.append(CtfMicrographStructXmd
                                           (
                                            image=item.url,
                                            CTF_Sampling_rate=10000*item.scanner_pixel_size/item.magnification,
                                            CTF_Defocus_U=item.defocusU,
                                            CTF_Defocus_V=item.defocusV,
                                            CTF_Defocus_angle=item.astigmatism_angle,
                                            CTF_Voltage=item.voltage,
                                            CTF_Spherical_aberration=item.Cs,
                                            CTF_Q0=item.amplitude_contrast,
                                            BlockName = item.BlockName
                                           )
                                        )
        #here you may change the block name if needed
        for item in self.blockNameListEMX:
            self.blockNameListXMD.append(BlockNamesXMD(
#                                                       BlockName="kk"+item.BlockName,
                                                       BlockName=item.BlockName,
                                                       size=item.size)
                                         )

        

    def XMD2EMX(self):
        
        for item in self.itemNameListXMD:

            self.itemNameListEMX.append(CtfMicrographStructEmx
                                           (
                                            url = item.image,
                                            magnification=10000,#microns
                                            scanner_pixel_size=item.CTF_Sampling_rate,
                                            defocusU = item.CTF_Defocus_U,
                                            defocusV = item.CTF_Defocus_V,
                                            astigmatism_angle = item.CTF_Defocus_angle,
                                            voltage = item.CTF_Voltage,
                                            Cs = item.CTF_Spherical_aberration,
                                            amplitude_contrast= item.CTF_Q0,
                                            BlockName = item.BlockName
                                           )
                                          )
        #here you may change the block name if needed
        for item in self.blockNameListXMD:
            self.blockNameListEMX.append(BlockNamesEMX(
#                                                       BlockName="kk"+item.BlockName,
                                                       BlockName=item.BlockName,
                                                       size=item.size)
                                         )

#        for item in self.itemNameListEMX:
#             print item.magnification, item.defocusU, item.BlockName

    def XMD2startObject(self):
        #create blocks in empty star file
        for block in self.blockNameListXMD:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = CtfMicrographStructXmd()
        _image=xmdStruct.prefix + xmdStruct._fields_[0][0]
        _CTF_Sampling_rate=xmdStruct.prefix + xmdStruct._fields_[1][0]
        _CTF_Defocus_U=xmdStruct.prefix + xmdStruct._fields_[2][0]
        _CTF_Defocus_V=xmdStruct.prefix + xmdStruct._fields_[3][0]
        _CTF_Defocus_angle=xmdStruct.prefix + xmdStruct._fields_[4][0]
        _CTF_Voltage=xmdStruct.prefix + xmdStruct._fields_[5][0]
        _CTF_Spherical_aberration=xmdStruct.prefix + xmdStruct._fields_[6][0]
        _CTF_Q0=xmdStruct.prefix + xmdStruct._fields_[7][0]        #fill block
        
        for block in self.blockNameListXMD:
            imageList=[]
            CTF_Sampling_rateList=[]
            CTF_Defocus_UList=[]
            CTF_Defocus_VList=[]
            CTF_Defocus_angleList=[]
            CTF_VoltageList=[]
            CTF_Spherical_aberrationList=[]
            CTF_Q0List=[]
            for item in self.itemNameListXMD:
                if block.BlockName == item.BlockName:
                    imageList.append(str(item.image))
                    CTF_Sampling_rateList.append(str(item.CTF_Sampling_rate))
                    CTF_Defocus_UList.append(str(item.CTF_Defocus_U))
                    CTF_Defocus_VList.append(str(item.CTF_Defocus_V))
                    CTF_Defocus_angleList.append(str(item.CTF_Defocus_angle))
                    CTF_VoltageList.append(str(item.CTF_Voltage))
                    CTF_Spherical_aberrationList.append(str(item.CTF_Spherical_aberration))
                    CTF_Q0List.append(str(item.CTF_Q0))
            
            self.outMetadata[block.BlockName].AddCifItem(([[_image,
                                                            _CTF_Sampling_rate,
                                                            _CTF_Defocus_U,
                                                            _CTF_Defocus_V,
                                                            _CTF_Defocus_angle,
                                                            _CTF_Voltage,
                                                            _CTF_Spherical_aberration,
                                                            _CTF_Q0]],
                                                          [[imageList,
                                                            CTF_Sampling_rateList,
                                                            CTF_Defocus_UList,
                                                            CTF_Defocus_VList,
                                                            CTF_Defocus_angleList,
                                                            CTF_VoltageList,
                                                            CTF_Spherical_aberrationList,
                                                            CTF_Q0List
                                                            ]]
                                                          ))

    def EMX2startObject(self):
        #create blocks in empty star file
        for block in self.blockNameListEMX:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = CtfMicrographStructEmx()
        _url_ = xmdStruct.prefix + xmdStruct._fields_[0][0]
        _magnification_= xmdStruct.prefix + xmdStruct._fields_[1][0]
        _scanner_pixel_size_= xmdStruct.prefix + xmdStruct._fields_[2][0]
        _defocusU_ = xmdStruct.prefix + xmdStruct._fields_[3][0]
        _defocusV_ = xmdStruct.prefix + xmdStruct._fields_[4][0]
        _astigmatism__angle = xmdStruct.prefix + xmdStruct._fields_[5][0]
        _voltage_ = xmdStruct.prefix + xmdStruct._fields_[5][0]
        _Cs_ = xmdStruct.prefix + xmdStruct._fields_[6][0]
        _amplitude_contrast_= xmdStruct.prefix + xmdStruct._fields_[7][0]

        #fill block
        for block in self.blockNameListEMX:
             urlList = []
             magnificationList= []
             scanner_pixel_sizeList= []
             defocusUList = []
             defocusVList = []
             astigmatism_angleList = []
             voltageList = []
             CsList = []
             amplitude_contrastList= []

             for item in self.itemNameListEMX:
                if block.BlockName == item.BlockName:
                     urlList.append(str(item.url))
                     magnificationList.append(str(item.magnification))
                     scanner_pixel_sizeList.append(str(item.scanner_pixel_size))
                     defocusUList.append(str(item.defocusU))
                     defocusVList.append(str(item.defocusV))
                     astigmatism_angleList.append(str(item.astigmatism_angle))
                     voltageList.append(str(item.voltage))
                     CsList.append(str(item.Cs))
                     amplitude_contrastList.append(str(item.amplitude_contrast))

             self.outMetadata[block.BlockName].AddCifItem(([[
                                                    _url_ ,
                                                    _magnification_,
                                                    _scanner_pixel_size_,
                                                    _defocusU_ ,
                                                    _defocusV_ ,
                                                    _astigmatism__angle ,
                                                    _voltage_ ,
                                                    _Cs_ ,
                                                    _amplitude_contrast_
                                                    ]],[[
                                                     urlList ,
                                                     magnificationList,
                                                     scanner_pixel_sizeList,
                                                     defocusUList ,
                                                     defocusVList ,
                                                     astigmatism_angleList ,
                                                     voltageList ,
                                                     CsList ,
                                                     amplitude_contrastList
                                                         ]]))


###########################################################################
##   Class Related to CTF: defocus per image
###########################################################################
class CtfParticleConverter(EmxBase):    
    def __init__(self, inputFn, outputFn):
        print """Not implemented yet, EMX approach does not follows XMIPP philosophy
                This will not be implemented until XMIPP 3.0 has been released
              """
    def run(self):
        print """Not implemented yet, EMX approach does not follows XMIPP philosophy
                This will not be implemented until XMIPP 3.0 has been released
              """
        

#Example usage
#cat Test/particlePicking.xmd | ./batch_emx_coordinates.py
#cat Test/particlePicking.emx | ./batch_emx_coordinates.py        
if __name__ == '__main__':
    
    def command_line_options():
        """ add command line options here"""
        import optparse
        _usage = "usage: %prog [options] Example:   %prog -i input.xmd -o out.emx"
        parser = optparse.OptionParser(_usage)        
        parser.add_option("-i", "--input_filename", dest="inputFn",
                          default="/dev/stdin", type="string",
                          help="EMS or XMIPP Metadata file")                           
        parser.add_option("-o", "--output_filename", dest="outputFn",
                          default="/dev/stdout", type="string",
                          help="XMIPP or EMX Metadata file")   
        parser.add_option("-t", "--conversion_type", dest="type",
                          default="coordinates", type="string",
                          help="Possible types of conversion: coordinates, alignment, class, ctfMicrograph, ctfParticle")                          
       
        (options, args) = parser.parse_args()
        return(options.inputFn,options.outputFn,options.type)
        
    inputFn,outputFn,convType=command_line_options()
    if convType == 'coordinates':
        ParticlePickingConverter(inputFn, outputFn).run()  
    elif convType == 'alignment':
        ParticleAlignmentConverter(inputFn, outputFn).run()  
    elif convType == 'class':
        ParticleClassConverter(inputFn, outputFn).run()  
    elif convType == 'ctfMicrograph':
        CtfMicrographConverter(inputFn, outputFn).run()  
    elif convType == 'ctfParticle':
        CtfParticleConverter(inputFn, outputFn).run()  
    else:
        print >> sys.stderr,  "ERROR: Wrong mode: ", convType
        exit(0)
    
    
    
#        
