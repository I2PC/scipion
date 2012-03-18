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
import emx_struct
"""
1) read star file in star object
2) parse star object an assign to data structure EMX/XMIPP
3) convert EMX to XMIPP (or vice versa)
4) parse XMIPP/EMX data structure to star object 
5) save star file
"""
from emx_struct import ParticleAlignmentEmx,\
                       ParticleAlignmentXmd,\
                       BlockNamesEMX,\
                       BlockNamesXMD,\
                       prefix_micrograph
from protlib_emx import EmxBase
import CifFile
import StarFile

##########################################################################
#   Class Related to Particle Picking Conversion
##########################################################################
class ParticleAlignmentConverter(EmxBase):    
    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
                       
    def readBlocksEMX(self):
        self.myStructEMX = ParticleAlignmentEmx()
        label = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListEMX.append(BlockNamesEMX(BlockName=blockName,
                                                    size=len(_auxList)))

#        for item in self.blockNameListEMX:
#            print item.BlockName, item.size
#        exit(1)

    def readBlocksXMD(self):
        self.myStructXMD = ParticlePickingStructXmd()
        label = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListXMD.append(BlockNamesXMD(BlockName=blockName,
                                                    size=len(_auxList)))
    def startObject2EMX(self):
        for item in self.blockNameListEMX:
            blockName = item.BlockName
            url = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _url=self.inMetadata[blockName][url]
            
            transformation_matrix_1_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[1][0]
            _transformation_matrix_1_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_1)

            transformation_matrix_1_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[2][0]
            _transformation_matrix_1_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_2)
            transformation_matrix_1_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[3][0]
            _transformation_matrix_1_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_3)
            transformation_matrix_offset_x = self.myStructEMX.prefix + self.myStructEMX._fields_[4][0]
            _transformation_matrix_offset_x=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_x)

            transformation_matrix_2_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[5][0]
            _transformation_matrix_2_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_1)
            transformation_matrix_2_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[6][0]
            _transformation_matrix_2_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_2)
            transformation_matrix_2_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[7][0]
            _transformation_matrix_2_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_3)
            transformation_matrix_offset_y = self.myStructEMX.prefix + self.myStructEMX._fields_[8][0]
            _transformation_matrix_offset_y=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_y)

            transformation_matrix_3_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[9][0]
            _transformation_matrix_3_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_1)
            transformation_matrix_3_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[10][0]
            _transformation_matrix_3_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_2)
            transformation_matrix_3_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[11][0]
            _transformation_matrix_3_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_3)
            transformation_matrix_offset_z = self.myStructEMX.prefix + self.myStructEMX._fields_[12][0]
            _transformation_matrix_offset_z=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_z)

            enable = self.myStructEMX.prefix + self.myStructEMX._fields_[13][0]
            _enable=self.inMetadata[blockName].GetLoopItem(enable)
            FOM = self.myStructEMX.prefix + self.myStructEMX._fields_[14][0]
            _FOM=self.inMetadata[blockName].GetLoopItem(FOM)


            for i in range(len(_url)):
                self.itemNameListEMX.append(ParticleAlignmentEmx
                                (
                                url=_url[i],
                                transformation_matrix_1_1=float(_transformation_matrix_1_1[i]),
                                transformation_matrix_1_2=float(_transformation_matrix_1_2[i]),
                                transformation_matrix_1_3=float(_transformation_matrix_1_3[i]),
                                transformation_matrix_offset_x=float(_transformation_matrix_offset_x[i]),
                                transformation_matrix_2_1=float(_transformation_matrix_2_1[i]),
                                transformation_matrix_2_2=float(_transformation_matrix_2_2[i]),
                                transformation_matrix_2_3=float(_transformation_matrix_2_3[i]),
                                transformation_matrix_offset_y=float(_transformation_matrix_offset_y[i]),
                                transformation_matrix_3_1=float(_transformation_matrix_3_1[i]),
                                transformation_matrix_3_2=float(_transformation_matrix_3_2[i]),
                                transformation_matrix_3_3=float(_transformation_matrix_3_3[i]),
                                transformation_matrix_offset_z=float(_transformation_matrix_offset_z[i]),
                                enable=int(_enable[i]),
                                FOM=float(_FOM[i])
                                )
                              )
#        for item in self.itemNameListEMX:
#            print item.url, \
#                  item.transformation_matrix_1_1,\
#                  item.transformation_matrix_offset_x,\
#                  item.FOM
#        exit(1)

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
            image=item.url,
            angleRot,
            angleTilt,
            anglePsi,
            shiftX,
            shiftY,
            shiftZ,
            flip,
            scale,
            enabled,
            fom

                                            Xcoor=int(round(item.coordinate_x)),
                                            Ycoor=int(round(item.coordinate_y)),
                                            #BlockName = "kk"+item.BlockName
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

#        for item in self.itemNameListXMD:
#             print item.Xcoor, item.Ycoor, item.BlockName
#        exit(1)
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
    