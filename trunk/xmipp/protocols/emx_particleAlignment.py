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
import sys
import numpy
import math
from transformations import decompose_matrix
from transformations import compose_matrix
from numpy.matrixlib.defmatrix import matrix
_pi=math.pi

"""
1) read star file in star object
2) parse star object an assign to data structure EMX/XMIPP
3) convert EMX to XMIPP (or vice versa)
4) parse XMIPP/EMX data structure to star object 
5) save star file
"""
from emx_struct import ParticleAlignmentStructEmx,\
                       ParticleAlignmentStructXmd,\
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
        self.myStructEMX = ParticleAlignmentStructEmx()
        label = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListEMX.append(BlockNamesEMX(BlockName=blockName,
                                                    size=len(_auxList)))

#        for item in self.blockNameListEMX:
#            print item.BlockName, item.size
#        exit(1)

    def readBlocksXMD(self):
        self.myStructXMD = ParticleAlignmentStructXmd()
        label = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListXMD.append(BlockNamesXMD(BlockName=blockName,
                                                    size=len(_auxList)))
#        for item in self.blockNameListXMD:
#            print item.BlockName, item.size
#        exit(1)
    def startObject2EMX(self):
        #get block, all should have the same labels
        labelNames = self.inMetadata.first_block().keys() 
#        print labelNames
#        exit(1)

        for item in self.blockNameListEMX:
            blockName = item.BlockName
            url = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _url=self.inMetadata[blockName][url]
            length = len(_url)
            transformation_matrix_1_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[1][0]
            if (transformation_matrix_1_1 in labelNames):
                _transformation_matrix_1_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_1)
            else:
                _transformation_matrix_1_1=[1.]*length
#            print _url,_transformation_matrix_1_1,length
#            print labelNames, transformation_matrix_1_1
#            exit(1)

            transformation_matrix_1_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[2][0]
            if (transformation_matrix_1_2 in labelNames):
                _transformation_matrix_1_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_2)
            else:
                _transformation_matrix_1_2=[0.]*length

            transformation_matrix_1_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[3][0]
            if (transformation_matrix_1_3 in labelNames):
                _transformation_matrix_1_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_1_3)
            else:
                _transformation_matrix_1_3=[0.]*length
                
            transformation_matrix_offset_x = self.myStructEMX.prefix + self.myStructEMX._fields_[4][0]
            if (transformation_matrix_offset_x in labelNames):
                _transformation_matrix_offset_x=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_x)
            else:
                _transformation_matrix_offset_x=[0.]*length

            transformation_matrix_2_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[5][0]
            if (transformation_matrix_2_1 in labelNames):
                _transformation_matrix_2_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_1)
            else:
                _transformation_matrix_2_1=[0.]*length

            transformation_matrix_2_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[6][0]
            if (transformation_matrix_2_2 in labelNames):
                _transformation_matrix_2_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_2)
            else:
                _transformation_matrix_2_2=[1.]*length

            transformation_matrix_2_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[7][0]
            if (transformation_matrix_2_3 in labelNames):
                _transformation_matrix_2_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_2_3)
            else:
                _transformation_matrix_2_3=[0.]*length

            transformation_matrix_offset_y = self.myStructEMX.prefix + self.myStructEMX._fields_[8][0]
            if (transformation_matrix_offset_y in labelNames):
                _transformation_matrix_offset_y=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_y)
            else:
                _transformation_matrix_offset_y=[0.]*length

            transformation_matrix_3_1 = self.myStructEMX.prefix + self.myStructEMX._fields_[9][0]
            if (transformation_matrix_3_1 in labelNames):
                _transformation_matrix_3_1=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_1)
            else:
                _transformation_matrix_3_1=[0.]*length
                
            transformation_matrix_3_2 = self.myStructEMX.prefix + self.myStructEMX._fields_[10][0]
            if (transformation_matrix_3_2 in labelNames):
                _transformation_matrix_3_2=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_2)
            else:
                _transformation_matrix_3_2=[0.]*length
                
            transformation_matrix_3_3 = self.myStructEMX.prefix + self.myStructEMX._fields_[11][0]
            if (transformation_matrix_3_3 in labelNames):
                _transformation_matrix_3_3=self.inMetadata[blockName].GetLoopItem(transformation_matrix_3_3)
            else:
                _transformation_matrix_3_3=[1.]*length
                
            transformation_matrix_offset_z = self.myStructEMX.prefix + self.myStructEMX._fields_[12][0]
            if (transformation_matrix_offset_z in labelNames):
                _transformation_matrix_offset_z=self.inMetadata[blockName].GetLoopItem(transformation_matrix_offset_z)
            else:
                _transformation_matrix_offset_z=[0.]*length

            enable = self.myStructEMX.prefix + self.myStructEMX._fields_[13][0]
            if (enable in labelNames):
                _enable=self.inMetadata[blockName].GetLoopItem(enable)
            else:
                _enable=[1]*length
                
            FOM = self.myStructEMX.prefix + self.myStructEMX._fields_[14][0]
            if (FOM in labelNames):
                _FOM=self.inMetadata[blockName].GetLoopItem(FOM)
            else:
                _FOM=[1.]*length


            for i in range(len(_url)):
                self.itemNameListEMX.append(ParticleAlignmentStructEmx
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
#                  item.transformation_matrix_1_2,\
#                  item.transformation_matrix_offset_x,\
#                  item.FOM
#        exit(1)

    def startObject2XMD(self):
        labelNames = self.inMetadata.first_block().keys() 
        for blockName in self.inMetadata.keys():
            image= self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
            imageList=self.inMetadata[blockName].GetLoopItem(image)
            
            length = len(imageList)

            angleRot= self.myStructXMD.prefix + self.myStructXMD._fields_[1][0]
            if (angleRot in labelNames):
                angleRotList=self.inMetadata[blockName].GetLoopItem(angleRot)
            else:
                angleRotList=[0.]*length
            
            angleTilt= self.myStructXMD.prefix + self.myStructXMD._fields_[2][0]
            if (angleTilt in labelNames):
                angleTiltList=self.inMetadata[blockName].GetLoopItem(angleTilt)
            else:
                angleTiltList=[0.]*length

            anglePsi= self.myStructXMD.prefix + self.myStructXMD._fields_[3][0]
            if (anglePsi in labelNames):
                anglePsiList=self.inMetadata[blockName].GetLoopItem(anglePsi)
            else:
                anglePsiList=[0.]*length

            shiftX= self.myStructXMD.prefix + self.myStructXMD._fields_[4][0]
            if (shiftX in labelNames):
                shiftXList=self.inMetadata[blockName].GetLoopItem(shiftX)
            else:
                shiftXList=[0.]*length
            
            shiftY= self.myStructXMD.prefix + self.myStructXMD._fields_[5][0]
            if (shiftX in labelNames):
                shiftYList=self.inMetadata[blockName].GetLoopItem(shiftY)
            else:
                shiftYList=[0.]*length
            
            shiftZ= self.myStructXMD.prefix + self.myStructXMD._fields_[6][0]
            if (shiftZ in labelNames):
                shiftZList=self.inMetadata[blockName].GetLoopItem(shiftZ)
            else:
                shiftZList=[0.]*length
            
            flip= self.myStructXMD.prefix + self.myStructXMD._fields_[7][0]
            if (flip in labelNames):
                flipList=self.inMetadata[blockName].GetLoopItem(flip)
            else:
                flipList=[0]*length
            
            scale= self.myStructXMD.prefix + self.myStructXMD._fields_[8][0]
            if (scale in labelNames):
                scaleList=self.inMetadata[blockName].GetLoopItem(scale)
            else:
                scaleList=[1.]*length
            
            enabled= self.myStructXMD.prefix + self.myStructXMD._fields_[9][0]
            if (enabled in labelNames):
                enabledList=self.inMetadata[blockName].GetLoopItem(enabled)
            else:
                enabledList=[1]*length
            
            fom= self.myStructXMD.prefix + self.myStructXMD._fields_[10][0]
            if (enabled in labelNames):
                fomList=self.inMetadata[blockName].GetLoopItem(fom)
            else:
                fomList=[1.]*length
            
            for i in range(len(fomList)):
                self.itemNameListXMD.append(ParticleAlignmentStructXmd
                                               (
                                                image=imageList[i],
                                                angleRot=float(angleRotList[i]),
                                                angleTilt=float(angleTiltList[i]),
                                                anglePsi=float(anglePsiList[i]),
                                                shiftX=float(shiftXList[i]),
                                                shiftY=float(shiftYList[i]),
                                                shiftZ=float(shiftZList[i]),
                                                flip=int(flipList[i]),
                                                scale=float(scaleList[i]),
                                                enable=int(enabledList[i]),
                                                fom=float(fomList[i])
                                               )
                                           )
                
#        for item in self.itemNameListXMD:
#            print item.image, item.angleRot, item.fom
#        exit(1)



    def EMX2XMD(self):
        

        ListEulerMatrices = []
        counterBlock=0
        for item in self.itemNameListEMX:
            _matrix = numpy.matrix([
                                    [ 
                                     float(item.transformation_matrix_1_1),
                                     float(item.transformation_matrix_1_2),
                                     float(item.transformation_matrix_1_3),
                                     float(item.transformation_matrix_offset_x)
                                    ],
                                    [ 
                                     float(item.transformation_matrix_2_1),
                                     float(item.transformation_matrix_2_2),
                                     float(item.transformation_matrix_2_3),
                                     float(item.transformation_matrix_offset_y)
                                    ],
                                    [ 
                                     float(item.transformation_matrix_3_1),
                                     float(item.transformation_matrix_3_2),
                                     float(item.transformation_matrix_3_3),
                                     float(item.transformation_matrix_offset_z)
                                    ],
                                    [0.,0.,0.,1.]
                                ])
#            ListEulerMatrices.append(_matrix)
            _scale, _shear, _angles, _trans, _persp = decompose_matrix(_matrix)
            if (   (math.fabs(_scale[0] - _scale[1]) > 0.001) 
                or (math.fabs(_scale[1] - _scale[2]) > 0.001) 
                or (math.fabs(_scale[0] - _scale[2]) > 0.001)
                ):
                print >> sys.stderr, "Scale is not equal in all directions"
                print >> sys.stderr, "XMIPP cannot handle different scales in each direction"
                print >> sys.stderr, "scale in x direction will be used"
                print >> sys.stderr, _scale[0],_scale[1],_scale[2]
            self.itemNameListXMD.append(ParticleAlignmentStructXmd
                                           (image     = item.url,
                                            angleRot  = _angles[0]*180./_pi,
                                            angleTilt = _angles[1]*180./_pi,
                                            anglePsi  = _angles[2]*180./_pi,
                                            shiftX    = _trans[0],
                                            shiftY    = _trans[1],
                                            shiftZ    = _trans[2],
                                            flip      = 0,
                                            scale     = _scale[0],
                                            enabled    = item.enable,
                                            fom       = item.FOM
                                           )
                                          )
        self.blockNameListXMD.append(BlockNamesXMD(
                                                       BlockName="aligment",
                                                       size=len(self.itemNameListXMD)
                                                   )
                                    )

#

#        for item in self.itemNameListXMD:
#             print item.Xcoor, item.Ycoor, item.BlockName
#        exit(1)
    def XMD2EMX(self):
        
        angle=[0]*3
        scale=[1]*3
        translate=[0]*3
        for item in self.itemNameListXMD:
            translate[0]=item.shiftX
            translate[1]=item.shiftY
            translate[2]=item.shiftZ
            
            angle[0] =  item.angleRot  * _pi /180.
            angle[1] =  item.angleTilt * _pi /180.
            angle[2] =  item.anglePsi  * _pi /180.

            scale[0]  = item.scale
            scale[1]  = item.scale
            scale[2]  = item.scale
            
            matrix = compose_matrix(scale=scale,
                                    angles=angle,
                                   translate=translate
                                   )

            self.itemNameListEMX.append(ParticleAlignmentStructEmx
                                           (
                                            url=item.image,
                                            transformation_matrix_1_1=matrix[0][0],
                                            transformation_matrix_1_2=matrix[0][1],
                                            transformation_matrix_1_3=matrix[0][2],
                                            transformation_matrix_offset_x=matrix[0][3],
                                            transformation_matrix_2_1=matrix[1][0],
                                            transformation_matrix_2_2=matrix[1][1],
                                            transformation_matrix_2_3=matrix[1][2],
                                            transformation_matrix_offset_y=matrix[1][3],
                                            transformation_matrix_3_1=matrix[2][0],
                                            transformation_matrix_3_2=matrix[2][1],
                                            transformation_matrix_3_3=matrix[2][2],
                                            transformation_matrix_offset_z=matrix[2][3],
                                            enable=item.enable,
                                            FOM=item.fom
                                           )
                                          )
        #here you may change the block name if needed
#        for item in self.blockNameListXMD:
        self.blockNameListEMX.append(BlockNamesEMX(
#                                                       BlockName="kk"+item.BlockName,
                                                       BlockName="aligment",
                                                       size=len(self.itemNameListEMX))
                                         )

#        for item in self.itemNameListEMX:
#             print item.url, item.transformation_matrix_1_1, item.FOM
#        exit(1)
        
    def XMD2startObject(self):
        myblock = CifFile.CifBlock()
        self.outMetadata[self.blockNameListXMD[0].BlockName] = myblock
        #get labels
        xmdStruct = ParticleAlignmentStructXmd()
        _image_     = xmdStruct.prefix + xmdStruct._fields_[0][0]
        _angleRot_  = xmdStruct.prefix + xmdStruct._fields_[1][0]
        _angleTilt_ = xmdStruct.prefix + xmdStruct._fields_[2][0]
        _anglePsi_  = xmdStruct.prefix + xmdStruct._fields_[3][0]
        _shiftX_    = xmdStruct.prefix + xmdStruct._fields_[4][0]
        _shiftY_    = xmdStruct.prefix + xmdStruct._fields_[5][0]
        _shiftZ_    = xmdStruct.prefix + xmdStruct._fields_[6][0]
        _flip_      = xmdStruct.prefix + xmdStruct._fields_[7][0]
        _scale_     = xmdStruct.prefix + xmdStruct._fields_[8][0]
        _enabled_   = xmdStruct.prefix + xmdStruct._fields_[9][0]
        _fom_       = xmdStruct.prefix + xmdStruct._fields_[10][0]
        #fill block
        imageList    =[]
        angleRotList =[]
        angleTiltList=[]
        anglePsiList =[]
        shiftXList   =[]
        shiftYList   =[]
        shiftZList   =[]
        flipList     =[]
        scaleList    =[]
        enabledList  =[]
        fomList      =[]
        for item in self.itemNameListXMD:
            imageList.append(str(item.image))
            angleRotList.append(str(item.angleRot))
            angleTiltList.append(str(item.angleTilt))
            anglePsiList.append(str(item.anglePsi))
            shiftXList.append(str(item.shiftX))
            shiftYList.append(str(item.shiftY))
            shiftZList.append(str(item.shiftZ))
            flipList.append(str(item.flip))
            scaleList.append(str(item.scale))
            enabledList.append(str(item.enabled))
            fomList.append(str(item.fom))
        self.outMetadata[self.blockNameListXMD[0].BlockName].AddCifItem(([[
                                                        _image_,    
                                                        _angleRot_, 
                                                        _angleTilt_,
                                                        _anglePsi_, 
                                                        _shiftX_,   
                                                        _shiftY_,   
                                                        _shiftZ_,   
                                                        _flip_,     
                                                        _scale_,    
                                                        _enabled_,  
                                                        _fom_
                                                     ]],
                                                      [[
                                                        imageList,    
                                                        angleRotList, 
                                                        angleTiltList,
                                                        anglePsiList, 
                                                        shiftXList,   
                                                        shiftYList,   
                                                        shiftZList,   
                                                        flipList,     
                                                        scaleList,    
                                                        enabledList,  
                                                        fomList      
                                                     ]]
                                                    ))

    def EMX2startObject(self):
        #create blocks in empty star file
        for block in self.blockNameListEMX:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = ParticleAlignmentStructEmx()

        _url_                            = xmdStruct.prefix + xmdStruct._fields_[0][0]
        _transformation_matrix_1_1_      = xmdStruct.prefix + xmdStruct._fields_[1][0]
        _transformation_matrix_1_2_      = xmdStruct.prefix + xmdStruct._fields_[2][0]
        _transformation_matrix_1_3_      = xmdStruct.prefix + xmdStruct._fields_[3][0]
        _transformation_matrix_offset_x_ = xmdStruct.prefix + xmdStruct._fields_[4][0]
        _transformation_matrix_2_1_      = xmdStruct.prefix + xmdStruct._fields_[5][0]
        _transformation_matrix_2_2_      = xmdStruct.prefix + xmdStruct._fields_[6][0]
        _transformation_matrix_2_3_      = xmdStruct.prefix + xmdStruct._fields_[7][0]
        _transformation_matrix_offset_y_ = xmdStruct.prefix + xmdStruct._fields_[8][0]
        _transformation_matrix_3_1_      = xmdStruct.prefix + xmdStruct._fields_[9][0]
        _transformation_matrix_3_2_      = xmdStruct.prefix + xmdStruct._fields_[10][0]
        _transformation_matrix_3_3_      = xmdStruct.prefix + xmdStruct._fields_[11][0]
        _transformation_matrix_offset_z_ = xmdStruct.prefix + xmdStruct._fields_[12][0]
        _enable_                         = xmdStruct.prefix + xmdStruct._fields_[13][0]
        _FOM_                            = xmdStruct.prefix + xmdStruct._fields_[14][0]

        
        #fill block
        for block in self.blockNameListEMX:
            urlList                            = []
            transformation_matrix_1_1List      = []
            transformation_matrix_1_2List      = []
            transformation_matrix_1_3List      = []
            transformation_matrix_offset_xList = []
            transformation_matrix_2_1List      = []
            transformation_matrix_2_2List      = []
            transformation_matrix_2_3List      = []
            transformation_matrix_offset_yList = []
            transformation_matrix_3_1List      = []
            transformation_matrix_3_2List      = []
            transformation_matrix_3_3List      = []
            transformation_matrix_offset_zList = []
            enableList                         = []
            FOMList                            = []
            for item in self.itemNameListEMX:
                    urlList.append(str(item.url))                            
                    transformation_matrix_1_1List.append(str(item.transformation_matrix_1_1))      
                    transformation_matrix_1_2List.append(str(item.transformation_matrix_1_2))      
                    transformation_matrix_1_3List.append(str(item.transformation_matrix_1_3))      
                    transformation_matrix_offset_xList.append(str(item.transformation_matrix_offset_x)) 
                    transformation_matrix_2_1List.append(str(item.transformation_matrix_2_1))      
                    transformation_matrix_2_2List.append(str(item.transformation_matrix_2_2))      
                    transformation_matrix_2_3List.append(str(item.transformation_matrix_2_3))      
                    transformation_matrix_offset_yList.append(str(item.transformation_matrix_offset_y)) 
                    transformation_matrix_3_1List.append(str(item.transformation_matrix_3_1))      
                    transformation_matrix_3_2List.append(str(item.transformation_matrix_3_2))      
                    transformation_matrix_3_3List.append(str(item.transformation_matrix_3_3))      
                    transformation_matrix_offset_zList.append(str(item.transformation_matrix_offset_z)) 
                    enableList.append(str(item.enable))                     
                    FOMList.append(str(item.FOM))                        
            self.outMetadata[block.BlockName].AddCifItem(([[
                                                        _url_,                           
                                                        _transformation_matrix_1_1_,     
                                                        _transformation_matrix_1_2_,     
                                                        _transformation_matrix_1_3_,     
                                                        _transformation_matrix_offset_x_,
                                                        _transformation_matrix_2_1_,     
                                                        _transformation_matrix_2_2_,     
                                                        _transformation_matrix_2_3_,     
                                                        _transformation_matrix_offset_y_,
                                                        _transformation_matrix_3_1_,     
                                                        _transformation_matrix_3_2_,     
                                                        _transformation_matrix_3_3_,     
                                                        _transformation_matrix_offset_z_,
                                                        _enable_,                        
                                                        _FOM_                                                                                       ]],
                                                          [[
                                                            urlList,                            
                                                            transformation_matrix_1_1List,      
                                                            transformation_matrix_1_2List,      
                                                            transformation_matrix_1_3List,      
                                                            transformation_matrix_offset_xList, 
                                                            transformation_matrix_2_1List,      
                                                            transformation_matrix_2_2List,      
                                                            transformation_matrix_2_3List,      
                                                            transformation_matrix_offset_yList, 
                                                            transformation_matrix_3_1List,      
                                                            transformation_matrix_3_2List,      
                                                            transformation_matrix_3_3List,      
                                                            transformation_matrix_offset_zList, 
                                                            enableList,                         
                                                            FOMList                   
                                                            ]]))
