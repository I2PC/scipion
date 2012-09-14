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
from emx_struct import CtfMicrographStructEmx,\
                       CtfMicrographStructXmd,\
                       BlockNamesEMX,\
                       BlockNamesXMD,\
                       prefix_micrograph
from protlib_emx import EmxBase
import CifFile
import StarFile

###########################################################################
##   Class Related to CTF conversion: defocus per micrograph
###########################################################################
class CtfMicrographConverter(EmxBase):    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
    
    def run(self):
        EmxBase.run(self)
        if self.emx2xmipp:
            #each ctf in a file
            myStructXmd = CtfMicrographStructXmd()
            label = myStructXmd.prefix + myStructXmd._fields_[0][0]
            self.saveFilePerLineXMD(label)

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
        for item in self.blockNameListEMX:
            blockName = item.BlockName
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
        for item in self.blockNameListXMD:
            blockName = item.BlockName
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
#        exit(1)
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
#             print item.magnification, item.defocusU, item.astigmatism_angle
#        exit(1)

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
        _voltage_ = xmdStruct.prefix + xmdStruct._fields_[6][0]
        _Cs_ = xmdStruct.prefix + xmdStruct._fields_[7][0]
        _amplitude_contrast_= xmdStruct.prefix + xmdStruct._fields_[8][0]

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

