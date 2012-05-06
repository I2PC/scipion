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
from emx_struct import ParticlePickingStructEmx,\
                       ParticlePickingStructXmd,\
                       BlockNamesEMX,\
                       BlockNamesXMD,\
                       prefix_micrograph
from protlib_emx import EmxBase
import CifFile
import StarFile

##########################################################################
#   Class Related to Particle Picking Conversion
##########################################################################
class ParticlePickingConverter(EmxBase):    
    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
    
                       
    def readBlocksEMX(self):
        """ This function is pretty useless, think in removing it.
        It may be useful if only some blocks need to be processed
        """
        self.myStructEMX = ParticlePickingStructEmx()
        label = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListEMX.append(BlockNamesEMX(BlockName=blockName,
                                                    size=len(_auxList)))
#        for item in self.blockNameListEMX:
#            print item.BlockName, item.size

    def readBlocksXMD(self):
        """ This function is pretty useless, think in removing it.
        It may be useful if only some blocks need to be processed
        """
        self.myStructXMD = ParticlePickingStructXmd()
        label = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListXMD.append(BlockNamesXMD(BlockName=blockName,
                                                    size=len(_auxList)))
#        for item in self.blockNameListXMD:
#            print item.BlockName, item.size

    def startObject2EMX(self):
        for item in self.blockNameListEMX:
            blockName = item.BlockName
           #for i in range(len(self.myStruct._fields_)-1):
            coordinate_x = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
            _auxCoordenateX=self.inMetadata[blockName].GetLoopItem(coordinate_x)

            coordinate_y = self.myStructEMX.prefix + self.myStructEMX._fields_[1][0]
            _auxCoordenateY=self.inMetadata[blockName].GetLoopItem(coordinate_y)

            micrograph = emx_struct.prefix_micrograph + self.myStructEMX._fields_[2][0]
            _micrograph=self.inMetadata[blockName][micrograph]
            for i in range(len(_auxCoordenateX)):
                self.itemNameListEMX.append(ParticlePickingStructEmx
                                           (
                                            coordinate_x=float(_auxCoordenateX[i]),
                                            coordinate_y=float(_auxCoordenateY[i]),
                                            url = _micrograph
                                           )
                                        )
#        for item in self.itemNameListEMX:
#            print item.coordinate_x, item.coordinate_y, item.url

    def startObject2XMD(self):
        for item in self.blockNameListXMD:
            blockName = item.BlockName   

            self.Xcoor = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
            _auxCoordenateX=self.inMetadata[blockName].GetLoopItem(self.Xcoor)
            self.Ycoor = self.myStructXMD.prefix + self.myStructXMD._fields_[1][0]
            _auxCoordenateY=self.inMetadata[blockName].GetLoopItem(self.Ycoor)
            
            for i in range(len(_auxCoordenateX)):
                self.itemNameListXMD.append(ParticlePickingStructXmd
                                           (
                                            Xcoor=int(_auxCoordenateX[i]),
                                            Ycoor=int(_auxCoordenateY[i]),
                                            MicName = blockName
                                           )
                                        )
#        for item in self.itemNameListXMD:
#            print item.Xcoor, item.Ycoor, item.MicName
#        exit(1)

    def EMX2XMD(self):
        
        for item in self.itemNameListEMX:
            self.itemNameListXMD.append(ParticlePickingStructXmd
                                           (
                                            Xcoor=int(round(item.coordinate_x)),
                                            Ycoor=int(round(item.coordinate_y)),
                                            #BlockName = "kk"+item.BlockName
                                            MicName = item.url
                                           )
                                          )
        #here you may change the block name if needed
        #we need a blockname with the micrograph name, that is, MicName
        #
        micrograph = emx_struct.prefix_micrograph + self.myStructEMX._fields_[2][0]
        for item in self.blockNameListEMX:
            _micrograph=self.inMetadata[item.BlockName][micrograph]
            self.blockNameListXMD.append(BlockNamesXMD(
#                                                       BlockName="kk"+item.BlockName,
                                                       BlockName=_micrograph,
                                                       size=item.size)
                                         )

#        for item in self.itemNameListXMD:
#             print item.Xcoor, item.Ycoor, item.MicName
#        for item in self.blockNameListXMD:
#             print item.BlockName, item.size
#        exit(1)

    def XMD2EMX(self):
        
        for item in self.itemNameListXMD:
            self.itemNameListEMX.append(ParticlePickingStructEmx
                                           (
                                            coordinate_x=float(item.Xcoor),
                                            coordinate_y=float(item.Ycoor),
                                            #BlockName = "kk"+item.BlockName
                                            url = item.MicName
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
#             print item.coordinate_x, item.coordinate_y, item.url
#        exit(1)
        
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
                if block.BlockName == item.MicName:
                    XcoorList.append(str(item.Xcoor))
                    YcoorList.append(str(item.Ycoor))
            self.outMetadata[block.BlockName].AddCifItem(([[_Xcorr_,_Ycoor_]],[[XcoorList,YcoorList]]))

    def EMX2startObject(self):
        #create blocks in empty star file
        #block names in emx are irrelevant but I will use the same than in xmipp
        #remember than blockname in xmipp is micrograph_url
        for block in self.blockNameListEMX:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = ParticlePickingStructEmx()
        emdStruct = ParticlePickingStructXmd()
        _Xcorr_ = xmdStruct.prefix + xmdStruct._fields_[0][0]
        _Ycoor_ = xmdStruct.prefix + xmdStruct._fields_[1][0]
        _url_   = prefix_micrograph + emdStruct._fields_[2][0]
        #fill block
        for block in self.blockNameListEMX:
            XcoorList =[]
            YcoorList =[]
            self.outMetadata[block.BlockName][_url_] = block.BlockName
            for item in self.itemNameListEMX:
                if block.BlockName == item.url:
                    XcoorList.append(str(item.coordinate_x))
                    YcoorList.append(str(item.coordinate_y))
            self.outMetadata[block.BlockName].AddCifItem(([[_Xcorr_,_Ycoor_]],[[XcoorList,YcoorList]]))

