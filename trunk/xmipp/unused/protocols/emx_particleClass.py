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
from emx_struct import ParticleClassStructEmx,\
                       ParticleClassStructXmd,\
                       BlockNamesEMX,\
                       BlockNamesXMD,\
                       prefix_particle
from protlib_emx import EmxBase
import CifFile
import StarFile

##########################################################################
#   Class Related to Particle Picking Conversion
##########################################################################
class ParticleClassConverter(EmxBase):    
    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
    
                       
    def readBlocksEMX(self):
        """ This function is pretty useless, think in removing it.
        It may be useful if only some blocks need to be processed
        """
        self.myStructEMX = ParticleClassStructEmx()
        label = prefix_particle + self.myStructEMX._fields_[1][0]
        #label = self.myStructEMX.prefix + self.myStructEMX._fields_[0][0]
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListEMX.append(BlockNamesEMX(BlockName=blockName,
                                                    size=len(_auxList)))
#        for item in self.blockNameListEMX:
#            print item.BlockName, item.size
#        exit(1)
    def readBlocksXMD(self):
        """ This function is pretty useless, think in removing it.
        It may be useful if only some blocks need to be processed
        """
        self.myStructXMD = ParticleClassStructXmd()
        label = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
        
        for blockName in self.inMetadata.keys():
            _auxList=self.inMetadata[blockName].GetLoopItem(label)
            self.blockNameListXMD.append(BlockNamesXMD(BlockName=blockName,
                                                    size=len(_auxList)))
#        for item in self.blockNameListXMD:
#            print item.BlockName, item.size
#        exit(1)

    def startObject2EMX(self):
        classId=[]
        url = prefix_particle + self.myStructEMX._fields_[1][0]
        for item in self.blockNameListEMX:
            classId=item.BlockName
            _url_=self.inMetadata[item.BlockName].GetLoopItem(url)
            for i in range(len(_url_)):
                self.itemNameListEMX.append(ParticleClassStructEmx
                                           (
                                            id=classId,
                                            url=_url_[i]
                                           )
                                        )
#        for item in self.itemNameListEMX:
#            print "item.id, item.url
#        exit(1)
    def startObject2XMD(self):
        for item in self.blockNameListXMD:
            blockName = item.BlockName   

            image = self.myStructXMD.prefix + self.myStructXMD._fields_[0][0]
            _image_=self.inMetadata[blockName].GetLoopItem(image)
            ref = self.myStructXMD.prefix + self.myStructXMD._fields_[1][0]
            _ref_=self.inMetadata[blockName].GetLoopItem(ref)
            
            for i in range(len(_image_)):
                self.itemNameListXMD.append(ParticleClassStructXmd
                                           (
                                            image=_image_[i],
                                            ref=int(_ref_[i])
                                           )
                                        )
#        for item in self.itemNameListXMD:
#            print item.image, item.ref
#        exit(1)

    def EMX2XMD(self):
        
        for item in self.itemNameListEMX:
            self.itemNameListXMD.append(ParticleClassStructXmd
                                           (
                                            image=item.url,
                                            #this is overkill but convert strings to ints
                                            ref=hash(item.id) & 0xffffffff,
                                           )
                                          )
        #here you may change the block name if needed
        #we need a blockname with the micrograph name, that is, MicName
        #
        self.blockNameListXMD.append(BlockNamesXMD(
                                       BlockName="class",
                                       size=len(self.itemNameListXMD))
                                         )

#        for item in self.itemNameListXMD:
#             print item.image, item.ref
#        for item in self.blockNameListXMD:
#             print item.BlockName, item.size
#        exit(1)

    def XMD2EMX(self):
        block=[]
        for item in self.itemNameListXMD:
            self.itemNameListEMX.append(ParticleClassStructEmx
                                           (
                                            id=str(item.ref),
                                            url=item.image
                                           )
                                          )
            block.append(str(item.ref))
        #here you may change the block name if needed
        blockDifferent=list(set(block))
        for item in blockDifferent:
            self.blockNameListEMX.append(BlockNamesEMX(
                                                       BlockName=item,
                                                       size=block.count(item))
                                         )

#        for item in self.itemNameListEMX:
#             print item.id, item.url
#        for item in self.blockNameListEMX:
#             print item.BlockName, item.size
#        exit(1)
        
    def XMD2startObject(self):
        #create blocks in empty star file
        for block in self.blockNameListXMD:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = ParticleClassStructXmd()
        _image_ = xmdStruct.prefix + xmdStruct._fields_[0][0]
        _ref_   = xmdStruct.prefix + xmdStruct._fields_[1][0]
        #fill block
        for block in self.blockNameListXMD:
            imageList =[]
            refList =[]
            for item in self.itemNameListXMD:
                imageList.append(str(item.image))
                refList.append(str(item.ref))
            self.outMetadata[block.BlockName].AddCifItem(([[
                                                _image_,
                                                _ref_
                                                ]],
                                                [[
                                                  imageList,
                                                  refList
                                                  ]]))

    def EMX2startObject(self):
        #create blocks in empty star file
        #block names in emx are irrelevant but I will use the same than in xmipp
        #remember than blockname in xmipp is micrograph_url
        for block in self.blockNameListEMX:
            myblock = CifFile.CifBlock()
            self.outMetadata[block.BlockName] = myblock
        #get labels
        xmdStruct = ParticleClassStructEmx()
        emdStruct = ParticleClassStructXmd()
        _id_ = xmdStruct.prefix + emdStruct._fields_[0][0]
        _url_   = prefix_particle + emdStruct._fields_[1][0]
        #fill block
        for block in self.blockNameListEMX:
            idList =[]
            urlList =[]
            self.outMetadata[block.BlockName][_url_] = block.BlockName
            for item in self.itemNameListEMX:
                if block.BlockName == item.id:
                    urlList.append(item.url)
            self.outMetadata[block.BlockName].AddCifItem(([[
                                                            _url_
                                                            ]],
                                                          [[
                                                            urlList
                                                            ]]
                                                          ))

