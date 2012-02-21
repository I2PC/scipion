#!/usr/bin/env xmipp_python
#cat Test/particlePicking.xmd | ./batch_emx_coordinates.py
#cat Test/particlePicking.emx | ./batch_emx_coordinates.py
import CifFile
from protlib_emx import emxBase

class convertParticlePickingClass(emxBase):    
    needed_itemsXMIPP =  (
          "_Xcoor",
          "_Ycoor"
          )
    needed_itemsEMX = (
          "_emx_particle.coordinate_x",
          "_emx_particle.coordinate_y"
          )

    
    def __init__(self):

        self.command_line_options()
        emx2xmipp=self.checkVersion()
        if emx2xmipp:
            self.convertAllBlocksEMX2XMIPP()
            self.saveFileEMX2XMIPP()
        else:
            self.convertAllBlocksXMIPP2EMX()
            self.saveFileXMIPP2EMX()                     
    def convertAllBlocksEMX2XMIPP(self):
        """loop over the blocks and write them"""
        for blockName in self.inMetadata.keys():
            #create output block
            myblock = CifFile.CifBlock()
            #read micrograph.url field
            micrographName = self.inMetadata[blockName]['_emx_micrograph.url'] 
            self.outMetadata[micrographName] = myblock
            self.cbOut = self.outMetadata[micrographName]#alias
            self.convertLoopEMX2XMIPP(blockName)
            
    def convertAllBlocksXMIPP2EMX(self):
        """loop over the blocks and write them"""
        for micrographName in self.inMetadata.keys():
            #create output block
            myblock = CifFile.CifBlock()
            self.outMetadata[micrographName] = myblock
            self.cbOut = self.outMetadata[micrographName]
            self.cbOut['_emx_micrograph.url'] = micrographName
            self.createDataHeaderXMIPP2EMX(micrographName)
            self.convertLoopXMIPP2EMX(micrographName)
            
    def createDataHeaderXMIPP2EMX(self,micrographName):
        """Data header is the xmipp data block name with the right label"""
        self.cbOut['_emx_micrograph.url'] = micrographName

    def convertLoopEMX2XMIPP(self,micrographName):
    
        _auxXList=self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[0])
        _auxYList=self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[1])
        _XList =[]
        _YList =[]
        for _item in range (len(_auxXList)): 
            _XList.append(str(int(round(float(_auxXList[_item]),0))))
            _YList.append(str(int(round(float(_auxYList[_item]),0))))
            
        self.cbOut.AddCifItem(([self.needed_itemsXMIPP],[[_XList,_YList]])) 
        
    def convertLoopXMIPP2EMX(self,micrographName):
        loopitems = self.inMetadata[micrographName].GetLoop((self.needed_itemsXMIPP[0]))  #get item names and values
        self.cbOut.AddCifItem(([self.needed_itemsEMX],\
              [[loopitems[self.needed_itemsXMIPP[0]],loopitems[self.needed_itemsXMIPP[1]]]]))
        
if __name__ == '__main__':

    convertParticlePickingClass()
