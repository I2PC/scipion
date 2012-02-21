#!/usr/bin/env xmipp_python
import CifFile
from protlib_emx import *
class convertParticlePickingClass:    
    needed_itemsXMIPP =  [
          "_Xcoor",
          "_Ycoor"
          ]
    needed_itemsEMX = [
          "_emx_particle.coordinate_x",
          "_emx_particle.coordinate_y"
          ]

#    needed_itemsXMIPP = (
#          "_Xcoor",
#          "_Ycoor"
#      )
#    needed_itemsEMX = [[
#          "_emx_particle.coordinate_x",
#          "_emx_particle.coordinate_y"
#      ]]

    
    def __init__(self,
                inputFileName, outputFileName):

        self.inputFileName  = inputFileName
        self.outputFileName = outputFileName
        emx2xmipp=checkVersion(self.inputFileName)
        self.inMetadata     = CifFile.CifFile(self.inputFileName)
        self.outMetadata    = CifFile.CifFile()
        if emx2xmipp:
            self.convertAllBlocksEMX2XMIPP()
        else:
            self.convertAllBlocksXMIPP2EMX()
        self.saveFile()
                     
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
            self.convertLoop(micrographName)

    def convertLoopEMX2XMIPP(self,micrographName):
    
        #loopitems = self.inMetadata[micrographName].GetLoop(self.needed_itemsEMX[0][0])  #get item names and values
        # It is possible to create the metadata without loops but...
        # since xmipp can only use integer values we need to loop and cast
        # if no casting is need the next to lines are the way to go
        #self.cbOut.AddCifItem((self.needed_itemsXMIPP,\
        #      [[loopitems[self.needed_itemsEMX[0][0]],loopitems[self.needed_itemsEMX[0][1]]]])) 
        _auxXList=self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[0])
        _auxYList=self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[1])
        _XList =[]
        _YList =[]
        for _item in range (len(_auxXList)): 
            _XList.append(str(int(round(float(_auxXList[_item]),0))))
            _YList.append(str(int(round(float(_auxYList[_item]),0))))
            
        self.cbOut.AddCifItem(([self.needed_itemsXMIPP],[[_XList,_YList]])) 
        
    def convertLoopXMIPP2EMX(self,micrographName):
        loopitems = self.inMetadata[micrographName].GetLoop(self.needed_itemsXMIPP)  #get item names and values
        self.cbOut.AddCifItem(([self.needed_itemsEMX],\
              [[loopitems[[self.needed_itemsXMIPP[0]]],loopitems[[self.needed_itemsXMIPP[1]]]]]))

    def saveFile(self):
        comment  =   "# XMIPP_STAR_1 *"
        comment += "\n##########################################################################"         
        comment +=  "\n#  Converted from " + emxVersion + " to " + xmippStartVersion
        comment +=  "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################" 
        outfile = open(self.outputFileName,"w")
        outfile.write(self.outMetadata.WriteOut(comment=comment,_email=contactMail))
        

        
if __name__ == '__main__':

    
    inputFileName,outputFileName  = command_line_options()
    instance = convertParticlePickingClass(inputFileName,outputFileName)
