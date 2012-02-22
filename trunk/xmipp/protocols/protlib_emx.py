import CifFile

class EmxBase:
    xmippStartVersion = 'XMIPP_STAR_1'
    emxVersion        = 'EMX1.0'
    contactMail       = 'xmipp@cnb.csic.es'
    smallNumber       = 0.00001
    
#    def __init__(self, runWithoutArgs=False):
#        a=0

    
    def checkVersion(self):
        """read first line: check IPP magic word. Abort otherwise"""
        import os
        if not os.path.exists(self.inputFileName):
            print "File: ", self.inputFileName, "does not exists."
            exit(0)
        #set buffering size to 0 otherwise stdin is read in a buffer
        #and it is unavailable for the next open
        fin = open(self.inputFileName, "r",0)
        firstLine = fin.readline()
        print "filename first line", self.inputFileName, firstLine
        
        result = firstLine.find(self.emxVersion)
        if result != -1:
            fin.close()
            #is EMX
            return (True)
        
        result = firstLine.find(self.xmippStartVersion)
        if result != -1:
            fin.close()
            return (False) # isEMX
        
        print "Error: Metadata Files should contain the string ",\
               self.emxVersion, "or", self.xmippStartVersion, \
            " in the first line.  Exiting program.\n", \
            "First line: ", firstLine
        exit(1)
    
    def saveFileEMX2XMIPP(self):
        comment  =   "# XMIPP_STAR_1 *"
        comment += "\n##########################################################################"         
        comment +=  "\n#  Converted from " + self.emxVersion + " to " + self.xmippStartVersion
        comment +=  "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################" 
        outfile = open(self.outputFileName,"w")
        outfile.write(self.outMetadata.WriteOut(comment=comment,_email=self.contactMail))
        
        
    def saveFileXMIPP2EMX(self):
        outfile = open(self.outputFileName,"w")
        comment =  "#  Converted from " + self.xmippStartVersion + " to " + self.emxVersion
        comment += "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################"
        outfile.write(self.outMetadata.WriteOut(_add=True,comment=comment,_email=self.contactMail))


class ParticlePickingConverter(EmxBase):    
    needed_itemsXMIPP =  (
          "_Xcoor",
          "_Ycoor"
          )
    needed_itemsEMX = (
          "_emx_particle.coordinate_x",
          "_emx_particle.coordinate_y"
          )

    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
        #next two lines should go after checkVersion
    
    def run(self):
        emx2xmipp = self.checkVersion()
        self.inMetadata     = CifFile.CifFile(self.inputFileName)
        self.outMetadata    = CifFile.CifFile()

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
                          help="Possible types of conversion: coordinates, alignment, class, ctf")                          
       
        (options, args) = parser.parse_args()
        return(options.inputFn,options.outputFn,options.type)
        
    inputFn,outputFn,type=command_line_options()
    if type == 'coordinates':
        ParticlePickingConverter(inputFn, outputFn).run()  
    elif type == 'alignment':
        ParticleAlignmentConverter(inputFn, outputFn).run()  
    elif type == 'class':
        ParticleClassConverter(inputFn, outputFn).run()  
    elif type == 'ctf':
        CtfConverter(inputFn, outputFn).run()  
    else:
        print "ERROR: Wrong mode: ", type
        exit(0)
    
    
    
        
