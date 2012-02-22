import CifFile
class emxBase:
    xmippStartVersion = 'XMIPP_STAR_1'
    emxVersion        = 'EMX1.0'
    contactMail       = 'xmipp@cnb.csic.es'
    smallNumber       = 0.00001
    
#    def __init__(self, runWithoutArgs=False):
#        a=0
    def command_line_options(self):
        """ add command line options here"""
        import optparse
        _usage = "usage: %prog [options] Example:   %prog -i input.xmd -o out.emx"
        parser = optparse.OptionParser(_usage)        
        parser.add_option("-i", "--inputFileName", dest="inputFileName",
                          default="/dev/stdin", type="string",
                          help="EMS or XMIPP Metadata file")                           
        parser.add_option("-o", "--outputFileName", dest="outputFileName",
                          default="/dev/stdout", type="string",
                          help="XMIPP or EMX Metadata file")                           
       
        (options, args) = parser.parse_args()
        
        self.inputFileName  = options.inputFileName
        self.outputFileName = options.outputFileName

    
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
