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



###############################################
# The following code may be used to launch de program
# This program is available through xmipp package as:
# xmipp_emx_convert
###############################################
#from emx_coordinates    import ParticlePickingConverter
#from emx_ctf_Micrograph import CtfMicrographConverter
#from emx_ctf_Particle   import CtfParticleConverter
#from emx_particleAlignment import ParticleAlignmentConverter
#if __name__ == '__main__':
#    
#    def command_line_options():
#        """ add command line options here"""
#        import optparse
#        _usage = "usage: %prog [options] Example:   %prog -i input.xmd -o out.emx"
#        parser = optparse.OptionParser(_usage)        
#        parser.add_option("-i", "--input_filename", dest="inputFn",
#                          default="/dev/stdin", type="string",
#                          help="EMS or XMIPP Metadata file")                           
#        parser.add_option("-o", "--output_filename", dest="outputFn",
#                          default="/dev/stdout", type="string",
#                          help="XMIPP or EMX Metadata file")   
#        parser.add_option("-t", "--conversion_type", dest="type",
#                          default="coordinates", type="string",
#                          help="Possible types of conversion: coordinates, alignment, class, ctfMicrograph, ctfParticle")                          
#       
#        (options, args) = parser.parse_args()
#        return(options.inputFn,options.outputFn,options.type)
#        
#    inputFn,outputFn,convType=command_line_options()
#    if convType == 'coordinates':
#        ParticlePickingConverter(inputFn, outputFn).run()  
#    elif convType == 'alignment':
#        ParticleAlignmentConverter(inputFn, outputFn).run()  
#    elif convType == 'class':
#        ParticleClassConverter(inputFn, outputFn).run()  
#    elif convType == 'ctfMicrograph':
#        CtfMicrographConverter(inputFn, outputFn).run()  
#    elif convType == 'ctfParticle':
#        CtfParticleConverter(inputFn, outputFn).run()  
#    else:
#        print >> sys.stderr,  "ERROR: Wrong mode: ", convType
#        exit(0)
