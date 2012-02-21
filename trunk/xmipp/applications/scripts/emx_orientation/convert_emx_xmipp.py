#!/usr/bin/env python
import CifFile
from transformations import *
smallNumber=0.00001
""" Test 4x4 matrix"""

class convertParticleAligmentClass:    
    needed_itemsXMIPP = [[
         "_image",
         "_angleRot",
         "_angleTilt",
         "_anglePsi",
         "_shiftX",
         "_shiftY",
         "_shiftZ",
         "_flip",
         "_scale",
         "_enable",
         "_fom"
        ]]
    needed_itemsEMX = [[
            "_emx_particle.url",
            "_emx_particle.transformation_matrix_1_1",
            "_emx_particle.transformation_matrix_1_2",
            "_emx_particle.transformation_matrix_1_3",
            "_emx_particle.transformation_matrix_offset_x",
            "_emx_particle.transformation_matrix_2_1",
            "_emx_particle.transformation_matrix_2_2",
            "_emx_particle.transformation_matrix_2_3",
            "_emx_particle.transformation_matrix_offset_y",
            "_emx_particle.transformation_matrix_3_1",
            "_emx_particle.transformation_matrix_3_2",
            "_emx_particle.transformation_matrix_3_3",
            "_emx_particle.transformation_matrix_offset_z",
            "_emx_particle.enable",
            "_emx_particle.FOM"
                   ]]
    xmippStartVersion = 'XMIPP_STAR_1'
    emxVersion = 'EMX1.0'
    contactMail = 'xmipp@cnb.csic.es'
    def __init__(self,
                inputFileName, outputFileName):
        self.inputFileName = inputFileName
        self.outputFileName = outputFileName
        self.checkVersion()
        self.inMetadata = CifFile.CifFile(self.inputFileName)
        self.outMetadata = CifFile.CifFile()
        self.convertAllBlocks()
        self.saveFile()
        
    def checkVersion(self):
        """read first 7 characters. If different from #EMX1.0
           then abort """
        import os
        if not os.path.exists(self.inputFileName):
            print "File: ", inputFileName, "does not exists."
            exit(0)
        #set buffering size to 0 otherwise stdin is read in a buffer
        #and it is inavailable dor the next open
        fin = open(self.inputFileName, "r", 0)
        firstLine = fin.readline()
        result = firstLine.find(self.emxVersion)
        if (result == -1):
            print "Error: EMX Metadata Files should contain the string ", \
                   self.emxVersion, \
                    " in the first line.  Exiting program."
            exit(1)
        fin.close()
            
                     
    def convertAllBlocks(self):
        """loop over the blocks and write them"""
        for blockName in self.inMetadata.keys():
            #create output block
            myblock = CifFile.CifBlock()
            #read micrograph.url field
            #micrographName = self.inMetadata[blockName]['_emx_micrograph.url'] 
            self.outMetadata[blockName] = myblock
            self.cbOut = self.outMetadata[blockName]
            #self.createDataHeader(blockName)
            self.convertLoop(blockName)

    def createDataHeader(self, micrographName):
        """Data header is the xmipp data block name with the right label"""
        self.cbOut['_emx_micrograph.url'] = micrographName    
    
    def convertLoop(self, blockName):
        
        #get item names and values
        loopitems = self.inMetadata[blockName].GetLoop(self.needed_itemsEMX[0][0])  
        #get block list
        ListBlockKeys = loopitems.keys()
        #get loop length
        blockLenght = len (self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][1]))
        #init output, check if 2x2 or 3x3 rotation
        
        #image is compulsory 
        _imageUrl   = self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][0])
        #store matrices here
        ListEulerMatrices = []
        for loopitem in range(blockLenght):
            ListEulerMatrices.append(numpy.identity(4))
        #fill shifts
        if (self.needed_itemsEMX[0][4]  in ListBlockKeys):
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][4])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][0][3]=aux[loopitem]
        
        if (self.needed_itemsEMX[0][8]  in ListBlockKeys):
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][8])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][1][3]=aux[loopitem]
        
        if (self.needed_itemsEMX[0][12]  in ListBlockKeys):
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][12])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][3]=aux[loopitem]
        # fill 3D rotations
        _2Dmatrix=True
        if (self.needed_itemsEMX[0][3]  in ListBlockKeys):
            _2Dmatrix=False
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][3])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][0][2]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][7])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][1][2]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][9])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][0]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][10])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][1]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][11])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][2]=aux[loopitem]
        #2D rotations
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][1])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][0][0]=aux[loopitem]
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][2])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][0][1]=aux[loopitem]
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][5])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][1][0]=aux[loopitem]
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][6])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][1][1]=aux[loopitem]

        #handle enable
        _enable = [1]*blockLenght
        if (self.needed_itemsEMX[0][13]  in ListBlockKeys):
            _enable= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][13])#enable
        #handle FOM
        _fom = [1.]*blockLenght
        if (self.needed_itemsEMX[0][14]  in ListBlockKeys):
            _fom= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0][14])#enable
        
        #availabel in _imageUrl
        #_image     = []
        _angleRot  = []
        _angleTilt = []
        _anglePsi  = []
        _shiftX    = []
        _shiftY    = []
        _shiftZ    = []
        _flip      = [0]*blockLenght
        _scale     = []

        for loopitem in range(blockLenght):
            #create 4x4 matrix
            scale, shear, angles, trans, persp = decompose_matrix(ListEulerMatrices[loopitem])
            #images should not present shear or persp
            if (math.fabs(numpy.linalg.norm(shear))>smallNumber):
                print "The input matrix ", ListEulerMatrices[loopitem],\
                "presents shear. This is not supported"
                exit(1)
            if (math.fabs(numpy.linalg.norm(persp)-1.)>smallNumber):
                print "The input matrix ", ListEulerMatrices[loopitem],\
                "presents perpective. This is not supported"
                exit(1)

            _angleRot.append ("%0.4f" % (angles[0]* 180./math.pi))
            _angleTilt.append("%0.4f" % (angles[1]* 180./math.pi))
            _anglePsi.append ("%0.4f" % (angles[2]* 180./math.pi))
            _shiftX.append ("%0.4f" % (trans[0]))
            _shiftY.append ("%0.4f" % (trans[1]))
            _shiftZ.append ("%0.4f" % (trans[2]))
            #flip only in 2D
            if(_2Dmatrix):
                if (numpy.linalg.det()< 0):
                    _flip[loopitem] = 1
            #check if different in each direction
            scaleAverage = (scale[0]+scale[1]+scale[2])/3.
            if(math.fabs(scaleAverage-scale[0])>smallNumber):
                print "Reading Iamge:", _imageUrl[loopitem], 
                print "Different scale factor in each axis.", scale, "This is not supported", 
                print "scale along x axis will be assigned"
            _scale.append ("%0.4f" % (scaleAverage))
            
        self.cbOut.AddCifItem((self.needed_itemsXMIPP, [
                               [
                                _imageUrl,
                                _angleRot,
                                _angleTilt,
                                _anglePsi,
                                _shiftX,
                                _shiftY,
                                _shiftZ,
                                _flip,
                                _scale,
                                _enable,
                                _fom
                                ]
                              ])) 
    def saveFile(self):
        comment = "# XMIPP_STAR_1 *"
        comment += "\n##########################################################################"
        comment += "\n#  Converted from " + self.emxVersion + " to " + self.xmippStartVersion
        comment += "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################" 
        outfile = open(self.outputFileName, "w")
        outfile.write(self.outMetadata.WriteOut(comment=comment, _email=self.contactMail))
        

        
if __name__ == '__main__':

    def command_line_options():
        """ add command line options here"""
        import optparse
        _usage = "usage: %prog [options] Example:   %prog -i input.xmd -o out.emx"
        parser = optparse.OptionParser(_usage)        
        parser.add_option("-i", "--inputFileName", dest="inputFileName",
                          default="/dev/stdin", type="string",
                          help="EMS Metadata file (particle aligment)")                           
        parser.add_option("-o", "--outputFileName", dest="outputFileName",
                          default="/dev/stdout", type="string",
                          help="XMIPP   Metadata file (particle aligment)")                           
       
        (options, args) = parser.parse_args()
        return(options.inputFileName, options.outputFileName)
    
    inputFileName, outputFileName = command_line_options()
    instance = convertParticleAligmentClass(inputFileName, outputFileName)
