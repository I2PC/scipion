#!/usr/bin/env python
import CifFile
from transformations import *

class convertParticlePickingClass:    
    needed_itemsXMIPP = (
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
        )

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
    emxVersion        = 'EMX1.0'
    contactMail       = 'xmipp@cnb.csic.es'

    def __init__(self,
                inputFileName, outputFileName):

        self.inputFileName  = inputFileName
        self.outputFileName = outputFileName
        self.checkVersion()
        self.inMetadata     = CifFile.CifFile(self.inputFileName)
        self.outMetadata    = CifFile.CifFile()
        self.convertAllBlocks()
        self.saveFile()
    
    def checkVersion(self):
        """read first 16 characters. If different from # XMIPP_STAR_1 *
           then abort """
        import os
        if not os.path.exists(self.inputFileName):
            print "File: ", inputFileName, "does not exists."
            exit(0)
        #set buffering size to 0 otherwise stdin is read in a buffer
        #and it is unavailable for the next open
        fin = open(self.inputFileName, "r",0)
        firstLine = fin.readline()
        result = firstLine.find(self.xmippStartVersion)
        if (result==-1):
            print "Error: Xmipp Metadata Files should contain the string ",\
                   self.xmippStartVersion,\
                " in the first line.  Exiting program."
            exit(1)
        fin.close()
        
                 
    def convertAllBlocks(self):
        """loop over the blocks and write them"""
        for blockName in self.inMetadata.keys():
            #create output block
            myblock = CifFile.CifBlock()
            self.outMetadata[blockName] = myblock
            self.cbOut = self.outMetadata[blockName]
            #self.createDataHeader(blockName)
            self.convertLoop(blockName)
         
    #def createDataHeader(self,blockName):
    #    """Data header is the xmipp data block name with the right label"""
    #    self.cbOut['_emx_micrograph.url'] = blockName    
    
    def convertLoop(self,blockName):
        #get item names and values
        loopitems = self.inMetadata[blockName].GetLoop(self.needed_itemsXMIPP[0])  
        #get block list
        ListBlockKeys = loopitems.keys()
        #get loop length
        blockLenght = len (self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[0]))        
        #image is compulsory 
        _imageUrl   = self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[0])
        #store intermediate results here
        shift=[]
        angle=[]
        scale=[]

        _x=[]
        _y=[]
        _z=[]
        for loopitem in range(blockLenght):
            _x.append(0)
            _y.append(0)
            _z.append(0)
        
        #fill shifts
        if (self.needed_itemsXMIPP[4]  in ListBlockKeys):
            _x= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[4])
        if (self.needed_itemsXMIPP[5]  in ListBlockKeys):
            _y= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[5])
        if (self.needed_itemsXMIPP[6]  in ListBlockKeys):
            _z= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[6])
        for loopitem in range(blockLenght):
            shift.append((_x[loopitem],_y[loopitem],_z[loopitem]))

        _rot=[]
        _tilt=[]
        _psi=[]
        
        for loopitem in range(blockLenght):
            _rot.append(0)
            _tilt.append(0)
            _psi.append(0)
        #fill rotations
        if (self.needed_itemsXMIPP[1]  in ListBlockKeys):
            _rot= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[1])
        if (self.needed_itemsXMIPP[2]  in ListBlockKeys):
            _tilt= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[2])
        if (self.needed_itemsXMIPP[3]  in ListBlockKeys):
            _psi= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[3])
        for loopitem in range(blockLenght):
            angle.append( 
                          (
                           (float(_rot[loopitem] )) *math.pi/180.,\
                           (float(_tilt[loopitem])) *math.pi/180.,\
                           (float(_psi[loopitem] )) *math.pi/180.
                          )  
                        )
                        
                        

        _scale=[]
        
        for loopitem in range(blockLenght):
            _scale.append(1.)
        #fill scale
        if (self.needed_itemsXMIPP[8]  in ListBlockKeys):
            _scale= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[8])
        for loopitem in range(blockLenght):
            scale.append((_scale[loopitem],_scale[loopitem],_scale[loopitem]))
        #enable
        _enable=[]
        for loopitem in range(blockLenght):
            _enable.append(1)
        if (self.needed_itemsXMIPP[9]  in ListBlockKeys):
            _enable= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[9])
        #fom
        _fom=[]
        for loopitem in range(blockLenght):
            _fom.append(1.)
        if (self.needed_itemsXMIPP[10]  in ListBlockKeys):
            _fom= self.inMetadata[blockName].GetLoopItem(self.needed_itemsXMIPP[10])
        
        #convert this in matrices
        #store matrices here
        ListEulerMatrices = []

        for loopitem in range(blockLenght):
            ListEulerMatrices.append(compose_matrix(scale=scale[loopitem],
                    angles=angle[loopitem],
                    translate=shift[loopitem]))

        #Convert matrices to emx
        _emx_particle_url = _imageUrl
        _emx_particle_transformation_matrix_offset_x = _x
        _emx_particle_transformation_matrix_offset_y = _y
        _emx_particle_transformation_matrix_offset_z = _z
        _emx_particle_transformation_matrix_1_1=[]
        _emx_particle_transformation_matrix_1_2=[]
        _emx_particle_transformation_matrix_1_3=[]
        _emx_particle_transformation_matrix_2_1=[]
        _emx_particle_transformation_matrix_2_2=[]
        _emx_particle_transformation_matrix_2_3=[]
        _emx_particle_transformation_matrix_3_1=[]
        _emx_particle_transformation_matrix_3_2=[]
        _emx_particle_transformation_matrix_3_3=[]
        
        for loopitem in range(blockLenght):
            _emx_particle_transformation_matrix_1_1.append(ListEulerMatrices[loopitem][0][0])
            _emx_particle_transformation_matrix_1_2.append(ListEulerMatrices[loopitem][0][1])
            _emx_particle_transformation_matrix_1_3.append(ListEulerMatrices[loopitem][0][2])
            _emx_particle_transformation_matrix_2_1.append(ListEulerMatrices[loopitem][1][0])
            _emx_particle_transformation_matrix_2_2.append(ListEulerMatrices[loopitem][1][1])
            _emx_particle_transformation_matrix_2_3.append(ListEulerMatrices[loopitem][1][2])
            _emx_particle_transformation_matrix_3_1.append(ListEulerMatrices[loopitem][2][0])
            _emx_particle_transformation_matrix_3_2.append(ListEulerMatrices[loopitem][2][1])
            _emx_particle_transformation_matrix_3_3.append(ListEulerMatrices[loopitem][2][2])
        
        _emx_particle_enable  = _enable
        _emx_particle_FOM     = _fom
        self.cbOut.AddCifItem((self.needed_itemsEMX, [
                               [
                                _emx_particle_url,
                                _emx_particle_transformation_matrix_1_1,
                                _emx_particle_transformation_matrix_1_2,
                                _emx_particle_transformation_matrix_1_3,
                                _emx_particle_transformation_matrix_offset_x,
                                _emx_particle_transformation_matrix_2_1,
                                _emx_particle_transformation_matrix_2_2,
                                _emx_particle_transformation_matrix_2_3,
                                _emx_particle_transformation_matrix_offset_y,
                                _emx_particle_transformation_matrix_3_1,
                                _emx_particle_transformation_matrix_3_2,
                                _emx_particle_transformation_matrix_3_3,
                                _emx_particle_transformation_matrix_offset_z,
                                _emx_particle_enable,
                                _emx_particle_FOM
                                ]
                              ])) 
    def saveFile(self):
        outfile = open(self.outputFileName,"w")
        comment =  "#  Converted from " + self.xmippStartVersion + " to " + self.emxVersion
        comment += "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################" 
        outfile.write(self.outMetadata.WriteOut(_add=True,comment=comment,_email=self.contactMail))
    

    
if __name__ == '__main__':

    def command_line_options():
        """ add command line options here"""
        import optparse
        _usage = "usage: %prog [options] Example:   %prog -i input.xmd -o out.emx"
        parser = optparse.OptionParser(_usage)        
        parser.add_option("-i", "--inputFileName", dest="inputFileName",
                          default="/dev/stdin", type="string",
                          help="XMIPP Metadata file (particle picking)")                           
        parser.add_option("-o", "--outputFileName", dest="outputFileName",
                          default="/dev/stdout", type="string",
                          help="EMX   Metadata file (particle picking)")                           
       
        (options, args) = parser.parse_args()
        return(options.inputFileName,options.outputFileName)
    
    inputFileName,outputFileName  = command_line_options()
    instance = convertParticlePickingClass(inputFileName,outputFileName)
