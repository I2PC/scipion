import CifFile
import StarFile
from transformations import *
import numpy
smallNumber       = 0.00001

##########################################################################
#   General Class for Star data handling
##########################################################################

class EmxBase:
    """some constants"""
    xmippStartVersion = 'XMIPP_STAR_1'
    emxVersion        = 'EMX1.0'
    contactMail       = 'xmipp@cnb.csic.es'
    
#    def __init__(self, runWithoutArgs=False):
#        a=0

    
    def checkVersion(self):
        """ Check first line for EMX or XMIPP magic word. Abort if
        neither of these two words are available"""
        import os
        if not os.path.exists(self.inputFileName):
            print "File: ", self.inputFileName, "does not exists."
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
        print "Error: Metadata Files should contain the string ",\
               self.emxVersion, "or", self.xmippStartVersion, \
            " in the first line.  Exiting program.\n", \
            "First line: ", firstLine
        exit(1)
    
    def saveFileEMX2XMIPP(self):
        """Auxiliary funtion for saving metadata as xmipp file"""
        comment  =   "# XMIPP_STAR_1 *"
        comment += "\n##########################################################################"         
        comment +=  "\n#  Converted from " + self.emxVersion + " to " + self.xmippStartVersion
        comment +=  "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################" 
        outfile = open(self.outputFileName,"w")
        outfile.write(self.outMetadata.WriteOut(comment=comment,_email=self.contactMail))
        
        
    def saveFileXMIPP2EMX(self):
        """Auxiliary funtion for saving metadata as emx file"""
        outfile = open(self.outputFileName,"w")
        comment =  "#  Converted from " + self.xmippStartVersion + " to " + self.emxVersion
        comment += "\n#  Inputfile: " + self.inputFileName
        comment += "\n##########################################################################"
        outfile.write(self.outMetadata.WriteOut(_add=True,comment=comment,_email=self.contactMail))


##########################################################################
#   Class Related to Particle Picking Conversion
##########################################################################
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
    
    def run(self):
        emx2xmipp = self.checkVersion()
        #do not change the order: first checkversion then ciffile
        #otherwise stdin will be lost
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
        for blockName in self.inMetadata.keys():
            #create output block
            myblock = CifFile.CifBlock()
            self.outMetadata[blockName] = myblock
            self.cbOut = self.outMetadata[blockName]
            self.createDataHeaderXMIPP2EMX(blockName)
            self.convertLoopXMIPP2EMX(blockName)
            
    def createDataHeaderXMIPP2EMX(self,micrographName):
        """emx requires micrograph.url"""
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
        loopitems = self.inMetadata[micrographName].GetLoop((self.needed_itemsXMIPP[0])) #get item names and values
        self.cbOut.AddCifItem(([self.needed_itemsEMX],\
        [[loopitems[self.needed_itemsXMIPP[0]],loopitems[self.needed_itemsXMIPP[1]]]]))    


##########################################################################
#   Class Related to Alignment Conversion
##########################################################################
class ParticleAlignmentConverter(EmxBase):    

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
    
    needed_itemsEMX = (
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
            self.outMetadata[blockName] = myblock
            self.cbOut = self.outMetadata[blockName]
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

    def convertLoopEMX2XMIPP(self,blockName):
        #get item names and values
        loopitems = self.inMetadata[blockName].GetLoop(self.needed_itemsEMX[0])  
        #get block list
        ListBlockKeys = loopitems.keys()
        #get loop length
        blockLenght = len (self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[1]))
        #image is compulsory 
        _imageUrl   = self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[0])
        #store matrices here
        ListEulerMatrices = []
        for loopitem in range(blockLenght):
            ListEulerMatrices.append(numpy.identity(4))
        #fill shifts
        if (self.needed_itemsEMX[4]  in ListBlockKeys):
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[4])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][0][3]=aux[loopitem]
        
        if (self.needed_itemsEMX[8]  in ListBlockKeys):
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[8])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][1][3]=aux[loopitem]
        
        if (self.needed_itemsEMX[12]  in ListBlockKeys):
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[12])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][3]=aux[loopitem]
        # fill 3D rotations
        _2Dmatrix=True
        if (self.needed_itemsEMX[3]  in ListBlockKeys):
            _2Dmatrix=False
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[3])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][0][2]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[7])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][1][2]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[9])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][0]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[10])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][1]=aux[loopitem]
            aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[11])
            for loopitem in range(blockLenght):
                ListEulerMatrices[loopitem][2][2]=aux[loopitem]
        #2D rotations
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[1])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][0][0]=aux[loopitem]
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[2])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][0][1]=aux[loopitem]
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[5])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][1][0]=aux[loopitem]
        aux= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[6])
        for loopitem in range(blockLenght):
            ListEulerMatrices[loopitem][1][1]=aux[loopitem]

        #handle enable
        _enable = [1]*blockLenght
        if (self.needed_itemsEMX[13]  in ListBlockKeys):
            _enable= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[13])#enable
        #handle FOM
        _fom = [1.]*blockLenght
        if (self.needed_itemsEMX[14]  in ListBlockKeys):
            _fom= self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX[14])#enable
        
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

#        self.cbOut.AddCifItem(([self.needed_itemsXMIPP],[[_XList,_YList]])) 
            
        self.cbOut.AddCifItem(( [self.needed_itemsXMIPP], 
                                [[
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
                                ]]
                              )) 
       
    def convertLoopXMIPP2EMX(self,blockName):
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
        self.cbOut.AddCifItem(( [self.needed_itemsEMX], [
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

##########################################################################
#   Class Related to Particle Picking Conversion
##########################################################################
class ParticleClassConverter(EmxBase):    
    needed_itemsXMIPP =  (
          "_image",
          "_ref"#int
          )
    needed_itemsEMX = (
          "_emx_particle.url"
          )

    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn

        self.imageList = []
        self.refList   = []
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
        ref=1
        myblock = CifFile.CifBlock()
        self.outMetadata['class'] = myblock
        self.cbOut = self.outMetadata['class']#alias
        #self.cbOut.AddCifItem(self.needed_itemsEMX)

        for blockName in self.inMetadata.keys():
            self.convertLoopEMX2XMIPP(blockName,ref)
            ref += 1
        self.cbOut.AddCifItem(( [self.needed_itemsXMIPP],[[self.imageList,self.refList]] )) 
    
    def convertAllBlocksXMIPP2EMX(self):
        """loop over the blocks and write them"""
        #get all different references from xmipp
        #just one block
        _referenceList = self.inMetadata.first_block()
        blockNameXmipp = self.inMetadata.keys()[0]
        _referenceList = _referenceList.GetLoopItem(self.needed_itemsXMIPP[1])
        #get a list with unique classes
        #_referenceList = self.inMetadata[blockNameXmipp].GetLoopItem(self.needed_itemsXMIPP[0])
        _referenceList = list(set(_referenceList))
        for blockNameEmx in _referenceList:
            #create output block
            myblock = CifFile.CifBlock()
            self.outMetadata[blockNameEmx] = myblock
            self.cbOut = self.outMetadata[blockNameEmx]#for emx
            self.createDataHeaderXMIPP2EMX(blockNameEmx)
            self.convertLoopXMIPP2EMX(blockNameXmipp,blockNameEmx)#for xmipp

    def createDataHeaderXMIPP2EMX(self,class_id):
        """Data header is the xmipp data block name with the right label"""
        self.cbOut['_emx_class.id'] = class_id
            
    def convertLoopEMX2XMIPP(self,blockName,ref):
         _imageList = self.inMetadata[blockName].GetLoopItem(self.needed_itemsEMX)
         _ref=[ref]*len(_imageList)
         self.imageList += _imageList
         self.refList   += _ref
         #add to cif should be done at the end
         
    def convertLoopXMIPP2EMX(self,blockNameXmipp,blockNameEmx):
        #get item names and values
        loopimages  = self.inMetadata[blockNameXmipp].GetLoopItem(self.needed_itemsXMIPP[0])  
        loopclasses = self.inMetadata[blockNameXmipp].GetLoopItem(self.needed_itemsXMIPP[1])  
        print "loopimages",loopimages
        print "loopclasses",loopclasses
        print "blockNameEmx",blockNameEmx
        #get block list
        _imageList = []
        for imageItem in range (len(loopclasses)):
            if(loopclasses[imageItem]==blockNameEmx):
                _imageList.append(loopimages[imageItem])
        print "_imageList",_imageList
        
#        print [(self.needed_itemsEMX)], [
#                               [_imageList
#                                ]
#                              ]

        self.cbOut.AddCifItem((self.needed_itemsEMX, _imageList )) 

##########################################################################
#   Class Related to CTF conversion, second case with CTF per image is pending
##########################################################################
class CtfConverter(EmxBase):    
    needed_itemsXMIPP =  (
            "image",
            "CTF_Sampling_rate",
            "CTF_Defocus_U",
            "CTF_Defocus_V",
            "CTF_Defocus_angle",
            "CTF_Voltage",
            "CTF_Spherical_aberration",
            "CTF_Q0"
          )
    needed_itemsEMX = (
            "_emx_micrograph.url",                     
            "_emx_micrograph.magnification",    
            "_emx_micrograph.scanner_pixel_size",    
            "_emx_micrograph.defocusU",
            "_emx_micrograph.defocusV",
            "_emx_micrograph.astigmatism_angle",
            "_emx_micrograph.voltage",
            "_emx_micrograph.Cs",
            "_emx_micrograph.amplitude_contrast"
          )

    
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
    
    def run(self):
        emx2xmipp = self.checkVersion()
        #do not change the order: first checkversion then ciffile
        #otherwise stdin will be lost
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
            blockname = self.inMetadata[blockName]['_emx_micrograph.url'] 
            self.outMetadata[blockname] = myblock
            self.cbOut = self.outMetadata[blockname]#alias
            self.convertLoopEMX2XMIPP(blockName)
            
    def convertAllBlocksXMIPP2EMX(self):
        """loop over the blocks and write them"""
        for blockName in self.inMetadata.keys():
            #create output block
            myblock = CifFile.CifBlock()
            self.outMetadata[blockName] = myblock
            self.cbOut = self.outMetadata[blockName]
            self.createDataHeaderXMIPP2EMX(blockName)
            self.convertLoopXMIPP2EMX(blockName)
            
    def createDataHeaderXMIPP2EMX(self,micrographName):
        """emx requires micrograph.url"""
        self.cbOut['_emx_micrograph.url'] = micrographName

    def convertLoopEMX2XMIPP(self,micrographName):
        _emx_micrograph.url = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[0])
        _emx_micrograph.magnification = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[1])
        _emx_micrograph.scanner_pixel_size = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[2])
        _emx_micrograph.defocusU = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[3])
        _emx_micrograph.defocusV = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[4])
        _emx_micrograph.astigmatism_angle = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[5])
        _emx_micrograph.voltage = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[6])
        _emx_micrograph.Cs = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[7])
        _emx_micrograph.amplitude_contrast = self.inMetadata[micrographName].GetLoopItem(self.needed_itemsEMX[8])

        CTF_Sampling_rate =[]
        for _item in range (len(_emx_micrograph.magnification)): 
            CTF_Sampling_rate.append(str( 10000000*float(_emx_micrograph.scanner_pixel_size[_item])
                                        / float(_emx_micrograph.magnification[_item]))
                                    )
                        
        self.cbOut.AddCifItem(([self.needed_itemsXMIPP],[[
                                                          _emx_micrograph.url,
                                                          CTF_Sampling_rate,
                                                          _emx_micrograph.defocusU,
                                                         _emx_micrograph.defocusV,
                                                         _emx_micrograph.astigmatism_angle,
                                                         _emx_micrograph.voltage,
                                                         _emx_micrograph.Cs,
                                                         _emx_micrograph.amplitude_contrast
                                                          ]])) 

    def convertLoopXMIPP2EMX(self,micrographName):
        loopitems = self.inMetadata[micrographName].GetLoop((self.needed_itemsXMIPP[0])) #get item names and values
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
        
    inputFn,outputFn,convType=command_line_options()
    if convType == 'coordinates':
        ParticlePickingConverter(inputFn, outputFn).run()  
    elif convType == 'alignment':
        ParticleAlignmentConverter(inputFn, outputFn).run()  
    elif convType == 'class':
        ParticleClassConverter(inputFn, outputFn).run()  
    elif convType == 'ctf':
        CtfConverter(inputFn, outputFn).run()  
    else:
        print "ERROR: Wrong mode: ", convType
        exit(0)
    
    
    
        
