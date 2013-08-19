from pyworkflow.em import *  
import xmipp
from data import *
import xmipp, xmipp3
from protocol_kerdensom import *



class XmippDefRotSpectra(XmippDefKerdensom):
    """Create the definition of parameters for
    the Rotational Spectra protocol"""
    def _addParams(self):
        self.addParam('howCenter', EnumParam, choices=['Middle of the image', 'Minimize first harmonic'], 
                      default=xmipp3.ROTSPECTRA_CENTER_MIDDLE, important=True, 
                      label='How to find the center of rotation', display=EnumParam.DISPLAY_COMBO, 
                      help='Select how to find the center of rotation.')  
        self.addParam('spectraInnerRadius', IntParam, default=15,
                      label='Inner radius for rotational harmonics (%)',
                      help='A percentage of the image radius', expertLevel=LEVEL_ADVANCED)
        self.addParam('spectraOuterRadius', IntParam, default=80,
                      label='Outer radius for rotational harmonics (%)',
                      help='A percentage of the image radius', expertLevel=LEVEL_ADVANCED)
        self.addParam('spectraLowHarmonic', IntParam, default=1,
                      label='Lowest harmonic to calculate',
                      expertLevel=LEVEL_ADVANCED)
        self.addParam('spectraHighHarmonic', IntParam, default=15,
                      label='Highest harmonic to calculate',
                      expertLevel=LEVEL_ADVANCED) 
        self.addSection(label='Kerdensom')
        XmippDefKerdensom._addParams(self)
          
#         self.addSection(label='Classification: classify kerdensom')
#         self.addParam('xDimension', IntParam, default=7,
#                       label='X-dimension of the self-organizing map')
#         self.addParam('yDimension', IntParam, default=7,
#                       label='Y-dimension of the self-organizing map')   
#         self.addParam('initialRegFactor', FloatParam, default=1000,
#                       label='Initial regularization factor',
#                       help='The kerdenSOM algorithm anneals from an initial high regularization factor \n' +
#                             'to a final lower one, in a user-defined number of steps. \n' +
#                             'If the output map is too smooth, lower the regularization factors \n' +
#                             'If the output map is not organized, higher the regularization factors \n' +
#                             'See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/KerDenSOM', expertLevel=LEVEL_ADVANCED)
#         self.addParam('finalRegFactor', FloatParam, default=200,
#                       label='Final regularization factor',
#                       expertLevel=LEVEL_ADVANCED)   
#         self.addParam('lowerRegFactorSteps', FloatParam, default=5,
#                       label='Steps to lower the regularization facor',
#                       expertLevel=LEVEL_ADVANCED) 
#         self.addParam('kerdenSOMParamters', IntParam,
#                       label='Additional kerdenSOM parameters',
#                       help='For a complete description. \n' +
#                       'See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/KerDenSOM', expertLevel=LEVEL_ADVANCED)
        
class XmippProtRotSpectra(XmippProtKerdensom):
    """Protocol to compute the rotational spectrum of the given particles"""
    
    _definition = XmippDefRotSpectra()
        
    def _prepareParams(self):
        XmippProtKerdensom._prepareParams(self)
        self._params['extraDir'] = self._getExtraPath()
        self._params['R1'] = self.spectraInnerRadius.get()
        self._params['R2'] = self.spectraOuterRadius.get()
        self._params['spectraLowHarmonic'] = self.spectraLowHarmonic.get()
        self._params['spectraHighHarmonic'] = self.spectraHighHarmonic.get()
        self._params['vectors'] = self._getExtraPath("rotSpectra.xmd")
    
    def _defineSteps(self):
        self._prepareParams()  
        imagesFn = self._params['imgsFn']
        centerFn = self._getExtraPath("center2d_center.xmd")
        # After any of this steps the file "center2d_center.xmd" should be produced
        if self.howCenter == xmipp3.ROTSPECTRA_CENTER_MIDDLE:
            self._insertMiddle(imagesFn, centerFn)
        else:
            self._insertFunctionStep('centerFirstHarmonic', imagesFn, centerFn)
        # Produce "rotSpectra.xmd" vectors
        self._insertFunctionStep('calculateSpectra', imagesFn, centerFn, self._params['vectors'])
        # Call kerdensom for classification
        self._insertKerdensom()
        # Store outputs in db
        self._insertFunctionStep('createOutput')
        
    def _insertMiddle(self, inputImages, outputCenter):
        R2 = self._params['R2']
        
        if R2 + 20 > 100:
            R3 = R2 + (100 - R2) / 2
            R4 = 100
        else:
            R3 = R2 + 10
            R4 = R2 + 20
        self._params['R3'] = R3
        self._params['R4'] = R4
        
        program = 'xmipp_image_find_center'
        args = '-i ' + inputImages
        args += ' --oroot %(extraDir)s/center2d --r1 %(R1)d --r2 %(R2)d --r3 %(R3)d --r4 %(R4)d'
        # Run the command with formatted parameters
        self._insertRunJobStep(program, args % self._params)
            
    def centerFirstHarmonic(self, inputImages, outputCenter):
        dims = xmipp.MetaDataInfo(str(inputImages))
        md = xmipp.MetaData()
        objId = md.addObject()
        md.setValue(xmipp.MDL_X, float(dims[0] / 2), objId)
        md.setValue(xmipp.MDL_Y, float(dims[1] / 2), objId)
        md.write(outputCenter)
        return [outputCenter] # this file should exists after the step
            
    def calculateSpectra(self, inputImages, inputCenter, outputSpectra):     
        md = xmipp.MetaData(inputCenter)
        objId = md.firstObject()
        self._params['xOffset'] = md.getValue(xmipp.MDL_X, objId)
        self._params['yOffset'] = md.getValue(xmipp.MDL_Y, objId)
        
        program = 'xmipp_image_rotational_spectra'
        args = "-i %s -o %s" % (inputImages, outputSpectra)
        args += ' --x0 %(xOffset)d --y0 %(yOffset)d --r1 %(R1)d --r2 %(R2)d' + \
                     ' --low %(spectraLowHarmonic)d --high %(spectraHighHarmonic)d'
        # Run the command with formatted parameters
        self.runJob(None, program, args % self._params)
        return [outputSpectra]     
    

#         
#     def _validate(self):
#         validateMsgs = []
#         if self.initialRegFactor < self.finalRegFactor:
#             validateMsgs.append("Regularization must decrease over iterations:")
#             validateMsgs.append("    Initial regularization must be larger than final")
#         return validateMsgs
        