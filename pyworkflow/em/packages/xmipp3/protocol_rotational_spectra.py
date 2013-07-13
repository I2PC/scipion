from pyworkflow.em import *  
import xmipp
from data import *
from xmipp3 import XmippProtocol
from pyworkflow.em.packages.xmipp3.protocol_kerdensom import *

CENTER_MIDDLE = 0
CENTER_FIRST_HARMONIC = 1


class XmippDefRotSpectra(XmippDefKerdensom):
    """Create the definition of parameters for
    the Rotational Spectra protocol"""
    def _addParams(self):
        self.addSection(label='Rotational spectra parameters')
        self.addParam('howCenter', EnumParam, choices=['Use the middle of the image', 'Minimize first harmonic'], 
                      default=CENTER_MIDDLE, important=True, label='How to find the center of rotation', display=EnumParam.DISPLAY_COMBO, 
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
        self.addSection(label='Classification: classify kerdensom')
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
    
    def _defineSteps(self):
        
#         # Convert input images if necessary
#         self.inputParts = self.inputParticles.get()        
#         partsFn = self._insertConvertStep('inputParts', XmippSetOfParticles,
#                                          self._getPath('input_particles.xmd'))
        
        XmippProtKerdensom._prepareDefinition(self)
        
        # Prepare arguments to call programs
        self._params = self._params + { 'extraDir': self._getExtraPath(),
#                                         'partsFn': partsFn,
                                        'howCenter' : self.howCenter.get(),
                                        'spectraInnerRadius' : self.spectraInnerRadius.get(),
                                        'spectraOuterRadius' : self.spectraOuterRadius.get(),
                                        'spectraLowHarmonic' : self.spectraLowHarmonic.get(),
                                        'spectraHighHarmonic' : self.spectraHighHarmonic.get()
                #                         'xDimension' : xDimension,
                #                         'yDimension' : yDimension,
                #                         'initialRegFactor' : initialRegFactor,
                #                         'finalRegFactor' : finalRegFactor,
                #                         'lowerRegFactorSteps' : lowerRegFactorSteps,
                #                         'kerdenSOMParamters' : kerdenSOMParamters                        
                                        }
        
        self._insertFindCenter()
        self._insertFunctionStep('calculateSpectra')
        XmippDefKerdensom._insertSteps()
        
    def _insertFindCenter(self):
        if self.howCenter == CENTER_MIDDLE:
            self._insertMiddle()
        else:
            self._insertFunctionStep('centerFirstHarmonic', )
        
    def _insertMiddle(self):
        centerFinder = self_params['howCenter']
        spectraOuterRadius = self_params['spectraOuterRadius']
        inputImages = self_params['inputImages']
        extraDir = self_params['extraDir']
        
        if spectraOuterRadius + 20 > 100:
            R3 = spectraOuterRadius + (100 - spectraOuterRadius) / 2
            R4 = 100
        else:
            R3 = spectraOuterRadius + 10
            R4 = spectraOuterRadius + 20
        program = 'xmipp_image_find_center'
        params = '-i %(inputImages)s --oroot %(extraDir)s/center2d --r1 %(spectraInnerRadius)i --r2 %(outerRadius)i' + \
                 ' --r3 ' + str(R3) + ' --r4 ' + str(R4)
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self._insertRunJobStep(program, arguments % self._params)
            
    def centerFirstHarmonic(self, inputImages, outputCenter):
            from xmipp import MetaDataInfo
            dims = MetaDataInfo(inputImages)
            md = xmipp.MetaData()
            id = md.addObject()
            md.setValue(xmipp.MDL_X, float(dims[0] / 2), id)
            md.setValue(xmipp.MDL_Y, float(dims[1] / 2), id)
            md.write(join(extraDir), "center2d_center.xmd")
            
    def calculateSpectra(self):
        extraDir = self_params['extraDir']        
        md = xmipp.MetaData(join(extraDir, "center2d_center.xmd"))
        id = md.firstObject()
        xOffset=md.getValue(MDL_X, id)
        yOffset=md.getValue(MDL_Y, id)
        
        program = 'xmipp_image_rotational_spectra'
        params = '-i %(inputImages)s -o %(extraDir)s/rotSpectra.xmd' + \
                 ' --x0 ' + str(xOffset) + \
                 ' --y0 ' + str(yOffset) + \
                 ' --r1 %(spectraInnerRadius)i --r2 %(outerRadius)i' + \
                 ' --low %(spectraLowHarmonic)i' + \
                 ' --high %(spectraHighHarmonic)i'
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self._insertRunJobStep(program, arguments % self._params)       
    
    
#     def kerdensom(self):
#         extraDir = self_params['extraDir']        
#         program = 'xmipp_classify_kerdensom'
#         params = '-i %(extraDir)s/rotSpectra.xmd' + \
#                  ' --oroot %(extraDir)s/kerdensom' + \
#                  ' --xdim %(xDimension)i' + \
#                  ' --ydim %(yDimension)i' + \
#                  ' --deterministic_annealing %(lowerRegFactorSteps)f %(initialRegFactor)f %(initialRegFactor)f' + \
#                  ' %(kerdenSOMParamters)s' 
#         # Run the command with formatted parameters
#         self._log.info('Launching: ' + program + ' ' + arguments % self._params)
#         self._insertRunJobStep(program, arguments % self._params)   
#         
#         deleteFiles(log, [os.path.join(extraDir,"rotSpectra.xmd"),os.path.join(extraDir,"rotSpectra.vec")], True)
#     
#         
#     def _validate(self):
#         validateMsgs = []
#         if self.initialRegFactor < self.finalRegFactor:
#             validateMsgs.append("Regularization must decrease over iterations:")
#             validateMsgs.append("    Initial regularization must be larger than final")
#         return validateMsgs
        