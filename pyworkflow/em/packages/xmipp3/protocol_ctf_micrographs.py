# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from convert import *
from pyworkflow.em.packages.xmipp3.utils import isMdEmpty
from pyworkflow.protocol.params import RelationParam



class XmippProtCTFMicrographs(ProtCTFMicrographs):
    """ Protocol to estimate CTF on a set of micrographs using Xmipp. """
    _label = 'ctf estimation'
    
    _criterion = ("ctfCritFirstZero<5 OR ctfCritMaxFreq>20 OR "
                  "ctfCritfirstZeroRatio<0.9 OR ctfCritfirstZeroRatio>1.1 OR "
                  "ctfCritFirstMinFirstZeroRatio>10 OR ctfCritCorr13<0 OR "
                  "ctfCritCtfMargin<0 OR ctfCritNonAstigmaticValidty<0.3 OR "
                  "ctfCritNonAstigmaticValidty>25")

    def __init__(self, **args):
        ProtCTFMicrographs.__init__(self, **args)
    
    def _createFilenameTemplates(self):
        prefix = join('%(micDir)s','xmipp_ctf')
        _templateDict = {
                         # This templates are relative to a micDir
                         'micrographs': 'micrographs.xmd',
                         'prefix': prefix,
                         'ctfparam': prefix +  '.ctfparam',
                         'ctfErrorParam': prefix +  '_error.ctfparam',
                         'psd': prefix + '.psd',
                         'enhanced_psd': prefix + '_enhanced_psd.xmp',
                         'ctfmodel_quadrant': prefix + '_ctfmodel_quadrant.xmp',
                         'ctfmodel_halfplane': prefix + '_ctfmodel_halfplane.xmp',
                         'ctf': prefix + '.xmd'
                         }
        self._updateFilenamesDict(_templateDict)
    
    def _defineProcessParams(self, form):
        form.addParam('doInitialCTF', BooleanParam, default=False, 
              label="Use defoci from a previous CTF estimation")
        form.addParam('ctfRelations', RelationParam, allowsNull=True,
                      condition='doInitialCTF',
                      relationName=RELATION_CTF,
                      attributeName='getInputMicrographs',
                      label='Previous CTF estimation',
                      help='Choose some CTF estimation related to input '
                           'micrographs, in case you want to use the defocus '
                           'values found previously')

        form.addParam('doCTFAutoDownsampling', BooleanParam, default=True,
                      label="Automatic CTF downsampling detection",
                      expertLevel=LEVEL_ADVANCED,
                      help='If this option is chosen, the algorithm '
                           'automatically tries by default the suggested '
                           'Downsample factor; and if it fails, +1; '
                           'and if it fails, -1.')
        form.addParam('doFastDefocus', BooleanParam, default=True,
                      label="Fast defocus estimate", expertLevel=LEVEL_ADVANCED,
                      help='Perform fast defocus estimate.')
        form.addParam('doAutomaticRejection', BooleanParam, default=False,
                      label="Automatic rejection", expertLevel=LEVEL_ADVANCED,
                      help='Automatically reject micrographs whose CTF looks '
                           'suspicious.')
    
    def getInputMicrographs(self):
        return self.inputMicrographs.get()

    #--------------------------- INSERT steps functions ------------------------
    def _insertFinalSteps(self, deps):
        stepId = self._insertFunctionStep('sortPSDStep', prerequisites=deps)
        return [stepId]
    
    def _getFirstJoinStepName(self):
        # This function will be used for streamming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (ej in Xmipp 'sortPSDStep')
        return 'sortPSDStep'
        
    #--------------------------- STEPS functions -------------------------------
    def _estimateCTF(self, micFn, micDir, micName):
        """ Run the estimate CTF program """
        localParams = self.__params.copy()

        if self.doInitialCTF:
            if self.ctfDict[micName] > 0:
                localParams['defocusU'] = self.ctfDict[micName]
                localParams['defocus_range'] = 0.01 * self.ctfDict[micName]
        else:
            ma = self._params['maxDefocus']
            mi = self._params['minDefocus']
            localParams['defocusU']= (ma + mi) / 2
            localParams['defocus_range'] = (ma - mi)/ 2
        
        # Create micrograph dir under extra directory
        makePath(micDir)
        if not exists(micDir):
            raise Exception("No created dir: %s " % micDir)

        finalName = micFn        
        ctfDownFactor = self.ctfDownFactor.get()
        downsampleList = [ctfDownFactor]

        if self.doCTFAutoDownsampling:
            downsampleList.append(ctfDownFactor+1)
            if ctfDownFactor >= 2:
                downsampleList.append(ctfDownFactor-1)
            else:
                if ctfDownFactor > 1:
                    downsampleList.append((ctfDownFactor+1)/2)
    
        deleteTmp = ""

        for downFactor in downsampleList:
            # Downsample if necessary
            if downFactor != 1:
                # Replace extension by 'mrc' cause there are some formats that
                # cannot be written (such as dm3)
                finalName = self._getTmpPath(replaceBaseExt(micFn, 'mrc'))
                self.runJob("xmipp_transform_downsample",
                            "-i %s -o %s --step %f --method fourier"
                            % (micFn, finalName, downFactor))
                deleteTmp=finalName
      
            # Update _params dictionary with mic and micDir
            localParams['micFn'] = finalName
            localParams['micDir'] = self._getFileName('prefix', micDir=micDir)
            localParams['samplingRate'] = self.inputMics.getSamplingRate() * downFactor
            
            # CTF estimation with Xmipp
            try:
                self.runJob(self._program,
                            self._args % localParams +
                            " --downSamplingPerformed %f" % downFactor)
            except Exception:
                break
            
            # Check the quality of the estimation and reject it necessary
            if self.evaluateSingleMicrograph(micFn,micDir):
                break
            
        if deleteTmp != "":
            cleanPath(deleteTmp)
    
    def _restimateCTF(self, micId):
        """ Run the estimate CTF program """
        ctfModel = self.recalculateSet[micId]
        self._prepareRecalCommand(ctfModel)
        # CTF estimation with Xmipp                
        self.runJob(self._program, self._args % self._params)
        
        mic=ctfModel.getMicrograph()
        micDir = self._getMicrographDir(mic)
        self.evaluateSingleMicrograph(mic.getFileName(),micDir)
    
    def sortPSDStep(self):
        # Gather all metadatas of all micrographs
        md=xmipp.MetaData()
        for _, micDir, mic in self._iterMicrographs():
            fnCTF=self._getFileName('prefix', micDir=micDir)+".xmd"
            enable=1
            if not os.path.exists(fnCTF):
                fnCTF = self._createErrorCtfParam(micDir)
                enable=-1
            mdCTF=xmipp.MetaData(fnCTF)
            mdCTF.setValue(xmipp.MDL_ENABLED,enable,mdCTF.firstObject())
            mdCTF.setValue(xmipp.MDL_MICROGRAPH_ID,long(mic.getObjId()),mdCTF.firstObject())
            md.unionAll(mdCTF)
        fnAllMicrographs=self._getPath("micrographs.xmd")
        md.write(fnAllMicrographs)
        
        # Now evaluate them
        self.runJob("xmipp_ctf_sort_psds","-i %s"%fnAllMicrographs)
    
        # And reject
        fnRejected=self._getPath("micrographs_rejected.xmd")
        self.runJob("xmipp_metadata_utilities",
                    '-i %s --query select "%s" -o %s'
                    % (fnAllMicrographs, self._criterion, fnRejected))
        mdRejected = xmipp.MetaData(fnRejected)
        self.someMicrographsRejected = mdRejected.size()>0

        if self.someMicrographsRejected:
            rejectedMicrographs = mdRejected.getColumnValues(xmipp.MDL_MICROGRAPH_ID)
            for mdid in md:
                micId = md.getValue(xmipp.MDL_MICROGRAPH_ID,mdid)
                if micId in rejectedMicrographs:
                    md.setValue(xmipp.MDL_ENABLED,-1,mdid)
                else:
                    md.setValue(xmipp.MDL_ENABLED,1,mdid)
        md.write(self._getPath("ctfs_selection.xmd"))
        cleanPath(fnRejected)
        
    def _createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        inputMics = self.inputMicrographs.get()
        ctfSet.setMicrographs(inputMics)
        defocusList = []
        if self.doAutomaticRejection:
            fn = self._getPath("micrographs.xmd")
        else:
            fn = self._getPath("ctfs_selection.xmd")
        mdFn = md.MetaData(fn)
        mdAux = md.MetaData()
        
        for _, micDir, mic in self._iterMicrographs():
            if exists(self._getFileName('ctfparam', micDir=micDir)):
                mdCTF = md.MetaData(self._getFileName('ctfparam', micDir=micDir))
                mdQuality = md.MetaData(self._getFileName('ctf', micDir=micDir))
                mdCTF.setValue(xmipp.MDL_CTF_CRIT_NONASTIGMATICVALIDITY,
                               mdQuality.getValue(xmipp.MDL_CTF_CRIT_NONASTIGMATICVALIDITY,mdQuality.firstObject()),
                               mdCTF.firstObject())
                mdCTF.setValue(xmipp.MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO,
                               mdQuality.getValue(xmipp.MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO,mdQuality.firstObject()),
                               mdCTF.firstObject())
            else:
                fnError=self._getFileName('ctfErrorParam', micDir=micDir)
                if not exists(fnError):
                    self._createErrorCtfParam(micDir)
                mdCTF = md.MetaData(fnError)
            mdAux.importObjects( mdFn, md.MDValueEQ(md.MDL_MICROGRAPH_ID, long(mic.getObjId())))
            mdAux.merge(mdCTF)
            
            mic.setSamplingRate(mdCTF.getValue(md.MDL_CTF_SAMPLING_RATE, 1))
            
            ctfModel = mdToCTFModel(mdAux, mic)
            self._setPsdFiles(ctfModel, micDir)
            ctfSet.append(ctfModel)
            
            # save the values of defocus for each micrograph in a list
            defocusList.append(ctfModel.getDefocusU())
            defocusList.append(ctfModel.getDefocusV())
                
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(inputMics, ctfSet)
        self._computeDefocusRange(ctfSet)

    def _checkNewCTFs(self, micSet):
        # In Xmipp, we need to sort the resulting CTF at the end,
        # so we can not update the output in streaming.
        return None, None
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        validateMsgs = []
        # downsampling factor must be greater than 1
        if self.ctfDownFactor.get() < 1:
            validateMsgs.append('Downsampling factor must be >=1.')
        if self.doInitialCTF:
            if not self.ctfRelations.hasValue():
                validateMsgs.append('If you want to use a previous estimation '
                                    'of the CTF, the corresponding set of CTFs '
                                    'is needed')
        return validateMsgs
    
    def _summary(self):
        summary = ProtCTFMicrographs._summary(self)
        if self.methodsVar.hasValue():
            summary.append(self.methodsVar.get())
        return summary
    
    def _methods(self):
        strMsg = "We calculated the CTF of micrographs %s using Xmipp [Sorzano2007a]" % self.getObjectTag('inputMicrographs')
        if self.doFastDefocus and not self.doInitialCTF:
            strMsg += " with a fast defocus estimate [Vargas2013a]"
        strMsg += "."

        if self.methodsVar.hasValue():
            strMsg += " " + self.methodsVar.get()
        
        if self.hasAttribute('outputCTF'):
            strMsg += '\nOutput set is %s.'%self.getObjectTag('outputCTF')
        
        return [strMsg]
    
    def _citations(self):
        papers = ['Sorzano2007a']
        if self.doFastDefocus and not self.doInitialCTF:
            papers.append('Vargas2013a')
        return papers
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareArgs(self,params):
        self._args = ("--micrograph %(micFn)s --oroot %(micDir)s --sampling_rate "
                      "%(samplingRate)s --defocusU %(defocusU)f --defocus_range "
                      "%(defocus_range)f --overlap 0.7 ")
        
        for par, val in params.iteritems():
            self._args += " --%s %s" % (par, str(val))
            
        if self.doFastDefocus and not self.doInitialCTF:
            self._args += " --fastDefocus"

    def _prepareCommand(self):
        self._createFilenameTemplates()
        self._program = 'xmipp_ctf_estimate_from_micrograph'       
        
        # Mapping between base protocol parameters and the package specific
        # command options
        self.__params = {'kV': self._params['voltage'],
                'Cs': self._params['sphericalAberration'],
                'ctfmodelSize': self._params['windowSize'],
                'Q0': self._params['ampContrast'],
                'min_freq': self._params['lowRes'],
                'max_freq': self._params['highRes'],
                'pieceDim': self._params['windowSize']
                }
        self._prepareArgs(self.__params)

        if self.ctfRelations.hasValue():
            self.ctfDict = {}
            for ctf in self.ctfRelations.get():
                ctfName = ctf.getMicrograph().getMicName()
                self.ctfDict[ctfName] = ctf.getDefocusU()
                    
    def _prepareRecalCommand(self, ctfModel):
        if self.recalculate:
            self._defineRecalValues(ctfModel)
            self._createFilenameTemplates()
            self._program = 'xmipp_ctf_estimate_from_psd'       
            self._args = "--psd %(psdFn)s "
            line = ctfModel.getObjComment().split()
    
            # get the size and the image of psd
    
            imgPsd = ctfModel.getPsdFile()
            psdFile = basename(imgPsd)
            imgh = ImageHandler()
            size, _, _, _ = imgh.getDimensions(imgPsd)
            
            mic = ctfModel.getMicrograph()
            micDir = self._getMicrographDir(mic)
            fnCTFparam = self._getFileName('ctfparam', micDir=micDir)
            mdCTFParam = xmipp.MetaData(fnCTFparam)
            downFactor = mdCTFParam.getValue(xmipp.MDL_CTF_DOWNSAMPLE_PERFORMED,
                                             mdCTFParam.firstObject())
            # cleanPath(fnCTFparam)
            
            params2 = {'psdFn': join(micDir, psdFile),
                       'defocusU': float(line[0]),
                       'defocusV': float(line[1]),
                       'angle': line[2],
                      }
            self._params = dict(self._params.items() + params2.items())
            
            # Mapping between base protocol parameters and the package specific command options
            self.__params = {'sampling_rate': self._params['samplingRate']*downFactor,
                             'downSamplingPerformed': downFactor,
                             'kV': self._params['voltage'],
                             'Cs': self._params['sphericalAberration'],
                             'min_freq': line[3],
                             'max_freq': line[4],
                             'defocusU': self._params['defocusU'],
                             'defocusV': self._params['defocusV'],
                             'azimuthal_angle': self._params['angle'],
                             'Q0': self._params['ampContrast'],
                             'defocus_range': 5000,
                             'ctfmodelSize': size
                            }
            
            for par, val in self.__params.iteritems():
                self._args += " --%s %s" % (par, str(val))
    
    def _setPsdFiles(self, ctfModel, micDir):
        ctfModel._psdFile = String(self._getFileName('psd', micDir=micDir))
        ctfModel._xmipp_enhanced_psd = String(self._getFileName('enhanced_psd', micDir=micDir))
        ctfModel._xmipp_ctfmodel_quadrant = String(self._getFileName('ctfmodel_quadrant', micDir=micDir))
        ctfModel._xmipp_ctfmodel_halfplane = String(self._getFileName('ctfmodel_halfplane', micDir=micDir))
    
    def evaluateSingleMicrograph(self,micFn,micDir):
        mdCTFparam=xmipp.MetaData(self._getFileName('ctfparam', micDir=micDir))
        ctfDownFactor=mdCTFparam.getValue(xmipp.MDL_CTF_DOWNSAMPLE_PERFORMED,mdCTFparam.firstObject())
        
        mdEval = md.MetaData()
        objId = mdEval.addObject()

# Check what happen with movies        
#             if md.MDL_TYPE == md.MDL_MICROGRAPH_MOVIE:
#                 mdEval.setValue(md.MDL_MICROGRAPH_MOVIE, micFn, objId)
#                 mdEval.setValue(md.MDL_MICROGRAPH,('%05d@%s'%(1, micFn)), objId)
#             else:
        mdEval.setValue(md.MDL_MICROGRAPH, micFn, objId)
        
        mdEval.setValue(md.MDL_PSD, str(self._getFileName('psd', micDir=micDir)), objId)
        mdEval.setValue(md.MDL_PSD_ENHANCED, str(self._getFileName('enhanced_psd', micDir=micDir)), objId)
        mdEval.setValue(md.MDL_CTF_MODEL, str(self._getFileName('ctfparam', micDir=micDir)), objId)
        mdEval.setValue(md.MDL_IMAGE1, str(self._getFileName('ctfmodel_quadrant', micDir=micDir)), objId)
        mdEval.setValue(md.MDL_IMAGE2, str(self._getFileName('ctfmodel_halfplane', micDir=micDir)), objId)
        mdEval.setValue(md.MDL_CTF_DOWNSAMPLE_PERFORMED,float(ctfDownFactor), objId)
        fnEval = self._getFileName('ctf', micDir=micDir)
        mdEval.write(fnEval)
        
        # Evaluate if estimated ctf is good enough
        # self.runJob("xmipp_ctf_sort_psds","-i %s --downsampling %f"%(fnEval,ctfDownFactor))   
        try:
            self.runJob("xmipp_ctf_sort_psds","-i %s"%(fnEval))
        except Exception:
            pass
#             print >> sys.stderr, "xmipp_ctf_sort_psds has been Failed!"
#             raise ex

        # Check if it is a good micrograph
        fnRejected = self._getTmpPath(basename(micFn +"_rejected.xmd"))
        self.runJob("xmipp_metadata_utilities",'-i %s --query select "%s" -o %s' % (fnEval, self._criterion, fnRejected))
        
        retval=True
        
        if not isMdEmpty(fnRejected):
            mdCTF = md.MetaData(fnEval)
            mdCTF.setValue(md.MDL_ENABLED, -1, mdCTF.firstObject())
            mdCTF.write(fnEval)
            retval=False
        cleanPath(fnRejected)
        return retval
    
    def _createCtfModel(self, mic):
        micDir = self._getMicrographDir(mic)
        ctfparam = self._getFileName('ctfparam', micDir=micDir)
        ctfModel2 = readCTFModel(ctfparam, mic)
        self._setPsdFiles(ctfModel2, micDir)
        return ctfModel2
    
    def _createErrorCtfParam(self, micDir):
                ctfparam = self._getFileName('ctfErrorParam', micDir=micDir)
                f = open(ctfparam, 'w+')
                lines = """# XMIPP_STAR_1 * 
# 
data_fullMicrograph
 _ctfSamplingRate -999
 _ctfVoltage -999
 _ctfDefocusU -999
 _ctfDefocusV -999
 _ctfDefocusAngle -999
 _ctfSphericalAberration -999
 _ctfChromaticAberration -999
 _ctfEnergyLoss -999
 _ctfLensStability -999
 _ctfConvergenceCone -999
 _ctfLongitudinalDisplacement -999
 _ctfTransversalDisplacement -999
 _ctfQ0 -999
 _ctfK -999
 _ctfBgGaussianK -999
 _ctfBgGaussianSigmaU -999
 _ctfBgGaussianSigmaV -999
 _ctfBgGaussianCU -999
 _ctfBgGaussianCV -999
 _ctfBgGaussianAngle -999
 _ctfBgSqrtK -999
 _ctfBgSqrtU -999
 _ctfBgSqrtV -999
 _ctfBgSqrtAngle -999
 _ctfBgBaseline -999
 _ctfBgGaussian2K -999
 _ctfBgGaussian2SigmaU -999
 _ctfBgGaussian2SigmaV -999
 _ctfBgGaussian2CU -999
 _ctfBgGaussian2CV -999
 _ctfBgGaussian2Angle -999
 _ctfX0 -999
 _ctfXF -999
 _ctfY0 -999
 _ctfYF -999
 _ctfCritFitting -999
 _ctfCritCorr13 -999
 _CtfDownsampleFactor -999
 _ctfCritPsdStdQ -999
 _ctfCritPsdPCA1 -999
 _ctfCritPsdPCARuns -999
"""
                f.write(lines)
                f.close()
                return ctfparam
