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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains the XmippCtfMicrographs protocol
"""

from pyworkflow.em import *  
from pyworkflow.utils.path import makePath, moveFile, removeBaseExt
from convert import *


class XmippCTFBase():
    """ This class cointains the common functionalities for all protocols to
estimate CTF on a set of micrographs using xmipp3 """

    def _createFilenameTemplates(self):
        __prefix = join('%(micDir)s','xmipp_ctf')
        _templateDict = {
                         # This templates are relative to a micDir
                         'micrographs': 'micrographs.xmd',
                         'prefix': __prefix,
                         'ctfparam': __prefix +  '.ctfparam',
                         'psd': __prefix + '.psd',
                         'enhanced_psd': __prefix + '_enhanced_psd.xmp',
                         'ctfmodel_quadrant': __prefix + '_ctfmodel_quadrant.xmp',
                         'ctfmodel_halfplane': __prefix + '_ctfmodel_halfplane.xmp'
                         }
        self._updateFilenamesDict(_templateDict)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _setPsdFiles(self, ctfModel, micDir):
        ctfModel._psdFile = String(self._getFileName('psd', micDir=micDir))
        ctfModel._xmipp_enhanced_psd = String(self._getFileName('enhanced_psd', micDir=micDir))
        ctfModel._xmipp_ctfmodel_quadrant = String(self._getFileName('ctfmodel_quadrant', micDir=micDir))
        ctfModel._xmipp_ctfmodel_halfplane = String(self._getFileName('ctfmodel_halfplane', micDir=micDir))
    
    def _citations(self):
        return ['Vargas2013']


class XmippProtCTFMicrographs(ProtCTFMicrographs, XmippCTFBase):
    """Protocol to estimate CTF on a set of micrographs using xmipp3"""
    _label = 'ctf estimation'
    
    def __init__(self, **args):
        ProtCTFMicrographs.__init__(self, **args)

    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Run the estimate CTF program """        
        # Create micrograph dir under extra directory
        print "creating path micDir=", micDir
        makePath(micDir)
        if not exists(micDir):
            raise Exception("No created dir: %s " % micDir)
        # Update _params dictionary with mic and micDir
        self._params['micFn'] = micFn
        self._params['micDir'] = self._getFileName('prefix', micDir=micDir)
        # CTF estimation with Xmipp  
        try:              
            self.runJob(self._program, self._args % self._params)
        except Exception:
            self._log.info("FAILED ESTIMATION FOR: " + micFn)    
    
    def createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(self.inputMicrographs.get())
        defocusList = []
        
        for _, micDir, mic in self._iterMicrographs():
            ctfparam = self._getFileName('ctfparam', micDir=micDir)
            
            if not os.path.exists(ctfparam):
                ctfparam = self._createErrorCtfParam(micDir)
                
            ctfModel = readCTFModel(ctfparam, mic)
            self._setPsdFiles(ctfModel, micDir)
            ctfSet.append(ctfModel)
            
            # save the values of defocus for each micrograph in a list
            defocusList.append(ctfModel.getDefocusU())
            defocusList.append(ctfModel.getDefocusV())
                
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(self.inputMics, ctfSet)
        self._defocusMaxMin(defocusList)
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareCommand(self):
        self._createFilenameTemplates()
        self._program = 'xmipp_ctf_estimate_from_micrograph'       
        self._args = "--micrograph %(micFn)s --oroot %(micDir)s --fastDefocus"
        
        # Mapping between base protocol parameters and the package specific command options
        self.__params = {'kV': self._params['voltage'],
                'Cs': self._params['sphericalAberration'],
                'sampling_rate': self._params['samplingRate'], 
                'ctfmodelSize': self._params['windowSize'],
                'Q0': self._params['ampContrast'],
                'min_freq': self._params['lowRes'],
                'max_freq': self._params['highRes'],
                'pieceDim': self._params['windowSize'],
                'defocusU': (self._params['maxDefocus']+self._params['minDefocus'])/2,
                'defocus_range': (self._params['maxDefocus']-self._params['minDefocus'])/2
                }
        
        for par, val in self.__params.iteritems():
            self._args += " --%s %s" % (par, str(val))
    
    def _createErrorCtfParam(self, micDir):
                ctfparam = join(micDir, 'xmipp_error_ctf.ctfparam')
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


class XmippProtRecalculateCTF(ProtRecalculateCTF, XmippCTFBase):
    """Protocol to recalculate CTF on a subSet of micrographs using xmipp3"""
    _label = 'ctf re-estimation'
    
    def __init__(self, **args):
        ProtRecalculateCTF.__init__(self, **args)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, line):
        """ Run the estimate CTF program """
        self._prepareCommand(line)
        # CTF estimation with Xmipp                
        self.runJob(self._program, self._args % self._params)    
    
    def createOutputStep(self):
        self._createFilenameTemplates()
        setOfMics = self.inputCtf.get().getMicrographs()
        outputMics = self._createSetOfMicrographs("_subset")
        outputMics.copyInfo(setOfMics)
        ctfSet = self._createSetOfCTF("_recalculated")
        defocusList = []
        for ctfModel in self.setOfCtf:
            mic = ctfModel.getMicrograph()
            outputMics.append(mic)
            for line in self.values:
                objId = self._getObjId(line)
                
                if objId == ctfModel.getObjId():
                    mic.setObjId(ctfModel.getObjId())
                    micDir = self._getMicrographDir(mic)
                    ctfparam = self._getFileName('ctfparam', micDir=micDir)
                    ctfModel2 = readCTFModel(ctfparam, mic)
                    self._setPsdFiles(ctfModel2, micDir)
                    ctfModel.copy(ctfModel2)
                    
                    # save the values of defocus for each micrograph in a list
                    defocusList.append(ctfModel2.getDefocusU())
                    defocusList.append(ctfModel2.getDefocusV())
                    break
            ctfSet.append(ctfModel)
            
        self.setOfCtf.close()
        ctfSet.setMicrographs(outputMics)
        
        self._defineOutputs(outputMicrographs=outputMics)
        self._defineOutputs(outputCTF=ctfSet)
        
        self._defineTransformRelation(setOfMics, outputMics)
        self._defineSourceRelation(self.inputCtf.get(), ctfSet)
        self._defocusMaxMin(defocusList)
        self._ctfCounter(defocusList)
        
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareCommand(self, line):
        
        self._defineValues(line)
        self._createFilenameTemplates()
        self._program = 'xmipp_ctf_estimate_from_psd'       
        self._args = "--psd %(psdFn)s "
        
        # get the size and the image of psd
        objId = self._getObjId(line)
        ctfModel = self.setOfCtf.__getitem__(objId)
        imgPsd = ctfModel.getPsdFile()
        psdFile = basename(imgPsd)
        imgh = ImageHandler()
        size, _, _, _ = imgh.getDimensions(imgPsd)
        
        mic = ctfModel.getMicrograph()
        micDir = self._getMicrographDir(mic)
        cleanPath(self._getFileName('ctfparam', micDir=micDir))
        
        params2 = {'psdFn': join(micDir, psdFile),
                   'defocusU': float(line[1]),
                   'defocusV': float(line[2]),
                   'angle': line[3],
                  }
        self._params = dict(self._params.items() + params2.items())
        
        # Mapping between base protocol parameters and the package specific command options
        self.__params = {'sampling_rate': self._params['samplingRate'],
                         'kV': self._params['voltage'],
                         'Cs': self._params['sphericalAberration'],
                         'min_freq': line[4],
                         'max_freq': line[5],
                         'defocusU': self._params['defocusU'],
                         'defocusV': self._params['defocusV'],
                         'azimuthal_angle': self._params['angle'],
                         'Q0': self._params['ampContrast'],
                         'defocus_range': 5000,
                         'ctfmodelSize': size
                        }
        
        for par, val in self.__params.iteritems():
            self._args += " --%s %s" % (par, str(val))
    