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
from pyworkflow.utils.path import makePath, moveFile, removeBaseExt, copyTree

from convert import *
from xmipp3 import XmippMdRow


class XmippCTFBase(ProtCTFMicrographs):
    """ This class cointains the common functionalities for all protocols to
estimate CTF on a set of micrographs using xmipp3 """

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
#        'ctffind_ctfparam': join('%(micDir)s', 'ctffind.ctfparam'),
#        'ctffind_spectrum': join('%(micDir)s', 'ctffind_spectrum.mrc')
        }
    
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
        self._params['micDir'] = self._getFilename('prefix', micDir=micDir)
        # CTF estimation with Xmipp                
        self.runJob(self._program, self._args % self._params)    
    
    def createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        defocusList = []
        
        for _, micDir, mic in self._iterMicrographs():
            ctfparam = self._getFilename('ctfparam', micDir=micDir)
            
            ctfModel = readCTFModel(ctfparam, mic)
            ctfSet.append(self._setPsdFiles(ctfModel, micDir))
            
            # save the values of defocus for each micrograph in a list
            defocusList.append(ctfModel.getDefocusU())
            defocusList.append(ctfModel.getDefocusV())
        
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(self.inputMics, ctfSet)
        self._defocusMaxMin(defocusList)
        
#         # Write as a Xmipp metadata
#         mdFn = self._getPath('micrographs_ctf.xmd')
#         writeSetOfMicrographs(micSet, mdFn, self.setupMicRow)
#                    
#         tmpFn = self._getPath('micrographs_backup.xmd')
#         writeSetOfMicrographs(micSet, tmpFn, self.setupMicRow)
#   
#         # Evaluate the PSD and add some criterias
#         auxMdFn = self._getTmpPath('micrographs.xmd')
#         self.runJob("xmipp_ctf_sort_psds","-i %s -o %s" % (mdFn, auxMdFn))
#         # Copy result to output metadata
#         moveFile(auxMdFn, mdFn)
        

    def _citations(self):
        return ['Vargas2013']
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareCommand(self):
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
    
    def _setPsdFiles(self, ctfModel, micDir):
        ctfModel._psdFile = String(self._getFilename('psd', micDir=micDir))
        ctfModel._xmipp_enhanced_psd = String(self._getFilename('enhanced_psd', micDir=micDir))
        ctfModel._xmipp_ctfmodel_quadrant = String(self._getFilename('ctfmodel_quadrant', micDir=micDir))
        ctfModel._xmipp_ctfmodel_halfplane = String(self._getFilename('ctfmodel_halfplane', micDir=micDir))
        return ctfModel


class XmippProtCTFMicrographs(XmippCTFBase):
    """Protocol to estimate CTF on a set of micrographs using xmipp3"""
    _label = 'ctf estimation'


class XmippProtRecalculateCTF(XmippCTFBase):
    """Protocol to recalculate CTF on a subSet of micrographs using xmipp3"""
    
    def __init__(self, **args):
        XmippCTFBase.__init__(self, **args)
    
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        
        form.addParam('inputCtf', PointerParam, important=True,
              label="input the SetOfCTF to recalculate", pointerClass='SetOfCTF')
        form.addParam('inputValues', TextParam)

        form.addParallelSection(threads=0, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform ctf estimation on a set of micrographs.
        """
        
        self._defineValues()
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each psd insert the steps to process it
        for line in self.values:
            # CTF Re-estimation with Xmipp
            copyId = self._insertFunctionStep('copyFiles', line, prerequisites=[])
            stepId = self._insertFunctionStep('_estimateCTF',line, prerequisites=[copyId]) # Make estimation steps independent between them
            deps.append(stepId)
        # Insert step to create output objects       
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def copyFiles(self, line):
        """Copy micrograph's directory tree"""
        objId = self._getObjId(line)
        ctfModel = self.setOfCtf.__getitem__(objId)
        mic = ctfModel.getMicrograph()
        
        prevDir = self._getPrevMicDir(mic)
        micDir = self._getMicrographDir(mic)
        # Create micrograph dir under extra directory
        print "creating path micDir=", micDir
        makePath(micDir)
        if not exists(micDir):
            raise Exception("No created dir: %s " % micDir)
        copyTree(prevDir, micDir)
    
    def _estimateCTF(self, line):
        """ Run the estimate CTF program """
        self._prepareCommand(line)
        # CTF estimation with Xmipp                
        self.runJob(self._program, self._args % self._params)    
    
    def createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        defocusList = []
        
        for ctfModel in self.setOfCtf:
            ctfNotUp = True
            
            for line in self.values:
                objId = self._getObjId(line)
                
                if objId == ctfModel.getObjId():
                    mic = ctfModel.getMicrograph()
                    mic.setObjId(ctfModel.getObjId())
                    micDir = self._getMicrographDir(mic)
                    ctfparam = self._getFilename('ctfparam', micDir=micDir)
                    ctfModel2 = readCTFModel(ctfparam, mic)
                    ctfModel2 = self._setPsdFiles(ctfModel2, micDir)
                    
                    # save the values of defocus for each micrograph in a list
                    defocusList.append(ctfModel2.getDefocusU())
                    defocusList.append(ctfModel2.getDefocusV())
                    ctfNotUp = False
                    break
            
            if ctfNotUp:
                ctfSet.append(ctfModel)
            else:
                ctfSet.append(ctfModel2)
        
        self._defineOutputs(outputCTF=ctfSet)
        self._defineSourceRelation(self.setOfCtf, ctfSet)
        self._defocusMaxMin(defocusList)
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputCTF'):
            summary.append(Message.TEXT_NO_CTF_READY)
        else:
            summary.append("CTF Re-estimation of *testing* micrographs.")
        return summary
    
    def _methods(self):
        methods = []
        
        if not hasattr(self, 'outputCTF'):
            methods.append(Message.TEXT_NO_CTF_READY)
        else:
            methods.append(self.methodsInfo.get())
            
        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _defineValues(self):
        """ This function get the acquisition info of the micrographs"""
        self.setOfCtf = self.inputCtf.get()
        self.values = []

        for l in self.inputValues.get().split(';'):
            list = []
            list = l.split(',')
            self.values.append(list)
            
        line = self.values[0]
        objId = self._getObjId(line)
        ctfModel = self.inputCtf.get().__getitem__(objId)
        mic = ctfModel.getMicrograph()
        
        acquisition = mic.getAcquisition()
        self._params = {'voltage': acquisition.getVoltage(),
                        'sphericalAberration': acquisition.getSphericalAberration(),
                        'magnification': acquisition.getMagnification(),
                        'ampContrast': acquisition.getAmplitudeContrast(),
                        'samplingRate': mic.getSamplingRate()
                       }
    
    def _prepareCommand(self, line):
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
    
    def _getPrevMicDir(self, mic):
        
        objFn = self.setOfCtf.getFileName()
        directory = dirname(objFn)
        return join(directory, "extra", removeBaseExt(mic.getFileName()))
    
    def _getObjId(self, list):
        return int(list[0])
    