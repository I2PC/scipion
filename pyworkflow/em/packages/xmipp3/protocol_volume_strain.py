# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
This sub-package contains wrapper to calculate the local strain of a volume
"""
import os

from pyworkflow.object import Float, String
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam
from pyworkflow.em.protocol.protocol import EMProtocol
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.utils import cleanPath
from xmipp3 import getMatlabEnviron
# import xmipp

        
class XmippProtVolumeStrain(ProtAnalysis3D):
    """Compare two states of a volume to analyze the local strains and rotations"""
    _label = 'calculate strain'
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume0', PointerParam, label="Initial state", important=True,
                      pointerClass='Volume',
                      help='Initial state of the structure, it will be deformed to fit into the final state')
        form.addParam('inputVolumeF', PointerParam, label="Final state", important=True,
                      pointerClass='Volume',
                      help='Initial state of the structure, it will be deformed to fit into the final state')
        form.addParam('inputMask', PointerParam, label="Mask for the final state", important=True,
                      pointerClass='VolumeMask',
                      help='Binary mask that defines where the strains and rotations will be calculated')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        fnVol0 = self.inputVolume0.get().getFileName()
        fnVolF = self.inputVolumeF.get().getFileName()
        fnMask = self.inputMask.get().getFileName()
        self._insertFunctionStep("calculateStrain",fnVol0,fnVolF,fnMask)
        self._insertFunctionStep("prepareOutput")
        self._insertFunctionStep("createChimeraScript")
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def calculateStrain(self, fnVol0, fnVolF, fnMask):
        fnRoot=self._getExtraPath('result')
        mirtDir = os.path.join(os.environ['XMIPP_HOME'], 'external', 'mirt')
        # -wait -nodesktop
        args='''-r "diary('%s'); xmipp_calculate_strain('%s','%s','%s','%s'); exit"'''%(fnRoot+"_matlab.log",fnVolF,fnVol0,fnMask,fnRoot)
        self.runJob("matlab", args, env=getMatlabEnviron(mirtDir))
    
    def prepareOutput(self):
        volDim = self.inputVolume0.get().getDim()[0]
        fnRoot=self._getExtraPath('result')
        self.runJob("xmipp_image_convert", "-i %s_initial.raw#%d,%d,%d,0,float -o %s_initial.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_image_convert", "-i %s_final.raw#%d,%d,%d,0,float -o %s_final.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_image_convert", "-i %s_initialDeformedToFinal.raw#%d,%d,%d,0,float -o %s_initialDeformedToFinal.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_image_convert", "-i %s_strain.raw#%d,%d,%d,0,float -o %s_strain.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("xmipp_image_convert", "-i %s_localrot.raw#%d,%d,%d,0,float -o %s_localrot.vol"%
                    (fnRoot,volDim,volDim,volDim,fnRoot))
        self.runJob("rm","-f "+self._getExtraPath('result_*.raw'))
    
    def createChimeraScript(self):
        fnRoot = "extra/result"
        scriptFile = self._getPath('result') + '_strain_chimera.cmd'
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % (fnRoot+"_final.vol"))
        fhCmd.write("open %s\n" % (fnRoot+"_strain.vol"))
        fhCmd.write("vol #1 hide\n")
        fhCmd.write("scolor #0 volume #1 cmap rainbow reverseColors True\n")
        fhCmd.close()

        scriptFile = self._getPath('result') + '_localrot_chimera.cmd'
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % (fnRoot+"_final.vol"))
        fhCmd.write("open %s\n" % (fnRoot+"_localrot.vol"))
        fhCmd.write("vol #1 hide\n")
        fhCmd.write("scolor #0 volume #1 cmap rainbow reverseColors True\n")
        fhCmd.close()

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
#        vol = self.inputVolume.get()
#        xDim = self._getDimensions()
#        volDim = vol.getDim()[0]
#        
#        if volDim != xDim:
#            errors.append("Make sure that the volume and the images have the same size")
        return errors    
    
    def _summary(self):
        summary = []
#        summary.append("Images evaluated: %i" % self.inputSet.get().getSize())
#        summary.append("Volume: %s" % self.inputVolume.getNameId())
#        summary.append("symmetry: %s" % self.symmetryGroup.get())
        return summary
    
    def _methods(self):
        methods = []
#        if hasattr(self, 'outputClasses') or hasattr(self, 'outputAverages'):
#            methods.append("We evaluated %i images regarding to volume %s"
#                           " using %s symmetry" %(self.inputSet.get().getSize(),\
#                                                  self.inputVolume.getNameId(), self.symmetryGroup.get()) )
        return methods
    