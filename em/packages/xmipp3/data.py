# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This modules contains basic hierarchy
for specific Xmipp3 EM data objects
"""

from pyworkflow.em import *   
from xmipp import *
    
class XmippSetOfMicrographs(SetOfMicrographs):
    """Represents a set of Micrographs for Xmipp"""
    def __init__(self, **args):
        SetOfMicrographs.__init__(self, **args)
        
    def __iter__(self):
        """Iterate over the set of micrographs in the MetaData"""
        md = MetaData(self.getFileName())
        
        for objId in md:    
            m = Micrograph()
            m.setFileName(md.getValue(MDL_MICROGRAPH, objId))
            if self.hasCTF():
                m.ctfModel = XmippCTFModel(md.getValue(MDL_CTF_MODEL, objId)) 
            yield m

        
class XmippCTFModel(CTFModel):
    
    ctfParams = {
                 "samplingRate":MDL_CTF_SAMPLING_RATE,
                 "voltage":MDL_CTF_VOLTAGE,
                 "defocusU":MDL_CTF_DEFOCUSU,
                 "defocusV":MDL_CTF_DEFOCUSV,
                 "defocusAngle":MDL_CTF_DEFOCUS_ANGLE,
                 "sphericalAberration":MDL_CTF_CS,
                 "chromaticAberration":MDL_CTF_CA,
                 "energyLoss":MDL_CTF_ENERGY_LOSS,
                 "lensStability":MDL_CTF_LENS_STABILITY,
                 "convergenceCone":MDL_CTF_CONVERGENCE_CONE,
                 "longitudinalDisplacement":MDL_CTF_LONGITUDINAL_DISPLACEMENT,
                 "transversalDisplacement":MDL_CTF_TRANSVERSAL_DISPLACEMENT,
                 "q0":MDL_CTF_Q0,
                 "k":MDL_CTF_K,
                 "bgGaussianK":MDL_CTF_BG_GAUSSIAN_K,
                 "bgGaussianSigmaU":MDL_CTF_BG_GAUSSIAN_SIGMAU,
                 "bgGaussianSigmaV":MDL_CTF_BG_GAUSSIAN_SIGMAV,
                 "bgGaussianCU":MDL_CTF_BG_GAUSSIAN_CU,
                 "bgGaussianCV":MDL_CTF_BG_GAUSSIAN_CV,
                 "bgGaussianAngle":MDL_CTF_BG_GAUSSIAN_ANGLE,
                 "bgSqrtK":MDL_CTF_BG_SQRT_K,
                 "bgSqrtU":MDL_CTF_BG_SQRT_U,
                 "bgSqrtV":MDL_CTF_BG_SQRT_V,
                 "bgSqrtAngle":MDL_CTF_BG_SQRT_ANGLE,
                 "bgBaseline":MDL_CTF_BG_BASELINE,
                 "bgGaussian2K":MDL_CTF_BG_GAUSSIAN2_K,
                 "bgGaussian2SigmaU":MDL_CTF_BG_GAUSSIAN2_SIGMAU,
                 "bgGaussian2SigmaV":MDL_CTF_BG_GAUSSIAN2_SIGMAV,
                 "bgGaussian2CU":MDL_CTF_BG_GAUSSIAN2_CU,
                 "bgGaussian2CV":MDL_CTF_BG_GAUSSIAN2_CV,
                 "bgGaussian2Angle":MDL_CTF_BG_GAUSSIAN2_ANGLE,
#                 "X0":MDL_CTF_X0,
#                 "XF":MDL_CTF_XF,
#                 "Y0":MDL_CTF_Y0,
#                 "YF":MDL_CTF_YF,
                 "critFitting":MDL_CTF_CRIT_FITTINGSCORE,
                 "critCorr13":MDL_CTF_CRIT_FITTINGCORR13,
#                 "downsampleFactor":MDL_CTF_DOWNSAMPLE_PERFORMED,
                 "critPsdStdQ":MDL_CTF_CRIT_PSDVARIANCE,
                 "critPsdPCA1":MDL_CTF_CRIT_PSDPCA1VARIANCE,
                 "critPsdPCARuns":MDL_CTF_CRIT_PSDPCARUNSTEST
                 }
    
    # Implementar el constructor para crear las variables del modelo usando el params de arriba
    # y leyendo del metadata. No hace falta el __getattr__
    
    def __init__(self, filename):
        md = MetaData(filename)
        objId = md.firstObject()
        
        for key, val in  self.ctfParams.iteritems():
            mdVal = md.getValue(val, objId)
            if not hasattr(self, key):
                setattr(self, key, Float(mdVal))
            else:
                getattr(self, key).set(mdVal)
                
                
        
        
    
