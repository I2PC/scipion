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
This module contains the protocol for 3d refinement with Relion.
"""

from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtRefine3D

from protocol_base import *



class ProtRelionRefine3D(ProtRefine3D, ProtRelionBase):
    """Protocol to refine a 3D map using Relion. Relion employs an empirical
Bayesian approach to refinement of (multiple) 3D reconstructions
or 2D class averages in electron cryo-microscopy (cryo-EM). Many
parameters of a statistical model are learned from the data,which
leads to objective and high-quality results.
    """    
    _label = '3D refine'
    IS_CLASSIFY = False
    CHANGE_LABELS = [xmipp.MDL_AVG_CHANGES_ORIENTATIONS, 
                     xmipp.MDL_AVG_CHANGES_OFFSETS]
    PREFIXES = ['half1_', 'half2_']
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'
        self.outputClasses = 'classes.xmd'
        self.outputVols = ''
        self.numberOfClasses.set(1) # For refinement we only need just one "class"
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _setSamplingArgs(self, args):
        """ Set sampling related params"""
        # Sampling stuff
        args['--healpix_order'] = self.angularSamplingDeg.get()
        args['--auto_local_healpix_order'] = self.localSearchAutoSamplingDeg.get()
        args['--auto_refine'] = ''
        args['--split_random_halves'] = ''
        #TODO: check why only for C*???
        # I have added by default for refine3d 
        # as extra parameters
        #if args['--sym'].startswith('C'):
        #    args['--low_resol_join_halves'] = "40";
        
    #--------------------------- STEPS functions --------------------------------------------     
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        vol = Volume()
        vol.setFileName(self._getExtraPath('relion_class001.mrc'))
        vol.setSamplingRate(imgSet.getSamplingRate())
        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(imgSet, vol)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validateNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _validateContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
    
    def _summaryNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
