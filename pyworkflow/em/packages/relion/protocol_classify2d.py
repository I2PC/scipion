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
This module contains the protocol base class for Relion protocols
"""

import os, re
from glob import glob

from protocol_base import *
from constants import MASK_FILL_ZERO



class ProtRelionClassify2D(ProtRelionBase, ProtClassify):
    """ This class implements the wrapper to Relion 2D - class averages program.
    """
    _label = '2d classify'
    IS_CLASSIFY = True
    IS_2D = True
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(rootDir)s/relion_it%(iter)03d_class%(ref)03d.mrc:mrc'
        self.outputClasses = 'classes.xmd'
        self.outputVols = ''
        

    #--------------------------- INSERT steps functions --------------------------------------------  
    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args.update({'--iter': self.numberOfIterations.get(),
                    '--tau2_fudge': self.regularisationParamT.get(),
                    '--flatten_solvent': '',
                    '--norm': '',
                    '--scale': '',
                    '--o': self._getExtraPath('relion'),
                    '--oversampling': '1',
                    '--j': self.numberOfThreads.get()
                    })
        self._setSamplingArgs(args)
        
    def _setCTFArgs(self, args):        
        # CTF stuff
        if self.doCTF:
            args['--ctf'] = ''
        
        if self.hasReferenceCTFCorrected:
            args['--ctf_corrected_ref'] = ''
            
        if self.haveDataPhaseFlipped:
            args['--ctf_phase_flipped'] = ''
            
        if self.ignoreCTFUntilFirstPeak:
            args['--ctf_intact_first_peak'] = ''
        
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        # Sampling stuff            
        args['--psi_step'] = self.inplaneAngularSamplingDeg.get() * 2
        args['--offset_range'] = self.offsetSearchRangePix.get()
        args['--offset_step']  = self.offsetSearchStepPix.get() * 2
        
    def _setMaskArgs(self, args):
        if self.referenceMask.hasValue():
            args['--solvent_mask'] = self.referenceMask.get().getFileName() #FIXE: CHANGE BY LOCATION
            
        if self.solventMask.hasValue():
            args['--solvent_mask2'] = self.solventMask.get().getFileName() #FIXME: CHANGE BY LOCATION
            
        
    def _setNormalArgs(self, args):
        args.update({'--i': self._getFileName('input_particles'),
                '--particle_diameter': self.maskRadiusA.get() * 2.0,
                '--angpix': self.inputParticles.get().getSamplingRate(),
                })
        self._setBasicArgs(args)
        self._setMaskArgs(args)
        
        if not self.isMapAbsoluteGreyScale:
            args[' --firstiter_cc'] = '' 
            
        if self.maskZero == MASK_FILL_ZERO:
            args['--zero_mask'] = ''
            
        args['--K'] = self.NumberOfClasses
        

    def _setContinueArgs(self):
        args = {}
        self.copyAttributes(self.continueRun.get(), 'regularisationParamT')
        self._setBasicArgs(args)
        args['--continue'] = self._getFileName('optimiser', 
                                               iter=self.continueIter.get())
        
    #--------------------------- STEPS functions --------------------------------------------       
    def createOutputStep(self):
        pass # should be implemented in subclasses
        
    
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
   
    


