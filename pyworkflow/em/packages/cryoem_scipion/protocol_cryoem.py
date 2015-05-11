# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains wrapper around cryoem algorithm
"""

import os
from glob import glob

from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.protocol.params as params
from pyworkflow.utils.path import cleanPath, cleanPattern

import pyworkflow.em as em
from pyworkflow.em.protocol import ProtInitialVolume
import cryoem_scipion


class ProtCryoem(ProtInitialVolume):
    """ Produces one or several initial volumes using cryoem pseudoatoms """
    _label = 'pseudoatoms'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputClasses', params.PointerParam, label="Input classes", 
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      help='Select the input classes2D from the project.\n'
                           'It should be a SetOfClasses2D class with class representative')
        form.addParam('resize', params.IntParam, default=32, expertLevel=LEVEL_ADVANCED,
                      label='Resize the input images',
                      help='In Angstrom.'
                      'Resize the input images to an image dxd in Angstrom.')
        form.addParam('stepscoarse', params.IntParam, default=2000, expertLevel=LEVEL_ADVANCED,
                      label='Number of steps for coarsing',
                      help=''
                      'It determinines the total number of sampling steps for the coarse phase.')
        form.addParam('stepsrefine', params.IntParam, default=1000, expertLevel=LEVEL_ADVANCED,
                      label='Number of steps for refining',
                      help=''
                      'It determinines the total number of sampling steps for the refine phase.')
        form.addParam('rotations', params.IntParam, default=5000, expertLevel=LEVEL_ADVANCED,
                      label='Number of random rotations',
                      help=''
                      'It determinines the number of random rotations.')
        form.addParam('counts', params.IntParam, default=10000, expertLevel=LEVEL_ADVANCED,
                      label='Counts per image',
                      help=''
                      'It determinines the number of counts per image.')
        form.addParam('mixture', params.IntParam, default=500, expertLevel=LEVEL_ADVANCED,
                      label='Components mixture model',
                      help=''
                      'It determinines the number of components in the mixture model.')
        
        
        
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for calling simple prime program"""
        
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runCryoemPseudoCoarse')
        self._insertFunctionStep('runCryoemPseudoRefine')
        self._insertFunctionStep('createOutputStep')        

    #--------------------------- STEPS functions --------------------------------------------        
    def convertInputStep(self):
        self.inputClasses.get().writeStack(self._getExtraPath("classes.mrc:mrcs"))
               
    def runCryoemPseudoCoarse(self):
        inputClasses = self.inputClasses.get()
        resz = self.resize.get()
        stps_coar = self.stepscoarse.get()
        rota = self.rotations.get()
                
        xdim, _, _ = inputClasses.getDimensions()
        Ts = float(inputClasses.getSamplingRate())
        args="%s/gmm-rec.py %s -v -o %s -s %d 50 50 -g -r %d -d %d -a %f -A %f"%(cryoem_scipion.CRYOEM_BIN,
                                                                       self._getExtraPath("classes.mrc"),\
                                                                       self._getExtraPath("coarseVolume.mrc"),\
                                                                       stps_coar,\
                                                                       rota,\
                                                                       resz,\
                                                                       Ts,Ts)
        
        #args="%s/gmm-rec.py %s -v -o %s -s 2000 50 50 -g -r 5000 -N 2000 -d 32 -K 200 -a %f -A %f"%(cryoem_scipion.CRYOEM_BIN,
        #                                                               self._getExtraPath("classes.mrc"),\
        #                                                               self._getExtraPath("coarseVolume.mrc"),\
        #                                                               Ts,Ts)
        
        #args="%s/gmm-rec.py %s -v -o %s -s 5000 100 100 -g -r 10000 -d 32 -N 2000 -K 200 -a 2.82 -A 2 -p 8"%(cryoem_scipion.CRYOEM_BIN,
        #                                                               self._getExtraPath("classes.mrc"),\
        #                                                               self._getExtraPath("coarseVolume.mrc"))
        #-r number of random rotations
        self.runJob("python", args)
        print 'ZERO APROXIMATION - ENDED'
        
 
    def runCryoemPseudoRefine(self):
        inputClasses = self.inputClasses.get()
        stps_ref = self.stepsrefine.get()
        cou = self.counts.get()
        mix = self.mixture.get()
        
        xdim, _, _ = inputClasses.getDimensions()

        Ts = float(inputClasses.getSamplingRate())
       
        args2 = "%s/gmm-rec.py %s -v -o %s -s %d -K %d -t %s -N %d -a %f -A %f"%(cryoem_scipion.CRYOEM_BIN,
                                                                                     self._getExtraPath("classes.mrc"),\
                                                                                     self._getExtraPath("refine.mrc"),\
                                                                                     stps_ref,\
                                                                                     mix,\
                                                                                     self._getExtraPath("coarseVolume_transformations.csv"),\
                                                                                     cou,\
                                                                                     Ts,Ts)
        #The number of pseudoatoms -K.
        #The output files *csv can be used as initial mixtures or transformations for subsequent reconstructions using 
                                                                                                        #**-m** or **-t**.
        #The number of counts per image is specified by **-N**.
        #The input pixel size of the original images is specified in Angstrom using **-a**.
        #After reconstruction, the final mixture is evaluated on a 3D grid with grid spacing given by **-A** in Angstrom, 
                                                                                       #to be viewed in chimera.
        self.runJob("python", args2)
        print 'REFINEMENT - ENDED'
        
        
    def createOutputStep(self):
        inputClasses = self.inputClasses.get()
        imgxdim, _, _ = inputClasses.getDimensions()
        fnVol = self._getExtraPath("refine.mrc")
        volxdim, _, _, _ = em.ImageHandler().getDimensions(fnVol)
        
        if imgxdim != volxdim:
            fnVolOut = self._getPath("volume.vol")
            args3 = "-i %s --dim %d -o %s"%(fnVol, xdim, fnVolOut)
            #TODO maybe change later using Python binding to resize images
            import pyworkflow.em.packages.xmipp3 as xmipp3
            self.runJob("xmipp_image_resize", args3, env=xmipp3.getEnviron())
        else:
            fnVolOut = fnVol
        
        vol = em.Volume()
        vol.setLocation(fnVolOut)
        vol.setSamplingRate(inputClasses.getSamplingRate())
        self._defineOutputs(outputVol=vol)
       
        self._defineSourceRelation(vol, self.inputClasses)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Input classes: %s" % self.getObjectTag('inputClasses'))
        return summary
    
    def _citations(self):
        return ['Joubert2015']
    
    def _methods(self):
        if self.inputClasses.get() is not None:
            retval="We used *gmm-rec* program [Joubert2015] to produce an initial volume from the set of classes %s."
            return [retval % self.getObjectTag('inputClasses')]
        else:
            return []
