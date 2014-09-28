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

from pyworkflow.protocol.params import (PointerParam, FloatParam, PathParam,
                                        BooleanParam, IntParam, LEVEL_EXPERT)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D

from protocol_base import ProtRelionBase


class ProtRelionPostprocess(ProtAnalysis3D, ProtRelionBase):
    """    
    Reconstruct a volume using Relion from a given set of particles.
    The alignment parameters will be converted to a Relion star file
    and used as direction projections to reconstruct.
    """
    _label = 'post-process'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        
        form.addSection(label='Input')
        form.addParam('protRelionRefine', PointerParam, pointerClass="ProtRelionRefine3D",
                      label='Select a previous relion refine3D protocol',
                      help='Select a previous relion refine3D run to get the 3D half maps.')
        
        form.addSection(label='Masking')
        form.addParam('doAutoMask', BooleanParam, default=True,
                       label='Perform automated masking?',
                       help='Perform automated masking, based on a density threshold')
        form.addParam('initMaskThreshold', FloatParam, default=0.02,
                      condition='doAutoMask',
                      label='Initial binarisation threshold',
                      help='This threshold is used to make an initial binary mask from the'
                           "average of the two unfiltered half-reconstructions. If you don't"
                           "know what value to use, display one of the unfiltered half-maps in a 3D"
                           "surface rendering viewer, like Chimera, and find the lowest threshold"
                           " that gives no noise peaks outside the reconstruction.")
        form.addParam('extendInitMask', IntParam, default=3,
                      label='Mask pixels extension (px)', condition='doAutoMask',
                      help='The initial binary mask is extended this number of pixels in all directions.') 
        form.addParam('addMaskEdge', IntParam, default=3,
                      label='add soft-edge width (px)', condition='doAutoMask',
                      help='The extended binary mask is further extended with a raised-cosine soft'
                           'edge of the specified width.')
        form.addParam('mask', PointerParam, pointerClass='Mask',
                      label='Provide a mask', allowsNull=True,
                      condition='not doAutoMask',
                      help='Use this to skip auto-masking by providing your own mask') 
        
        form.addSection(label='Sharpening')
        form.addParam('mtf', PathParam,
                      label='MTF-curve file',
                      help='User-provided STAR-file with the MTF-curve of the detector.'
                           'Relion param: <--mtf>')
        form.addParam('doAutoBfactor', BooleanParam, default=True,
                       label='Estimate B-factor automatically?',
                       help='If set to Yes, then the program will use the automated procedure'
                            'described by Rosenthal and Henderson (2003, JMB) to estimate an overall'
                            'B-factor for your map, and sharpen it accordingly.')
        line = form.addLine('B-factor resolution (A): ', condition='doAutoBfactor',
                            help='There are the frequency (in Angstroms), lowest and highest,'
                            'that will be included in the linear fit of the Guinier plot as described'
                            'in Rosenthal and Henderson (2003, JMB).')
        line.addParam('bfactorLowRes', FloatParam, default='10.0', label='low ')
        line.addParam('bfactorHighRes', FloatParam, default='0.0', label='high ')            
        form.addParam('bfactor', FloatParam, default=-350,
                      condition='not doAutoBfactor',
                      label='Provide B-factor:',
                      help= 'User-provided B-factor (in A^2) for map sharpening, e.g. -400.'
                            'Use negative values for sharpening. Be careful: if you over-sharpen\n'
                            'your map, you may end up interpreting noise for signal!\n'
                            'Relion param: *--adhoc_bfac*')
        
        form.addSection(label='Filtering')
        form.addParam('skipFscWeighting', BooleanParam, default=False,
                       label='Skip FSC-weighting for sharpening?',
                       help='If set to No (the default), then the output map will be low-pass'
                       ' filtered according to the mask-corrected, gold-standard FSC-curve.'
                       ' Sometimes, it is also useful to provide an ad-hoc low-pass filter'
                       ' (option below), as due to local resolution variations some parts of'
                       ' the map may be better and other parts may be worse than the overall'
                       ' resolution as measured by the FSC. In such cases, set this option to'
                       ' Yes and provide an ad-hoc filter as described below.')
        form.addParam('lowRes', FloatParam, default=5,
                      condition='skipFscWeighting',
                      label='Low-pass filter (A):',
                      help='This option allows one to low-pass filter the map at a user-provided'
                      ' frequency (in Angstroms). When using a resolution that is higher than the'
                      ' gold-standard FSC-reported resolution, take care not to interpret noise'
                      ' in the map for signal...')
        form.addParam('filterEdgeWidth', IntParam, default=2, expertLevel=LEVEL_EXPERT,
                      label='Low-pass filter edge width:',
                      help='Width of the raised cosine on the low-pass filter edge'
                           '(in resolution shells)\n'
                           'Relion param: *--filter_edge_width*')
        form.addParam('randomizeAtFsc', FloatParam, default=5, expertLevel=LEVEL_EXPERT,
                      label='Randomize phases threshold',
                      help='Randomize phases from the resolution where FSC drops below this value\n'
                           'Relion param: *--randomize_at_fsc*')
        
        form.addParallelSection(threads=0, mpi=0) 
            
    #--------------------------- INSERT steps functions --------------------------------------------  

    def _insertAllSteps(self):
        self._initialize()
        self._insertPostProcessingStep()
        self._insertFunctionStep('createOutputStep')
        
        
    def _initialize(self):
        # ToDo: implement a better way to get the pattern of unmasked maps.
        self.input = self.protRelionRefine.get()._getExtraPath('relion')
        print "INPUT: ", self.input
        
    def _insertPostProcessingStep(self):
        
        output = self._getExtraPath('postprocess')
        self.samplingRate = self.protRelionRefine.get().inputParticles.get().getSamplingRate()
        
        args = " --i %s --o %s --angpix %f" %(self.input, output, self.samplingRate)

        if self.doAutoMask:
            args += " --auto_mask --inimask_threshold %f" % self.initMaskThreshold.get()
            args += " --extend_inimask %d" % self.extendInitMask.get()
            args += " --width_mask_edge %d" % self.addMaskEdge.get()
        else:
            args += ' --mask %s' % self.mask.get().getFileName()
            
        mtfFile = self.mtf.get()
        if mtfFile:
            args += ' --mtf %s' % mtfFile
            
        if self.doAutoBfactor:
            args += ' --auto_bfac --autob_lowres %f' % self.bfactorLowRes.get()
            args += ' --autob_highres %f' % self.bfactorHighRes.get()
        else:
            args += ' --adhoc_bfac %f' % self.bfactor.get()
            
        if self.skipFscWeighting:
            args += ' --skip_fsc_weighting --low_pass %f' % self.lowRes.get()
            
        # Expert params
        args += ' --filter_edge_width %d' % self.filterEdgeWidth.get()
        args += ' --randomize_at_fsc %f' % self.randomizeAtFsc.get()
        
        self._insertFunctionStep('postProcessStep', args)
        
    #--------------------------- STEPS functions --------------------------------------------
    def postProcessStep(self, params):
        
        self.runJob('relion_postprocess', params)
        
    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getExtraPath('postprocess.mrc'))
        volume.setSamplingRate(self.samplingRate)
        vol = self.protRelionRefine.get().outputVolume
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(vol, volume)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
