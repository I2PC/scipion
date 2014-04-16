#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for wrapping the call to relion_postprocess program
#
#   Author:  J. M. de la Rosa Trevin   (jmdelarosa@cnb.csic.es)  Feb 2014
#------------------------------------------------------------------------------------------------

from os.path import dirname, join
import xmipp
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from protlib_utils import runShowJ
from protlib_xmipp import getSampling

VOLUMES = 'volumes.xmd'
                      
class ProtRelionPostProcess(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.relion_postprocess.name, scriptname, project)
        self.Import = 'from protocol_relion_postprocess import *'

    def summary(self):
        lines = []
        return lines
    
    def papers(self):
        return ['Bayesian view: Scheres, JMB (2012) [http://www.ncbi.nlm.nih.gov/pubmed/22100448]',
                'RELION implementation: Scheres, JSB (2012) [http://www.ncbi.nlm.nih.gov/pubmed/23000701]'
                  ]

    def validate(self):
        errors = []
        
        return errors 
    
    def defineSteps(self):
        samplingRate = getSampling(self.ModelStar)  
        args = ""
        args += " --i %s/relion" % dirname(self.ModelStar)
        args += " --angpix %f" % samplingRate
        args += " --o %s/relion_postprocess" % self.WorkingDir
        
        if self.AutoMask:
            args += ' --auto_mask'
        else:
            args += ' --inimask_threshold %f' % self.InimaskThreshold
            args += ' --extend_inimask %f' % self.ExtendInimask
            args += ' --width_mask_edge %f' % self.WidthMaskEdge
            args += ' --mask %s' % self.Mask
            
        mtfFile = self.Mtf.strip()
        if mtfFile:
            args += ' --mtf %s' % mtfFile
            
        if self.AutoBfac:
            args += ' --auto_bfac'
            args += ' --autob_lowres %f' % self.AutobLowres
            args += ' --autob_highres %f' % self.AutobHighres
        else:
            args += ' --adhoc_bfac %f' % self.AdhocBfac
            
        if not self.UseFscWeighting:
            args += ' --skip_fsc_weighting'
            
        if self.LowPass > 0.:
            args += ' --low_pass %f' % self.LowPass
            
        # Expert params
        args += ' --filter_edge_width %d' % self.FilterEdgeWidth
        args += ' --randomize_at_fsc %f' % self.RandomizeAtFsc
        args += ' --verb %d' % self.Verb   
        
        outputs = [self.workingDirPath('relion_postprocess%s.mrc' % s) 
                   for s in ['', '_masked', '_automask']]
        self.insertRunJobStep('relion_postprocess', args, verifyFiles=outputs)
        self.insertStep('createVolumesMd', 
                        verifyfiles=[self.workingDirPath(VOLUMES)],
                        volumeList=outputs)
        
    def visualize(self):
        runShowJ(self.workingDirPath(VOLUMES))
    
def createVolumesMd(log, volumeList):
    md = xmipp.MetaData()
    
    for v in volumeList:
        md.setValue(xmipp.MDL_IMAGE, '%s:mrc' % v,
                    md.addObject())
        
    root = dirname(volumeList[0])
    md.write('volumes@' + join(root, 'volumes.xmd'))
    
    
