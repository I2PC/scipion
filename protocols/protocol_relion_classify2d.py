#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Relion-based 2D classification (relion 1.2)
#
# Author: J. M. de la Rosa Trevin     (jmdelarosa@cnb.csic.es) Feb 2014
#

from os.path import join, exists
from os import remove, rename
from protlib_base import protocolMain
from config_protocols import protDict
from xmipp import *
from xmipp import MetaDataInfo
from protlib_utils import runShowJ, getListFromVector, getListFromRangeString, \
                          runJob, runChimera, which, runChimeraClient
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, cyanStr
from protlib_gui_ext import showWarning, showTable, showError
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile, \
                               replaceBasenameExt
from protocol_ml2d import lastIteration
from protlib_filesystem import createLink
from protlib_relion import *


class ProtRelion2DClassifier(ProtRelionBase):
    
    def __init__(self, scriptname, project):
        ProtRelionBase.__init__(self, protDict.relion_classify2d.name, scriptname, project)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'
        self.outputClasses = 'classes.xmd'
        self.outputVols = ''
        
        self.Import = 'from protocol_relion_classify2d import *'
        
        if self.DoContinue:
            #if optimizer has not been properly selected this will 
            #fail, let us go ahead and handle the situation in verify
            try:
                self.inputProperty('MaskRadiusA')
                self.inputProperty('NumberOfClasses')
                self.inputProperty('RegularisationParamT')

            except:
                print "Can not access the parameters from the original relion run"
        
    def _summary(self):
        lastIteration = self.lastIter()
        lines = ['Performed <%d/%d> iterations ' % (lastIteration, self.NumberOfIterations)]
        lines.append('Number of classes = <%d>' % self.NumberOfClasses)
        
        return lines

    def _summaryContinue(self):
        lastIteration = self.lastIter()
        firstIteration = getIteration(self.optimiserFileName)
        lines = ['Continuation from run: <%s>, iter: <%d>' % (self.PrevRunName, firstIteration)]
        if (lastIteration - firstIteration) < 0 :
            performedIteration = 0
        else:
            performedIteration = lastIteration - firstIteration
        lines.append('Performed <%d> iterations (number estimated from the files in working directory)' % performedIteration) 
        lines.append('Input fileName = <%s>'%self.optimiserFileName)
        return lines
    
    def _validateContinue(self):
        lastIterationPrecRun = self.PrevRun.lastIter()
        errors = []
        if '2D/RelionClass2D' not in self.optimiserFileName:
            errors.append('Please select an optimiser file produced by Relion 2D classification')
        return errors 
    
    def _insertSteps(self):
        args = {'--iter': self.NumberOfIterations,
                '--tau2_fudge': self.RegularisationParamT,
                '--flatten_solvent': '',
                #'--zero_mask': '',# set it to true using additional argument
                '--norm': '',
                '--scale': '',
                '--o': '%s/relion' % self.ExtraDir
                }
        if getattr(self, 'ReferenceMask', '').strip():
            args['--solvent_mask'] = self.ReferenceMask.strip()
            
        if getattr(self, 'SolventMask', '').strip():
            args['--solvent_mask2'] = self.SolventMask.strip()
            
        args.update({'--i': self.ImgStar,
                     '--particle_diameter': self.MaskRadiusA * 2.0,
                     '--angpix': self.SamplingRate,
                     '--oversampling': '1'
                     })
        
        if not self.IsMapAbsoluteGreyScale:
            args[' --firstiter_cc'] = '' 
            
        # CTF stuff
        if self.DoCTFCorrection:
            args['--ctf'] = ''
        
        if self.HasReferenceCTFCorrected:
            args['--ctf_corrected_ref'] = ''
            
        if self.HaveDataPhaseFlipped:
            args['--ctf_phase_flipped'] = ''
            
        if self.IgnoreCTFUntilFirstPeak:
            args['--ctf_intact_first_peak'] = ''
            
        MaskZero = getattr(self, 'MaskZero', '')
        if MaskZero.startswith('Yes'):
            args['--zero_mask'] = ''
            
        args['--K'] = self.NumberOfClasses
            
        # Sampling stuff
        args['--psi_step'] = self.InPlaneAngularSamplingDeg
            
        iover = 1 #TODO: check this DROP THIS
        args['--offset_range'] = self.OffsetSearchRangePix
        args['--offset_step']  = self.OffsetSearchStepPix * pow(2, iover)

        args['--j'] = self.NumberOfThreads
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        params += ' ' + self.AdditionalArguments

        self.insertRunJobStep(self.program, params, self._getIterFiles(self.NumberOfIterations))

    def _insertStepsContinue(self):
        args = {
                '--o': '%s/relion' % self.ExtraDir,
                '--continue': self.optimiserFileName,
                
                '--iter': self.NumberOfIterations,
                
                '--tau2_fudge': self.RegularisationParamT,# should not be changed 
                '--flatten_solvent': '',# use always
                #'--zero_mask': '',# use always. This is confussing since Sjors gui creates the command line
                #                  # with this option
                #                  # but then the program complains about it. 
                '--oversampling': '1',
                '--norm': '',
                '--scale': '',
                }
        iover = 1 #TODO: check this DROP THIS
        
         # Sampling stuff
        args['--psi_step'] = self.InPlaneAngularSamplingDeg * pow(2, iover)        
        args['--offset_range'] = self.OffsetSearchRangePix
        args['--offset_step']  = self.OffsetSearchStepPix * pow(2, iover)

        args['--j'] = self.NumberOfThreads
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        params += ' ' + self.AdditionalArguments

        self.insertRunJobStep(self.program, params, self._getIterFiles(self.NumberOfIterations))

    def _getExtraOutputs(self):
        return []
    
    def _getFinalSuffix(self):
        return ''       

    def _getPrefixes(self):
        """ This should be redefined in refine protocol. """
        return [''] 

    def _getChangeLabels(self):
        return [MDL_AVG_CHANGES_ORIENTATIONS, MDL_AVG_CHANGES_OFFSETS, MDL_AVG_CHANGES_CLASSES]  

    def _visualizeDisplayImagesClassification(self):
        """ Read Relion _data.star images file and 
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        for it in self._visualizeIterations:
            data_classes = self._getIterClasses(it)
            self.display2D(data_classes, extraParams='--mode metadata --render first')
            
    def _checkIterData(self):
        """ check that for last visualized iteration, the data is produced. """
        dataStar = self.getFilename('data', iter=self._visualizeLastIteration)
        if not xmippExists(dataStar):
            message = "No data available for <iteration %d>, file <%s> not found." % (self._visualizeLastIteration, dataStar)
            showError("File not Found", message, self.master)
            return False
        return True
            
             