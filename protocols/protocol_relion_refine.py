#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#Wrapper for relion 1.2

#   Author: Roberto Marabini
#

from os.path import join, exists
from os import remove, rename, environ
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
from protlib_gui_figure import XmippPlotter
from protlib_gui_ext import showWarning, showTable, showError

from protocol_relion_base import ProtRelionBase, runNormalizeRelion, convertImagesMd, renameOutput, \
                                 convertRelionBinaryData, convertRelionMetadata

class ProtRelionRefinner( ProtRelionBase):
    def __init__(self, scriptname, project):
        ProtRelionBase.__init__(self, protDict.relion_refine.name, scriptname, project)
        self.Import = 'from protocol_relion_refine import *'
        self.relionType='refine'

    def summary(self):
        lines = ProtRelionBase.summary(self)
        lastIteration=self.lastIter()
        lines += ['Performed <%d> iterations ' % lastIteration]
        return lines
    
    def validate(self):
        errors = ProtRelionBase.validate(self)
        if self.NumberOfMpi < 3:
            errors.append('''relion refine requires at least 3 mpi proesses to compute golden standard''')
        return errors 
    
    def defineSteps(self):
        #create new directories, normalize and convert metadata
        ProtRelionBase.defineSteps(self)
        # launch relion program
        self.insertRelionRefine()
        #horrible hack since I do not know a priory the number of iterations
        #I just check that the 0 and 1 are done (0 -> is the input after filtration
        lastIteration   = 1
        NumberOfClasses = 1
        extra = self.workingDirPath('extra')
        ExtraInputs  = [ join(extra,'relion_model.star')] 
        ExtraOutputs = [ join(extra,'relion_model.xmd')] 
        ProtRelionBase.defineSteps2(self,lastIteration,NumberOfClasses,ExtraInputs,ExtraOutputs)
                       
    def insertRelionRefine(self):
        args = {#'--iter': self.NumberOfIterations,
                #'--tau2_fudge': self.RegularisationParamT,
                '--flatten_solvent': '',
                #'--zero_mask': '',# this is an option but is almost always true
                '--norm': '',
                '--scale': '',
                '--o': '%s/relion' % self.ExtraDir
                }
        if len(self.ReferenceMask):
            args['--solvent_mask'] = self.ReferenceMask
            
        if self.doContinue:
            args['--continue'] = self.ContinueFrom
            #note: no movie realigment is included, since 1) it is not the core of relion and 2) by design continue
            #does not allow to change parameters in xmipp.        
        else: # Not continue
            args.update({'--i': self.ImgStar,
                         '--particle_diameter': self.MaskDiameterA,
                         '--angpix': self.SamplingRate,
                         '--ref': self.Ref3D,
                         '--oversampling': '1'
                         })
            
            if not self.IsMapAbsoluteGreyScale:
                args[' --firstiter_cc'] = '' 
                
            if self.InitialLowPassFilterA > 0:
                args['--ini_high'] = self.InitialLowPassFilterA
                
            # CTF stuff
            if self.DoCTFCorrection:
                args['--ctf'] = ''
            
            if self.HasReferenceCTFCorrected:
                args['--ctf_corrected_ref'] = ''
                
            if self.HaveDataPhaseFlipped:
                args['--ctf_phase_flipped'] = ''
                
            if self.IgnoreCTFUntilFirstPeak:
                args['--ctf_intact_first_peak'] = ''
                
            args['--sym'] = self.SymmetryGroup.upper()
            args['--auto_refine']=''
            args['--split_random_halves']=''
            
            if args['--sym'][:1] == 'C':
                args['--low_resol_join_halves'] = "40";
            #args['--K'] = self.NumberOfClasses
            
        # Sampling stuff
        # Find the index(starting at 0) of the selected
        # sampling rate, as used in relion program
        iover = 1 #TODO: check this DROP THIS
        index = ['30','15','7.5','3.7','1.8',
                 '0.9','0.5','0.2','0.1'].index(self.AngularSamplingDeg)
        args['--healpix_order'] = float(index + 1 - iover)
        args['--auto_local_healpix_order'] = float(index + 1 - iover)
        
        #if self.PerformLocalAngularSearch:
        #    args['--sigma_ang'] = self.LocalAngularSearchRange / 3.
            
        args['--offset_range'] = self.OffsetSearchRangePix
        args['--offset_step']  = self.OffsetSearchStepPix * pow(2, iover)

        args['--j'] = self.NumberOfThreads
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        params += self.AdditionalArguments
        verifyFiles=[]
        #relionFiles=['data','model','optimiser','sampling']
        #refine does not predefine the number of iterations so no verify is possible
        #let us check that at least iter 1 is done
        for v in self.relionFiles:
             verifyFiles += [self.getFilename(v+'Re', iter=1, workingDir=self.WorkingDir )]
             #verifyFiles += [self.getFilename(v+'Re', iter=self.NumberOfIterations )]
#        f = open('/tmp/myfile','w')
#        for item in verifyFiles:
#            f.write("%s\n" % item)
#        f.close
        self.insertRunJobStep(self.program, params,verifyFiles)
        ###################self.insertRunJobStep('echo shortcut', params,verifyFiles)

    def createFilenameTemplates(self):
        myDict=ProtRelionBase.createFilenameTemplates(self)
        self.relionFiles += ['half1_model','half2_model']
        #relionFiles=['data','half1_model', 'half2_model','optimiser','sampling']
        for v in self.relionFiles:
            myDict[v+'Re']=self.extraIter + v +'.star'
            myDict[v+'Xm']=self.extraIter + v +'.xmd'
        #myDict['volumeh1']    = self.extraIter + "half1_class%(ref3d)03d.spi"
        #myDict['volumeMRCh1'] = self.extraIter + "half1_class%(ref3d)03d.mrc:mrc"
        myDict['volumeh1'] = self.extraIter + "half1_class%(ref3d)03d.mrc:mrc"
        
        #myDict['volumeh2']    = self.extraIter + "half2_class%(ref3d)03d.spi"
        #myDict['volumeMRCh2'] = self.extraIter + "half2_class%(ref3d)03d.mrc:mrc"
        myDict['volumeh2'] = self.extraIter + "half2_class%(ref3d)03d.mrc:mrc"
        
        extra = self.workingDirPath('extra')
        self.extraIter2 = join(extra, 'relion_')
        #myDict['volumeFinal']      = self.extraIter2 + "class%(ref3d)03d.spi"
        #myDict['volumeMRCFinal']   = self.extraIter2 + "class%(ref3d)03d.mrc:mrc"
        myDict['volumeFinal']   = self.extraIter2 + "class%(ref3d)03d.mrc:mrc"
        myDict['modelXmFinal']=self.ExtraDir+'/relion_model.xmd'
        myDict['imagesAssignedToClass']=''

        return myDict

               
