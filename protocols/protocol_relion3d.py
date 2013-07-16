#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from os.path import join, exists
from os import remove
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from xmipp import MetaData, Image, MDL_IMAGE, MDL_ITER, MDL_LL, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, \
getBlocksInMetaDataFile, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDValueEQ, MDL_SAMPLINGRATE, MDL_CTF_MODEL
from xmipp import MetaDataInfo
from protlib_utils import runShowJ, getListFromVector, getListFromRangeString, runJob
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, cyanStr
from protlib_gui_ext import showWarning
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile, replaceBasenameExt
from protocol_ml2d import lastIteration
from protlib_filesystem import createLink
from protlib_import import exportReliontoMetadataFile

class ProtRelion3D(XmippProtocol):
    relionFiles=['data','model','optimiser','sampling']
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.relion3d.name, scriptname, project)
        self.Import = 'from protocol_relion3d import *'
        self.ParamsStr = ''
        if self.NumberOfMpi > 1:
            self.program = 'relion_refine_mpi'
        else:
            self.program = 'relion_refine'
        self.ParamsDict['program'] = self.program
        
        self.addParam('ORoot', self.WorkingDir + '/')        
        acquisionInfo = self.findAcquisitionInfo(self.ImgMd)

        if not acquisionInfo is None: 
            md = MetaData(acquisionInfo)
            self.addParam('SamplingRate', md.getValue(MDL_SAMPLINGRATE, md.firstObject()))
        
    def summary(self):
        md = MetaData(self.ImgMd)
        lines = ["Wrapper implemented for RELION **1.2**"]
        lines.append("Input images:  [%s] (<%u>)" % (self.ImgMd, md.size()))
       
        return lines
    
    def validate(self):
        errors = []
        # TODO: Check if relion is installed
        md = MetaData(self.ImgMd)
        if md.containsLabel(MDL_IMAGE):
            # Check that have same size as images:
            from protlib_xmipp import validateInputSize
            if not xmippExists(self.Ref3D):
               errors.append("Reference: <%s> doesn't exists" % ref)
            if len(errors) == 0:
                validateInputSize([self.Ref3D], self.ImgMd, md, errors)
        else:
            errors.append("Input metadata <%s> doesn't contains image label" % self.ImgMd)
            
        if self.DoCTFCorrection and not md.containsLabel(MDL_CTF_MODEL):
            errors.append("CTF correction selected and input metadata <%s> doesn't contains CTF information" % self.ImgMd)
            
        return errors 
    
    def defineSteps(self): 
        self.doContinue = len(self.ContinueFrom) > 0
        
        restart = False
        #temporary normalized files
        tmpFileNameSTK = join(self.ExtraDir, 'norRelion.stk')
        tmpFileNameXMD = join(self.ExtraDir, 'norRelion.xmd')
        if restart:            #Not yet implemented
            print 'NOT YET IMPLEMENTED'
            pass
        else:
            #create extra output directory
            self.insertStep('createDir', verifyfiles=[self.ExtraDir], path=self.ExtraDir)
            #normalize data

            if self.DoNormalizeInputImage:
                self.insertStep('runNormalizeRelion', verifyfiles=[tmpFileNameSTK,tmpFileNameXMD],
                                inputMd  = self.ImgMd,
                                outputMd = tmpFileNameSTK,
                                normType = 'NewXmipp',
                                bgRadius = int(self.MaskDiameterA/2.0/self.SamplingRate),
                                Nproc    = self.NumberOfMpi
                                )
            else:
                self.Db.insertStep('createLink',
                                   verifyfiles=[tmpFileName],
                                   source=self.ImgMd,
                                   dest=tmpFileNameXMD)
            # convert input metadata to relion model
            self.ImgStar = self.extraPath(replaceBasenameExt(tmpFileNameXMD, '.star'))
            self.insertStep('convertImagesMd', verifyfiles=[self.ImgStar],
                            inputMd=tmpFileNameXMD, 
                            outputRelion=self.ImgStar                            
                            )
            # launch relion program
            self.insertRelionRefine()
            # convert relion output to xmipp
            # relion_it00N_data.star, angular assigment
            self.ImgStar = self.extraPath(replaceBasenameExt(tmpFileNameXMD, '.star'))
            verifyFiles=[]
            inputs=[]
            for i in range (0,self.NumberOfIterations):
                for v in self.relionFiles:
                     verifyFiles += [self.getFilename(v+'Xm', iter=i )]
                for v in self.relionFiles:
                     inputs += [self.getFilename(v+'Re', iter=i )]
                self.insertStep('convertRelionMetadata', verifyfiles=verifyFiles,
                                inputs  = inputs,
                                outputs = verifyFiles
                                )
            
    def insertRelionRefine(self):
        args = {'--iter': self.NumberOfIterations,
                '--tau2_fudge': self.RegularisationParamT,
                '--flatten_solvent': '',
                '--zero_mask': '',# this is an option but is almost always true
                '--norm': '',
                '--scale': '',
                '--o': '%s/relion' % self.ExtraDir
                }
        if len(self.ReferenceMask):
            args['--solvent_mask'] = self.ReferenceMask
            
        if self.doContinue:
            args['--continue'] = self.ContinueFrom
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
            
            args['--K'] = self.NumberOfClasses
            
        # Sampling stuff
        # Find the index(starting at 0) of the selected
        # sampling rate, as used in relion program
        iover = 1 #TODO: check this DROP THIS
        index = ['30','15','7.5','3.7','1.8',
                 '0.9','0.5','0.2','0.1'].index(self.AngularSamplingDeg)
        args['--healpix_order'] = float(index + 1 - iover)
        
        if self.PerformLocalAngularSearch:
            args['--sigma_ang'] = self.LocalAngularSearchRange / 3.
            
        args['--offset_range'] = self.OffsetSearchRangePix
        args['--offset_step']  = self.OffsetSearchStepPix * pow(2, iover)

        args['--j'] = self.NumberOfThreads
        
        # Join in a single line all key, value pairs of the args dict    
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])
        params += self.AdditionalArguments
        verifyFiles=[]
        #relionFiles=['data','model','optimiser','sampling']
        for v in self.relionFiles:
             verifyFiles += [self.getFilename(v+'Re', iter=self.NumberOfIterations )]
#        f = open('/tmp/myfile','w')
#        for item in verifyFiles:
#            f.write("%s\n" % item)
#        f.close
        self.insertRunJobStep(self.program, params,verifyFiles)
        
            
    def createFilenameTemplates(self):
        extra = self.workingDirPath('extra')
        extraIter = join(extra, 'relion_it%(iter)03d_')
        myDict={}
        for v in self.relionFiles:
            myDict[v+'Re']=extraIter + v +'.star'
            myDict[v+'Xm']=extraIter + v +'.xmd'
        return myDict
#            
#
#def runRelion3D(log, program, params, mpi, threads):
#    print "program: ", program
#    print "params: ", params
#    runJob(log, program, params, mpi, threads)
    
      
        
def convertImagesMd(log, inputMd, outputRelion):
    """ Convert the Xmipp style MetaData to one ready for Relion.
    Main differences are: STAR labels are named different and
    in Xmipp the images rows contains a path to the CTFModel, 
    while Relion expect the explict values in the row.
    Params:
     input: input filename with the Xmipp images metadata
     output: output filename for the Relion metadata
    """
    from protlib_import import exportMdToRelion
    
    md = MetaData(inputMd)
    # Get the values (defocus, magnification, etc) from each 
    # ctfparam files and put values in the row
    md.fillExpand(MDL_CTF_MODEL)
    # Create the mapping between relion labels and xmipp labels
    exportMdToRelion(md, outputRelion)
    
def convertRelionMetadata(log, inputs,outputs):
    """ Convert the relion style MetaData to one ready for xmipp.
    Main differences are: STAR labels are named different and
    optimiser.star -> changes in orientation, offset. number images assigned to each class
    data.star -> loglikelihood (per image) may be used to delete worst images (10%)
                 orientation.shift per particle
    model.star:average_P_max (plot per iteration)
               block: model_classes
                  class distribution, number of particles per class
                  estimated error in orientation and translation
               block: model_class_N: 
                    resolution-dependent SNR: report resol where it drops below 1? (not that important?)
                                    (only in auto-refine) Gold-std FSC: make plot!
    """
    for i,o in zip(inputs,outputs):
        exportReliontoMetadataFile(i,o)
            
def renameOutput(log, WorkingDir, ProgId):
    ''' Remove ml2d prefix from:
        ml2dclasses.stk, ml2dclasses.xmd and ml2dimages.xmd'''
    prefix = '%s2d' % ProgId
    for f in ['%sclasses.stk', '%sclasses.xmd', '%simages.xmd']:
        f = join(WorkingDir, f % prefix)
        nf = f.replace(prefix, '')
        moveFile(log, f, nf)
                
def runNormalizeRelion(log,inputMd,outputMd,normType,bgRadius,Nproc,):
    stack = inputMd
    if exists(outputMd):
        remove(outputMd)
    program = "xmipp_transform_normalize"
    args = "-i %(stack)s -o %(outputMd)s "

    if bgRadius <= 0:
        particleSize = MetaDataInfo(stack)[0]
        bgRadius = int(particleSize/2)

    if normType=="OldXmipp":
        args += "--method OldXmipp"
    elif normType=="NewXmipp":
        args += "--method NewXmipp --background circle %(bgRadius)d"
    else:
        args += "--method Ramp --background circle %(bgRadius)d"
    runJob(log, program, args % locals(), Nproc)

