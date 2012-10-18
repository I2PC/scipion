#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, October 2011

#from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob, runShowJ
from protlib_filesystem import deleteFile, createLink2, exists, replaceFilenameExt
import xmipp
from protlib_gui_ext import showWarning

# The dictionary with specific filename templates 
# is defined here to allow use of it outside the protocol
_prefix = join('%(micrographDir)s','xmipp_ctf')
_templateDict = {
        # This templates are relative to a micrographDir
        'prefix': _prefix,
        'ctfparam': _prefix +  '.ctfparam',
        'psd': _prefix + '.psd',
        'enhanced_psd': _prefix + '_enhanced_psd.xmp',
        'ctfmodel_quadrant': _prefix + '_ctfmodel_quadrant.xmp',
        'ctfmodel_halfplane': _prefix + '_ctfmodel_halfplane.xmp'
        }

def _getFilename(key, **args):
    return _templateDict[key] % args

class ProtScreenMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_micrographs.name, scriptname, project)
        self.Import = "from protocol_screen_micrographs import *"
        self.setPreviousRun(self.ImportRun) 
        self.inputFilename('microscope', 'micrographs', 'acquisition')
        self.inputProperty('TiltPairs', 'MicrographsMd')
        self.micrographs = self.getFilename('micrographs')
        #TODO: check all the possible casses
        if not self.TiltPairs:
            self.MicrographsMd = self.micrographs

    def defineSteps(self):
        extraDir=self.workingDirPath('extra')
        parent_id = self.insertStep('createDir',verifyfiles=[extraDir],path=extraDir)

        filesToImport = [self.Input[k] for k in ['microscope', 'acquisition']]
        if self.TiltPairs:
            filesToImport.append(self.MicrographsMd)
        self.insertImportOfFiles(filesToImport)

        # Read Microscope parameters
        MD = xmipp.MetaData(self.Input['microscope'])
        objId = MD.firstObject()
        Voltage = MD.getValue(xmipp.MDL_CTF_VOLTAGE,objId)
        SphericalAberration = MD.getValue(xmipp.MDL_CTF_CS,objId)
        Magnification = MD.getValue(xmipp.MDL_MAGNIFICATION,objId)
        MD2 = xmipp.MetaData(self.Input['acquisition'])
        AngPix = MD2.getValue(xmipp.MDL_SAMPLINGRATE,objId)

        # Create verifyFiles for the MPI and output directories
        MD = xmipp.MetaData(self.Input['micrographs'])
        
        # Now the estimation actions
        for objId in MD:
            inputFile = MD.getValue(xmipp.MDL_MICROGRAPH, objId)
            micrographName = os.path.basename(inputFile)
            shortname = os.path.splitext(micrographName)[0]
            micrographDir = os.path.join(extraDir,shortname)                    
            parent_id = self.insertParallelStep('createDir',verifyfiles=[micrographDir],path=micrographDir,parent_step_id=XmippProjectDb.FIRST_STEP)

            # Downsample if necessary
            if self.DownsampleFactor != 1:
                finalname = self.tmpPath(shortname + "_tmp.mrc")
                parent_id = self.insertParallelRunJobStep("xmipp_transform_downsample",
                                                   "-i %s -o %s --step %f --method fourier" % (inputFile,finalname,self.DownsampleFactor),
                                                   [finalname],parent_step_id=parent_id)
            else:
                finalname = inputFile
            
            # CTF estimation with Xmipp
            args="--micrograph "+finalname+\
                 " --oroot " + _getFilename('prefix', micrographDir=micrographDir)+\
                 " --kV "+str(Voltage)+\
                 " --Cs "+str(SphericalAberration)+\
                 " --sampling_rate "+str(AngPix*self.DownsampleFactor)+\
                 " --downSamplingPerformed "+str(self.DownsampleFactor)+\
                 " --ctfmodelSize 256"+\
                 " --Q0 "+str(self.AmplitudeContrast)+\
                 " --min_freq "+str(self.LowResolCutoff)+\
                 " --max_freq "+str(self.HighResolCutoff)+\
                 " --pieceDim "+str(self.WinSize)+\
                 " --defocus_range "+str((self.MaxFocus-self.MinFocus)*10000/2)+\
                 " --defocusU "+str((self.MaxFocus+self.MinFocus)*10000/2)
            if (self.FastDefocus):
                args+=" --fastDefocus"
            self.insertParallelRunJobStep('xmipp_ctf_estimate_from_micrograph', args,
                                     verifyfiles=[_getFilename('ctfparam', micrographDir=micrographDir)],parent_step_id=parent_id)
        
        # Gather results after external actions
        self.insertStep('gatherResults',verifyfiles=[self.micrographs],
                           TmpDir=self.TmpDir,
                           WorkingDir=self.WorkingDir,
                           summaryFile=self.micrographs,
                           importMicrographs=self.Input['micrographs'],
                           Downsampling=self.DownsampleFactor,
                           NumberOfMpi=self.NumberOfMpi)
    
    def createFilenameTemplates(self):
        return _templateDict
    
    def validate(self):
        errors = []
        if self.DownsampleFactor<1:
            errors.append("Downsampling must be >=1");

        # Check that there are any micrograph to process
        if xmippExists(self.Input['micrographs']):
            md = xmipp.MetaData(self.Input['micrographs'])
            if md.isEmpty():
                errors.append("Imported micrographs file <%(micrographs)s> is empty" % self.Input)
        
        if self.AmplitudeContrast < 0:
            errors.append("Q0 should be positive")
        
        if self.MinFocus < 0 or self.MaxFocus<0:
            errors.append("Defoci range must be positive (minFocus, maxFocus)>0")
            
        if self.MaxFocus < self.MinFocus:
            errors.append("maxFocus must be larger than minFocus")
        
        return errors

    def summary(self):
        message = []
        md = xmipp.MetaData(self.Input['micrographs'])
        message.append("CTF screening of <%d> micrographs." % md.size())
        message.append("Input directory: [%s]" % self.PrevRun.WorkingDir)
        return message
    
    def visualize(self):
        summaryFile = self.getFilename('micrographs')
        
        if not exists(summaryFile): # Try to create partial summary file
            summaryFile = summaryFile.replace(self.WorkingDir, self.TmpDir)
            buildSummaryMetadata(self.WorkingDir, self.Input['micrographs'], summaryFile)
            
        if exists(summaryFile):
            runShowJ(summaryFile, extraParams = "--mode metadata")
        else:
            showWarning('Warning', 'There are not results yet',self.master)
    
def gatherResults(log,TmpDir,WorkingDir,summaryFile, importMicrographs,Downsampling,NumberOfMpi):
    buildSummaryMetadata(WorkingDir, importMicrographs, summaryFile)
    dirSummary,fnSummary=os.path.split(summaryFile)
    runJob(log,"xmipp_ctf_sort_psds","-i %s -o %s/aux_%s"%(summaryFile,dirSummary,fnSummary),NumberOfMpi=NumberOfMpi)
    runJob(log,"mv","-f %s/aux_%s %s"%(dirSummary,fnSummary,summaryFile))
    if Downsampling!=1:
        runJob(log,"rm","-f "+TmpDir+"/*")

def buildSummaryMetadata(WorkingDir,importMicrographs,summaryFile):
    md = xmipp.MetaData()
    importMd = xmipp.MetaData(importMicrographs)
    for id in importMd:
        inputFile = importMd.getValue(xmipp.MDL_MICROGRAPH,id)
        micrographName = os.path.basename(inputFile)
        shortname = replaceFilenameExt(micrographName, '')
        micrographDir = join(WorkingDir, 'extra', shortname)                    
        
        objId = md.addObject()
        md.setValue(xmipp.MDL_MICROGRAPH, inputFile, objId)
        ctfparam = _getFilename('ctfparam', micrographDir=micrographDir)
        labels = [xmipp.MDL_PSD, xmipp.MDL_PSD_ENHANCED, xmipp.MDL_CTF_MODEL,xmipp.MDL_IMAGE1, xmipp.MDL_IMAGE2]
        if exists(ctfparam): # Get filenames
            keys = ['psd', 'enhanced_psd', 'ctfparam', 'ctfmodel_quadrant', 'ctfmodel_halfplane']
            values = [_getFilename(key, micrographDir=micrographDir) for key in keys]
        else: # No files
            values = ['NA' for i in range(len(labels))]

        # Set values in metadata
        for label, value in zip(labels, values):
            md.setValue(label, value, objId)
    
    if not md.isEmpty():
        md.sort(xmipp.MDL_MICROGRAPH)
        md.write(summaryFile)
