#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, October 2011

#from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob, runShowJ,printLog
from protlib_filesystem import deleteFile, createLink2, exists, replaceFilenameExt, createDir
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
        'ctfmodel_halfplane': _prefix + '_ctfmodel_halfplane.xmp',
        'ctffind_ctfparam': join('%(micrographDir)s', 'ctffind.ctfparam'),
        'ctffind_spectrum': join('%(micrographDir)s', 'ctffind_spectrum.mrc')
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
        self.MicrographsMd = self.Input['micrographs']
        if self.TiltPairs:
            self.MicrographsMd='micrographPairs@'+self.MicrographsMd

    def defineSteps(self):
        extraDir=self.workingDirPath('extra')
        parent_id = self.insertStep('createDir',verifyfiles=[extraDir],path=extraDir)

        filesToImport = [self.Input[k] for k in ['microscope', 'acquisition']]
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
        MD = xmipp.MetaData(self.MicrographsMd)
        #if removed in import do not process them
        MD.removeDisabled()
	
        # Now the estimation actions
        for objId in MD:
            inputFile = MD.getValue(xmipp.MDL_MICROGRAPH, objId)
            micrographName = os.path.basename(inputFile)
            shortname = os.path.splitext(micrographName)[0]
            micrographDir = os.path.join(extraDir,shortname)                    

            self.insertParallelStep('estimateSingleCTF',verifyfiles=[_getFilename('ctfparam', micrographDir=micrographDir)],
                                    WorkingDir=self.WorkingDir, inputFile=inputFile, DownsampleFactor=self.DownsampleFactor,
                                    AutomaticDownsampling=self.AutomaticDownsampling, Voltage=Voltage,
                                    SphericalAberration=SphericalAberration, AngPix=AngPix, AmplitudeContrast=self.AmplitudeContrast,
                                    LowResolCutoff=self.LowResolCutoff, HighResolCutoff=self.HighResolCutoff,WinSize=self.WinSize,
                                    MaxFocus=self.MaxFocus,MinFocus=self.MinFocus,FastDefocus=self.FastDefocus,
                                    parent_step_id=XmippProjectDb.FIRST_STEP)
        
        # Gather results after external actions
        self.insertStep('gatherResults',verifyfiles=[self.micrographs],
                           TmpDir=self.TmpDir,
                           WorkingDir=self.WorkingDir,
                           summaryFile=self.micrographs,
                           importMicrographs=self.MicrographsMd,
                           Downsampling=self.DownsampleFactor,
                           NumberOfMpi=self.NumberOfMpi)
        if self.AutomaticRejection!="":
            self.insertStep("automaticRejection",WorkingDir=self.WorkingDir,condition=self.AutomaticRejection)
        if self.TiltPairs:
            self.insertStep("copyTiltPairs",WorkingDir=self.WorkingDir,inputMicrographs=self.MicrographsMd)
    
    def createFilenameTemplates(self):
        return _templateDict
    
    def validate(self):
        errors = []

        if self.DownsampleFactor<1:
            errors.append("Downsampling must be >=1");

        micsFn = self.Input['micrographs']
        # Check that there are any micrograph to process
        if xmippExists(micsFn):
            md = xmipp.MetaData(micsFn)
            if md.isEmpty():
                errors.append("Imported micrographs file <%(micrographs)s> is empty" % self.Input)
#        else:
#            errors.append("Expected micrographs file <%(micrographs)s> is missing" % self.Input)
            
        if self.AmplitudeContrast < 0:
            errors.append("Q0 should be positive")
        
        if self.MinFocus < 0 or self.MaxFocus<0:
            errors.append("Defoci range must be positive (minFocus, maxFocus)>0")
            
        if self.MaxFocus < self.MinFocus:
            errors.append("maxFocus must be larger than minFocus")

        return errors

    def summary(self):
        message = []
        from protlib_xmipp import getMdSize
        size = getMdSize(self.Input['micrographs'])
        message.append("CTF screening of <%d> micrographs." % size)
        message.append("Input directory: [%s]" % self.PrevRun.WorkingDir)
        return message
    
    def papers(self):
        papers=[]
        papers.append('Jonic, JSB (2007) [http://www.ncbi.nlm.nih.gov/pubmed/16987671]')
        papers.append('Sorzano, JSB (2007) [http://www.ncbi.nlm.nih.gov/pubmed/17911028]')
        papers.append('Sorzano, BMC SB (2009) [http://www.ncbi.nlm.nih.gov/pubmed/19321015]')
        papers.append('Sorzano, Meth.Mol.Biol. (2013) [http://www.ncbi.nlm.nih.gov/pubmed/23086876]')
        if self.FastDefocus:
            papers.append('Vargas, JSB (2013) [http://www.ncbi.nlm.nih.gov/pubmed/23261401]')
        return papers
    
    def visualize(self):
        summaryFile = self.getFilename('micrographs')
        
        if not exists(summaryFile): # Try to create partial summary file
            summaryFile = summaryFile.replace(self.WorkingDir, self.TmpDir)
            buildSummaryMetadata(self.WorkingDir, self.Input['micrographs'], summaryFile)
        
        if exists(summaryFile):
            self.regenerateSummary(summaryFile)
            runShowJ(summaryFile, extraParams = "--mode metadata")
        else:
            showWarning('Warning', 'There are not results yet',self.master)
    
    def regenerateSummary(self,summaryFile):
        import time
        summaryTime=time.ctime(os.path.getmtime(summaryFile))
        md=xmipp.MetaData(summaryFile)
        regenerate=False
        for objId in md:
            fnCTF=md.getValue(xmipp.MDL_CTF_MODEL,objId)
            if fnCTF!="NA":
                ctfTime=time.ctime(os.path.getmtime(fnCTF))
                if ctfTime>summaryTime:
                    regenerate=True
                    break
        if regenerate:
            print("Regenerating "+summaryFile+" because there are newer CTFs")
            gatherResults(self.Log,TmpDir=self.TmpDir,
                          WorkingDir=self.WorkingDir,
                          summaryFile=self.micrographs,
                          importMicrographs=self.MicrographsMd,
                          Downsampling=self.DownsampleFactor,
                          NumberOfMpi=self.NumberOfMpi)

def estimateSingleCTF(log, WorkingDir, inputFile, DownsampleFactor, AutomaticDownsampling,
                      Voltage, SphericalAberration, AngPix, AmplitudeContrast, LowResolCutoff, 
                      HighResolCutoff,WinSize,MaxFocus,MinFocus,FastDefocus):
    extraDir=os.path.join(WorkingDir,'extra')
    tmpDir=os.path.join(WorkingDir,'tmp')
    micrographName = os.path.basename(inputFile)
    shortname = os.path.splitext(micrographName)[0]
    micrographDir = os.path.join(extraDir,shortname)
    oroot="%s/xmipp_ctf"%micrographDir           

    if not os.path.exists(micrographDir):
        createDir(log,micrographDir)

    downsampleList=[DownsampleFactor]
    if AutomaticDownsampling:
        downsampleList.append(DownsampleFactor+1)
        if DownsampleFactor>=2:
            downsampleList.append(DownsampleFactor-1)
        else:
            downsampleList.append(DownsampleFactor/2)
    
    for DownsampleFactor in downsampleList:
        # Downsample if necessary
        deleteTmp=False
        if DownsampleFactor != 1:
            finalname = os.path.join(tmpDir,shortname + "_tmp.mrc")
            runJob(log,"xmipp_transform_downsample","-i %s -o %s --step %f --method fourier" % (inputFile,finalname,DownsampleFactor))
            deleteTmp=True
        else:
            finalname = inputFile
            
        # CTF estimation with Xmipp
        args="--micrograph "+finalname+\
             " --oroot "+oroot+\
             " --kV "+str(Voltage)+\
             " --Cs "+str(SphericalAberration)+\
             " --sampling_rate "+str(AngPix*DownsampleFactor)+\
             " --downSamplingPerformed "+str(DownsampleFactor)+\
             " --ctfmodelSize 256"+\
             " --Q0 "+str(AmplitudeContrast)+\
             " --min_freq "+str(LowResolCutoff)+\
             " --max_freq "+str(HighResolCutoff)+\
             " --pieceDim "+str(WinSize)+\
             " --defocus_range "+str((MaxFocus-MinFocus)*10000/2)+\
             " --defocusU "+str((MaxFocus+MinFocus)*10000/2)
        if (FastDefocus):
            args+=" --fastDefocus"
        runJob(log,'xmipp_ctf_estimate_from_micrograph', args)
        
        if deleteTmp:
            deleteFile(log, finalname)
        
        md = xmipp.MetaData()
        id = md.addObject()
        md.setValue(xmipp.MDL_MICROGRAPH,inputFile,id)
        md.setValue(xmipp.MDL_CTF_MODEL,oroot+".ctfparam",id)
        md.setValue(xmipp.MDL_PSD,oroot+".psd",id)
        fnEval=os.path.join(tmpDir,shortname+".xmd")
        md.write(fnEval)
        criterion="ctfCritFirstZero<5 OR ctfCritMaxFreq>20 OR ctfCritfirstZeroRatio<0.9 OR ctfCritfirstZeroRatio>1.1 OR "\
                  "ctfCritFirstMinFirstZeroRatio>10 OR ctfCritCorr13<0 OR ctfCritCtfMargin<0 OR ctfCritNonAstigmaticValidty<0.3 OR " \
                  "ctfCritNonAstigmaticValidty>25"
        runJob(log,"xmipp_ctf_sort_psds","-i %s --downsampling %f"%(fnEval,DownsampleFactor))
        fnRejected=os.path.join(tmpDir,shortname+"_rejected.xmd")
        runJob(log,"xmipp_metadata_utilities","-i %s --query select %s -o %s"%(fnEval,criterion,fnRejected))
        md.read(fnRejected)
        if md.size()==0:
            break

def gatherResults(log,TmpDir,WorkingDir,summaryFile, importMicrographs,Downsampling,NumberOfMpi):
    buildSummaryMetadata(WorkingDir, importMicrographs, summaryFile)
    dirSummary,fnSummary=os.path.split(summaryFile)
    runJob(log,"xmipp_ctf_sort_psds","-i %s -o %s/aux_%s --downsampling %f"%(summaryFile,dirSummary,fnSummary,Downsampling),
           NumberOfMpi=NumberOfMpi)
    runJob(log,"mv","-f %s/aux_%s %s"%(dirSummary,fnSummary,summaryFile))
    runJob(log,"touch",summaryFile)
    if Downsampling!=1:
        runJob(log,"find",TmpDir+" -type f -exec rm {} \;")

def buildSummaryMetadata(WorkingDir,importMicrographs,summaryFile):
    md = xmipp.MetaData()
    importMd = xmipp.MetaData(importMicrographs)
    importMd.removeDisabled()
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
            values = ['NA'] * len(labels)

        # Set values in metadata
        for label, value in zip(labels, values):
            md.setValue(label, value, objId)
    
    if not md.isEmpty():
        md.sort(xmipp.MDL_MICROGRAPH)
        md.write(summaryFile)

def automaticRejection(log,WorkingDir,condition):
    fnMic=os.path.join(WorkingDir,"micrographs.xmd")
    fnRejected=os.path.join(WorkingDir,"tmp/rejectedMicrographs.xmd")
    runJob(log,"xmipp_metadata_utilities",'-i %s --query select "%s" -o %s'%(fnMic,condition,fnRejected))
    md = xmipp.MetaData(fnRejected)
    fnRejectedMicrographs=md.getColumnValues(xmipp.MDL_MICROGRAPH)
    md.read(fnMic)
    for id in md:
        fnCurrentMicrograph=md.getValue(xmipp.MDL_MICROGRAPH,id)
        if fnCurrentMicrograph in fnRejectedMicrographs:
            md.setValue(xmipp.MDL_ENABLED,-1,id)
    md.write(fnMic)
    deleteFile(log,fnRejected)

def copyTiltPairs(log,WorkingDir,inputMicrographs):
    md=xmipp.MetaData(inputMicrographs)
    md.write("micrographPairs@"+os.path.join(WorkingDir,"micrographs.xmd"),xmipp.MD_APPEND)

