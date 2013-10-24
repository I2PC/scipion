#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, October 2011

#from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob, runShowJ,printLog
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
        if not self.TiltPairs:
            self.MicrographsMd = self.Input['micrographs']
        else:
            self.inputFilename('tilted_pairs')
            self.MicrographsMd = self.Input['tilted_pairs']

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
        MD = xmipp.MetaData(self.MicrographsMd)
        #if removed in import do not process them
        MD.removeDisabled()
	
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
                           importMicrographs=self.MicrographsMd,
                           Downsampling=self.DownsampleFactor,
                           NumberOfMpi=self.NumberOfMpi)
        if self.AutomaticRejection!="":
            self.insertStep("automaticRejection",WorkingDir=self.WorkingDir,condition=self.AutomaticRejection)
    
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
        summaryTime=os.path.getmtime(summaryFile)
        md=xmipp.MetaData(summaryFile)
        regenerate=False
        for objId in md:
            fnCTF=md.getValue(xmipp.MDL_CTF_MODEL,objId)
            ctfTime=os.path.getmtime(fnCTF)
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
    
def estimateCtfCtffind1(_log, micrograph,
                          micrographDir,
                          oroot,
                          kV,
                          Cs,
                          sampling_rate,
                          downSamplingPerformed,
                          ctfmodelSize,
                          Q0,
                          min_freq,
                          max_freq,
                          pieceDim,
                          MinFocus,
                          MaxFocus,
                          StepFocus
                          ):
        #def estimateCtfCtffind(log,CtffindExec,micrograph,micrographDir,tmpDir,Voltage,SphericalAberration,AngPix,Magnification,
        #                       DownsampleFactor,AmplitudeContrast,LowResolCutoff,HighResolCutoff,MinFocus,MaxFocus,StepFocus,WinSize):
        # Convert image to MRC
        printLog(_log)
        if not micrograph.endswith('.mrc'):
            from protlib_filesystem import uniqueRandomFilename
            deleteTempMicrograph = True
            fnMicrograph=os.path.split(micrograph)[1]
            mrcMicrograph =  join(tmpDir,uniqueRandomFilename(os.path.splitext(fnMicrograph)[0]+'.mrc'))
            runJob(log,'xmipp_image_convert','-i ' + micrograph + ' -o ' + mrcMicrograph + ' -v 0')
        else:
            deleteTempMicrograph = False
            mrcMicrograph = micrograph;

        #ctffind3.exe << eof
        #$1.mrc
        #test.mrc
        #2.26,200.0,0.10,60000.0,14.4
        #256,96,0.8,5000.0,100000,1000.0
        #eof

        # Prepare parameters for CTFTILT
        # multiply Q0 by -1?
        
        #since we only have sampling rate set magnification to a suitable constant
        Magnification=60000
        params = '  << eof > ' + micrographDir + '/ctffind.log\n'
        params += mrcMicrograph + "\n"
        params += micrographDir + '/ctffind_spectrum.mrc\n'
        params += str(Cs) + ',' + \
                  str(kV) + ',' + \
                  str(Q0 ) + ',' + \
                  str(Magnification/downSamplingPerformed) + ',' + \
                  str(Magnification/downSamplingPerformed*sampling_rate*1e-4) + "\n"
        params += str(pieceDim) + ',' + \
                  str(sampling_rate*downSamplingPerformed / min_freq) + ',' + \
                  str(sampling_rate*downSamplingPerformed / max_freq) + ',' + \
                  str(MinFocus*10000) + ',' + \
                  str(MaxFocus*10000) + ',' + \
                  str(StepFocus*10000) + "\n"

        CtffindExec =  which('ctffind3.exe')
        runJob(_log, "export NATIVEMTZ=kk ; "+CtffindExec,params)
    
        fnOut = _getFilename('ctffind_ctfparam', micrographDir=micrographDir)
    
        # Remove temporary files
        if deleteTempMicrograph:
            deleteFile(log, mrcMicrograph)
    
        # Pick values from ctffind
        ctffindLog = join(micrographDir,'ctffind.log')
        if not exists(ctffindLog):
            raise xmipp.XmippError("Cannot find "+ctffindLog)
        
        # Effectively pickup results
        fh = open(ctffindLog, 'r')
        lines = fh.readlines()
        fh.close()
        DF1 = 0.
        DF2 = 0.
        Angle = 0.
        found = False
        for i in range(len(lines)):
            if not (lines[i].find('Final Values') == -1):
                words = lines[i].split()
                DF1 = float(words[0])
                DF2 = float(words[1])
                Angle = float(words[2])
                found = True
                break
        if not found:
            raise xmipp.XmippError("Cannot find defocus values in "+ctffindLog)
        
        # Generate Xmipp .ctfparam file:
        MD = xmipp.MetaData()
        MD.setColumnFormat(False)
        objId = MD.addObject()
        MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE, float(sampling_rate), objId)
        MD.setValue(xmipp.MDL_CTF_VOLTAGE,       float(kV), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUSU,      float(DF2), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUSV,      float(DF1), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, float(Angle), objId)
        MD.setValue(xmipp.MDL_CTF_CS,            float(Cs), objId)
        MD.setValue(xmipp.MDL_CTF_Q0,            float(Q0), objId)
        MD.setValue(xmipp.MDL_CTF_K,             1.0, objId)
        MD.write(fnOut)

def gatherResults(log,TmpDir,WorkingDir,summaryFile, importMicrographs,Downsampling,NumberOfMpi):
    buildSummaryMetadata(WorkingDir, importMicrographs, summaryFile)
    dirSummary,fnSummary=os.path.split(summaryFile)
    runJob(log,"xmipp_ctf_sort_psds","-i %s -o %s/aux_%s"%(summaryFile,dirSummary,fnSummary),NumberOfMpi=NumberOfMpi)
    runJob(log,"mv","-f %s/aux_%s %s"%(dirSummary,fnSummary,summaryFile))
    runJob(log,"touch",summaryFile)
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
