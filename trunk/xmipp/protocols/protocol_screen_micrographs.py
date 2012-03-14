#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sorzano, October 2011

#from config_protocols import protDict
from protlib_base import *
from protlib_utils import which, runJob, runShowJ
from protlib_filesystem import deleteFile, exists, replaceFilenameExt
import xmipp
from protlib_gui_ext import showError

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
        'ctffind_ctfparam': join('%(micrographDir)s', 'ctfind.ctfparam'),
        'ctffind_spectrum': join('%(micrographDir)s', 'ctffind_spectrum.mrc')
        }

def _getFilename(key, **args):
    return _templateDict[key] % args

class ProtScreenMicrographs(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.screen_micrographs.name, scriptname, project)
        self.Import = "from protocol_screen_micrographs import *"
        importProt = self.getProtocolFromRunName(self.ImportRun) 
        self.importDir = importProt.WorkingDir
        self.TiltPairs = importProt.TiltPairs
        self.importMicroscope = importProt.getFilename('microscope')
        self.importMicrographs = importProt.getFilename('micrographs')
        self.CtffindExec =  which('ctffind3.exe')

    def defineSteps(self):
        self.micrographs = self.getFilename('micrographs')
        # Read Microscope parameters
        MD = xmipp.MetaData(self.importMicroscope)
        objId = MD.firstObject()
        Voltage = MD.getValue(xmipp.MDL_CTF_VOLTAGE,objId)
        SphericalAberration = MD.getValue(xmipp.MDL_CTF_CS,objId)
        AngPix = MD.getValue(xmipp.MDL_CTF_SAMPLING_RATE,objId)
        Magnification = MD.getValue(xmipp.MDL_MAGNIFICATION,objId)

        # Create verifyFiles for the MPI and output directories
        MD = xmipp.MetaData(self.importMicrographs)
        
        # Now the estimation actions
        CtfFindActions = []
        for objId in MD:
            inputFile = MD.getValue(xmipp.MDL_MICROGRAPH, objId)
            micrographName = os.path.basename(inputFile)
            shortname = os.path.splitext(micrographName)[0]
            micrographDir = self.workingDirPath(shortname)                    
            parent_id = self.insertParallelStep('createDir',verifyfiles=[micrographDir],path=micrographDir)

            # Downsample if necessary
            if self.DownsampleFactor != 1:
                finalname = self.tmpPath(shortname + "_tmp.mrc")
                parent_id = self.insertParallelRunJobStep("xmipp_transform_downsample",
                                                   "-i %s -o %s --step %f --method fourier" % (inputFile,finalname,self.DownsampleFactor),
                                                   [finalname],parent_step_id=parent_id)
            else:
                finalname = inputFile
            
            # CTF estimation with Xmipp
            self.insertParallelRunJobStep('xmipp_ctf_estimate_from_micrograph',
                                     "--micrograph "+finalname+\
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
                                     " --defocusU "+str((self.MaxFocus+self.MinFocus)*10000/2),
                                     verifyfiles=[_getFilename('ctfparam', micrographDir=micrographDir)],parent_step_id=parent_id)

            # CTF estimation with Ctffind
            if self.DoCtffind:
                CtfFindActions.append([dict(verifyfiles=[_getFilename('ctffind_ctfparam', micrographDir=micrographDir)],
                                     parent_step_id=parent_id,
                                     CtffindExec=self.CtffindExec,micrograph=finalname,micrographDir=micrographDir,
                                     tmpDir=self.TmpDir,
                                     Voltage=Voltage,SphericalAberration=SphericalAberration,
                                     AngPix=AngPix,Magnification=Magnification,AmplitudeContrast=self.AmplitudeContrast,
                                     LowResolCutoff=self.LowResolCutoff,HighResolCutoff=self.HighResolCutoff,
                                     MinFocus=self.MinFocus,MaxFocus=self.MaxFocus,
                                     StepFocus=self.StepFocus,WinSize=self.WinSize)])
        for action in CtfFindActions:
            self.insertParallelStep('estimateCtfCtffind',action)
        
        # Gather results after external actions
        self.insertStep('gatherResults',verifyfiles=[self.micrographs],
                           TmpDir=self.TmpDir,
                           WorkingDir=self.WorkingDir,
                           summaryFile=self.micrographs,
                           importMicrographs=self.importMicrographs,
                           DoCtffind=self.DoCtffind,
                           Downsampling=self.DownsampleFactor)
    
    def createFilenameTemplates(self):
        return _templateDict
    
    def validate(self):
        errors = []

        # Check that there are any micrograph to process
        if not exists(self.importMicrographs):
            errors.append("Cannot find imported micrographs file:\n   <%s>" % self.importMicrographs)
        else:
            md = xmipp.MetaData(self.importMicrographs)
            if md.isEmpty():
                errors.append("Imported micrographs file <%s> is empty" % self.importMicrographs)
        if not exists(self.importMicroscope):
            errors.append("Cannot find imported microscopy file:\n   <%s>" % self.importMicroscope)
        
        if self.AmplitudeContrast < 0:
            errors.append("Q0 should be positive")
        
        if self.MinFocus < 0 or self.MaxFocus<0:
            errors.append("Defoci range must be positive (minFocus, maxFocus)>0")
            
        if self.MaxFocus < self.MinFocus:
            errors.append("maxFocus must be larger than minFocus")
        
        # Check CTFFIND is available
        if self.DoCtffind and self.CtffindExec=="":
            errors.append("Cannot locate <ctffind3.exe>")
    
        return errors

    def summary(self):
        message = []
        md = xmipp.MetaData(self.importMicrographs)
        message.append("CTF screening of <%d> micrographs from <%s>" % (md.size(), self.importDir))
        if self.DoCtffind:
            message.append("CTF validated with <CTFFIND3>")
        return message
    
    def visualize(self):
        summaryFile = self.getFilename('micrographs')
        if exists(summaryFile):
            runShowJ(summaryFile,extraParams="--mode metadata")
        else:
            summaryFile = summaryFile.replace(self.WorkingDir, self.TmpDir)
            buildSummaryMetadata(self.WorkingDir, self.DoCtffind, self.importMicrographs, summaryFile)
            if exists(summaryFile):
                runShowJ(summaryFile,extraParams="--mode metadata")
            else:
                showError('Error', 'There are not results yet')
    
def estimateCtfCtffind(log,CtffindExec,micrograph,micrographDir,tmpDir,Voltage,SphericalAberration,AngPix,Magnification,
                       AmplitudeContrast,LowResolCutoff,HighResolCutoff,MinFocus,MaxFocus,StepFocus,WinSize):
    # Convert image to MRC
    if not micrograph.endswith('.mrc'):
        import uuid
        deleteTempMicrograph = True
        mrcMicrograph =  join(tmpDir,os.path.splitext(micrograph)[0]+"_"+str(uuid.uuid4())+'.mrc')
        runJob(log,'xmipp_image_convert','-i ' + micrograph + ' -o ' + mrcMicrograph + ' -v 0')
    else:
        deleteTempMicrograph = False
        mrcMicrograph = micrograph;

    # Prepare parameters for CTFTILT
    params = '  << eof > ' + micrographDir + '/ctffind.log\n'
    params += join(micrographDir,mrcMicrograph) + "\n"
    params += micrographDir + '/ctffind_spectrum.mrc\n'
    params += str(SphericalAberration) + ',' + \
              str(Voltage) + ',' + \
              str(AmplitudeContrast) + ',' + \
              str(Magnification) + ',' + \
              str(Magnification*AngPix*1e-4) + "\n"
    params += str(WinSize) + ',' + \
              str(AngPix / LowResolCutoff) + ',' + \
              str(AngPix / HighResolCutoff) + ',' + \
              str(MinFocus*10000) + ',' + \
              str(MaxFocus*10000) + ',' + \
              str(StepFocus*10000) + "\n"
    runJob(log, "export NATIVEMTZ=kk ; "+CtffindExec,params)

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
    MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE, float(AngPix), objId)
    MD.setValue(xmipp.MDL_CTF_VOLTAGE,       float(Voltage), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUSU,      float(DF2), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUSV,      float(DF1), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, float(Angle), objId)
    MD.setValue(xmipp.MDL_CTF_CS,            float(SphericalAberration), objId)
    MD.setValue(xmipp.MDL_CTF_Q0,            float(-AmplitudeContrast), objId)
    MD.setValue(xmipp.MDL_CTF_K,             1.0, objId)
    MD.write(fnOut)

def gatherResults(log,TmpDir,WorkingDir,summaryFile, importMicrographs,DoCtffind,Downsampling):
    buildSummaryMetadata(WorkingDir, DoCtffind, importMicrographs, summaryFile)
    runJob(log,"xmipp_ctf_sort_psds","-i " + summaryFile)
    if Downsampling!=1:
        runJob(log,"rm","-f "+TmpDir+"/*")

def buildSummaryMetadata(WorkingDir,DoCtffind,importMicrographs,summaryFile):
    md = xmipp.MetaData()
    importMd = xmipp.MetaData(importMicrographs)
    for id in importMd:
        inputFile = importMd.getValue(xmipp.MDL_MICROGRAPH,id)
        micrographName = os.path.basename(inputFile)
        shortname = replaceFilenameExt(micrographName, '')
        micrographDir = join(WorkingDir, shortname)                    
        
        objId = md.addObject()
        md.setValue(xmipp.MDL_MICROGRAPH, inputFile, objId)
        ctfparam = _getFilename('ctfparam', micrographDir=micrographDir)
        labels = [xmipp.MDL_PSD, xmipp.MDL_PSD_ENHANCED, xmipp.MDL_CTFMODEL,xmipp.MDL_ASSOCIATED_IMAGE1, xmipp.MDL_ASSOCIATED_IMAGE2]
        if exists(ctfparam): # Get filenames
            keys = ['psd', 'enhanced_psd', 'ctfparam', 'ctfmodel_quadrant', 'ctfmodel_halfplane']
            values = [_getFilename(key, micrographDir=micrographDir) for key in keys]
        else: # No files
            values = ['NA' for i in range(len(labels))]

        if DoCtffind:
            ctffindCTF = _getFilename('ctffind_ctfparam')
            labels += [xmipp.MDL_CTFMODEL2, xmipp.MDL_ASSOCIATED_IMAGE3]
            if exists(ctffindCTF):
                values += [ctffindCTF, _getFilename('ctffind_spectrum')]
            else:
                values += ['NA', 'NA']
        # Set values in metadata
        for label, value in zip(labels, values):
            md.setValue(label, value, objId)
            
    if not md.isEmpty():
        md.sort(xmipp.MDL_MICROGRAPH)
        md.write(summaryFile)