#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs: 
#  - downsampling
#  - power spectral density (PSD) and CTF estimation on the micrograph
#
# For each micrograph given, this script will perform 
# the requested operations below.
# For each micrograph a subdirectory will be created
#
# Author: Sjors Scheres, March 2007
#         Roberto Marabini (mpi extension)
#         Carlos Oscar Sorzano, November 2010


import glob, os, shutil, sys, time
import xmipp
from protlib_utils import which, runJob

def preprocessMicrograph(log,micrograph,micrographDir,DoPreprocess,Crop,Stddev,Down):
    # Decide name after preprocessing
    finalname=os.path.join(micrographDir,'micrograph')
    (filepath, micrographName)=os.path.split(relpath(micrograph))
    (shortname2, extension)=os.path.splitext(micrographName)        
    if DoPreprocess and (not Stddev == -1 or not Crop == -1 or not Down == 1):
        finalname += ".mrc"
    else:
        finalname += extension
        if not os.path.exists(finalname):
            runJob(log,"ln",'-s ' + relpath(micrograph, micrographDir) + ' ' + finalname)
            if micrograph.endswith(".raw"):
                runJob(log,"ln",'-s ' + relpath(micrograph, micrographDir)+'.inf' + ' ' + finalname+'.inf')
            return
    if not self.DoPreprocess:
        return
    
    # Crop
    iname=micrograph
    if not self.Crop == -1:
        runJob(log,"xmipp_transform_window"," -i " + iname + " -o " + finalname + " --crop " + str(Crop) + " -v 0")
        iname=finalname
    
    # Remove bad pixels
    if not Stddev == -1:
        params = " -i " + iname + " --bad_pixels outliers " + str(Stddev)+" -v 0"
        if not iname == finalname:
            params += " -o " + finalname
            iname=finalname
        runJob(log,"xmipp_transform_filter",params)
    
    # Downsample
    if not Down == 1:
        tmpFile=os.path.join(micrographDir,"tmp.mrc")
        runJob(log,"xmipp_transform_downsample","-i " + iname + " -o " + tmpFile + " --step "+ str(Down)+' --method fourier')
        runJob("mv -fi "+tmpFile+" "+finalname)

def estimateCtfXmipp(log,micrograph,micrographDir,Voltage,SphericalAberration,AngPix,
                     AmplitudeContrast,LowResolCutoff,HighResolCutoff,MinFocus,MaxFocus,WinSize):
    params="--micrograph "+os.path.join(micrographDir,micrograph)+" --oroot "+os.path.join(micrographDir,"xmipp_ctf")+\
           " --kV "+str(Voltage)+\
           " --Cs "+str(SphericalAberration)+\
           " --sampling_rate "+str(AngPix)+\
           " --ctfmodelSize 256"+\
           " --Q0 "+str(AmplitudeContrast)+\
           " --min_freq "+str(LowResolCutoff)+\
           " --max_freq "+str(HighResolCutoff)+\
           " --pieceDim "+str(WinSize)+\
           " --defocusRange "+str((MaxFocus-MinFocus)/2)+\
           " --defocusU "+str((MaxFocus+MinFocus)/2)
    runJob(log,"xmipp_ctf_estimate_from_micrograph",params)

def estimateCtfCtffind(log,CtffindExec,micrograph,micrographDir,Voltage,SphericalAberration,AngPix,Magnification,
                       AmplitudeContrast,LowResolCutoff,HighResolCutoff,MinFocus,MaxFocus,StepFocus,WinSize):
    # Convert image to MRC
    if not micrograph.endswith('.mrc'):
        mrcMicrograph= os.path.join(micrographDir,'tmp.mrc')
        runJob(log,'xmipp_image_convert','-i ' + micrograph + ' -o ' + mrcMicrograph + ' -v 0')
    else:
        mrcMicrograph = micrograph;

    # Prepare parameters for CTFTILT
    (filepath, micrographName)=os.path.split(micrograph)
    (fnRoot, extension)=os.path.splitext(micrographName)
    params += '  << eof > ' + micrographDir + '/ctffind.log\n'
    params += mrcMicrograph + "\n"
    params += micrographDir + '/ctffind_spectrum.mrc\n'
    params += str(SphericalAberration) + ',' + \
              str(Voltage) + ',' + \
              str(AmplitudeContrast) + ',' + \
              str(Magnification) + ',' + \
              str(Magnification*AngPix*1e-4) + "\n"
    params += str(WinSizeCTF) + ',' + \
              str(AngPix / LowResolCutoff) + ',' + \
              str(AngPix / HighResolCutoff) + ',' + \
              str(MinFocus) + ',' + \
              str(MaxFocus) + ',' + \
              str(StepFocus) + "\n"
    runJob(log, "export NATIVEMTZ=kk ; "+self.CtffindExec)

    fnOut=micrographDir + '/ctffind.ctfparam'

    # Remove temporary files
    if os.path.exists(micrographDir + '/tmp.mrc'):
        os.remove(micrographDir + '/tmp.mrc')

    # Pick values from ctffind
    if not os.path.exists(micrographDir + '/ctffind.log'):
        return
    
    # Effectively pickup results
    fh=open(micrographDir + '/ctffind.log', 'r')
    lines=fh.readlines()
    fh.close()
    DF1=0.
    DF2=0.
    Angle=0.
    found=False
    for i in range(len(lines)):
        if not (lines[i].find('Final Values') == -1):
            words=lines[i].split()
            DF1=float(words[0])
            DF2=float(words[1])
            Angle=float(words[2])
            found=True
            break

    if not found:
        return
    
    # Generate Xmipp .ctfparam file:
    MD=xmipp.MetaData()
    MD.setColumnFormat(False)
    objId=MD.addObject()
    MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE, AngPix, objId)
    MD.setValue(xmipp.MDL_CTF_VOLTAGE,       float(self.Voltage), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUSU,      float(-DF2), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUSV,      float(-DF1), objId)
    MD.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, float(Angle), objId)
    MD.setValue(xmipp.MDL_CTF_CS,            float(SphericalAberration), objId)
    MD.setValue(xmipp.MDL_CTF_Q0,            float(-AmplitudeContrast), objId)
    MD.setValue(xmipp.MDL_CTF_K,             1.0, objId)
    MD.write(fnOut)

#        
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtPreprocessMicrographs)