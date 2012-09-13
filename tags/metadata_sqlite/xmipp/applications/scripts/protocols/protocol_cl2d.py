#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
#
# Example use:
# ./xmipp_protocol_cl2d.py
#
# Author: Carlos Oscar Sanchez Sorzano, July 2009
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {file} Selfile with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
InSelFile='all.sel'
# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='CL2D/classes64'
# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/home/coss/temp/CL2Dsort/Johan'
# {expert} Directory name for logfiles:
""" All logfiles will be stored here
"""
LogDir='Logs'

#------------------------------------------------------------------------------------------------
# {section} Preprocessing parameters
#------------------------------------------------------------------------------------------------
# Sampling rate
""" Sampling rate (Angstroms/Pixel) """
SamplingRate = 1

# Highpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Highpass =0.01

# Lowpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Lowpass =0.4

#------------------------------------------------------------------------------------------------
# {section} Class averages parameters
#------------------------------------------------------------------------------------------------
# Perform CL2D
DoCL2D=False

# Number of references (or classes) to be used:
NumberOfReferences=64

# {expert} Number of initial references
""" Initial number of initial models
"""
NumberOfReferences0=4

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=20

# {expert} Also include mirror transformation in the alignment?
"""  Including the mirror transformation is useful if your particles have a handedness and may fall either face-up or face-down on the grid.
 Note that when you want to use this CL2D run for later RCT reconstruction, you can NOT include the mirror transformation here.
"""
DoMirror=True

# {expert} Use the fast version of this algorithm?
DoFast=True

# {expert}{list}|correntropy|correlation| Comparison method
""" Use correlation or correntropy """
ComparisonMethod='correlation'

# {expert}{list}|classical|robust| Clustering criterion
""" Use the classical clustering criterion or the robust clustering criterion """
ClusteringMethod='classical'

# {expert} Additional parameters for classify_CL2D
""" -verbose, -alignImages, -corr_split, ...
"""
AdditionalParameters='-verbose'

#------------------------------------------------------------------------------------------------
# {section} Core analysis
#------------------------------------------------------------------------------------------------
# Perform core analysis
DoCoreAnalysis=True

# Good class core size (%)
""" A class is a good class if at least this percentage (around 50%) of the
    images assigned to it have been together in all the previous levels.
    Larger values of this parameter tend to keep few good classes. Smaller
    values of this parameter tend to consider more classes as good ones."""
thGoodClass=50

# Junk Zscore
""" Which is the average Z-score to be considered as junk. Typical values
    go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is
    within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores
    accept more images."""
thJunkZscore=2

# PCA Zscore
""" Which is the PCA Z-score to be considered as junk. Typical values
    go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is
    within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores
    accept more images.
    
    This Z-score is measured after projecting onto the PCA space."""
thPCAZscore=2

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of MPI processes to use:
NumberOfMpiProcesses=128

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. 
"""
SystemFlavour='TORQUE-OPENMPI'

# {hidden} This protocol only works in parallel
DoParallel=True

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_cl2d.py'
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

import os,sys,shutil

class CL2D_class:

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 SamplingRate,
                 Highpass,
                 Lowpass,
                 DoCL2D,
                 NumberOfReferences,
                 NumberOfReferences0,
                 NumberOfIterations,
                 DoMirror,
                 DoFast,
                 ComparisonMethod,
                 ClusteringMethod,
                 AdditionalParameters,
                 DoCoreAnalysis,
                 thGoodClass,
                 thJunkZscore,
                 thPCAZscore,
                 NumberOfMpiProcesses,
                 SystemFlavour):
	     
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,selfile

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.SamplingRate=SamplingRate
        self.Highpass=Highpass
        self.Lowpass=Lowpass
        self.DoCL2D=DoCL2D
        self.NumberOfReferences=NumberOfReferences
        self.NumberOfReferences0=NumberOfReferences0
        self.NumberOfIterations=NumberOfIterations
        self.DoMirror=DoMirror
        self.DoFast=DoFast
        self.ComparisonMethod=ComparisonMethod
        self.ClusteringMethod=ClusteringMethod
        self.AdditionalParameters=AdditionalParameters
        self.DoCoreAnalysis=DoCoreAnalysis
        self.thGoodClass=thGoodClass
        self.thJunkZscore=thJunkZscore
        self.thPCAZscore=thPCAZscore
        self.thr=2
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour
   
        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Delete working directory if it exists, make a new one
        if (DoDeleteWorkingDir): 
            if os.path.exists(self.WorkingDir):
                shutil.rmtree(self.WorkingDir)
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Create a selfile with absolute pathname in the WorkingDir
        mysel=selfile.selfile()
        mysel.read(InSelFile)
        newsel=mysel.make_abspath()
        InSelFile=InSelFile.rsplit("/")[-1]
        self.InSelFile=os.path.abspath(self.WorkingDir+'/'+InSelFile)
        newsel.write(self.InSelFile)

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
            os.path.abspath(self.WorkingDir))

        # Execute CL2D in the working directory
        if self.DoCL2D:
            self.execute_CLalign2D()
     
        # Execute CL2D core in the working directory
        if self.DoCoreAnalysis:
            self.execute_core_analysis()
     
        # Finish
        self.close()

    def execute_CLalign2D(self):
        import selfile,launch_job
        if not self.Highpass==0 or not self.Lowpass==0:
            highCutoff=self.Highpass/self.SamplingRate
            lowCutoff=self.Lowpass/self.SamplingRate
            slope=max(self.Highpass/2,0.01)
            preprocessedDir=self.WorkingDir+"/Imgs"
            
            if os.path.exists(preprocessedDir):
                shutil.rmtree(preprocessedDir)
            
            print '*********************************************************************'
            print '*  Copying input images' 
            mysel=selfile.selfile()
            mysel.read(self.InSelFile)
            newsel=mysel.copy_sel(preprocessedDir)
            newsel.write(self.InSelFile)
            params= '-i '+str(self.InSelFile)+\
                    ' -fourier_mask raised_cosine '+str(slope)
            if self.Highpass>0 and self.Lowpass>0:
                params+=" -band_pass "+str(highCutoff)+" "+str(lowCutoff)
            elif self.Highpass>0:
                params+=" -high_pass "+str(highCutoff)
            elif self.Lowpass>0:
                params+=" -low_pass "+str(lowCutoff)
            
            print '*********************************************************************'
            print '*  Executing Fourier filtering program :' 
            launch_job.launch_job("xmipp_fourier_filter",
                                  params,
                                  self.log,
                                  False,
                                  1,
                                  1,
                                  self.SystemFlavour)

        print '*********************************************************************'
        print '*  Executing CL2D program :' 
        if (self.NumberOfReferences0>self.NumberOfReferences):
            self.NumberOfReferences0=self.NumberOfReferences
        params= '-i '+str(self.InSelFile)+' -o '+WorkingDir+'/class '+\
                ' -codes '+str(self.NumberOfReferences)+\
                ' -codes0 '+str(self.NumberOfReferences0)+\
                ' -iter '+str(self.NumberOfIterations)
        params+=' '+self.AdditionalParameters
        if (self.DoFast):
            params+= ' -fast '
        if (not self.DoMirror):
            params+= ' -no_mirror '
        if (self.ComparisonMethod=='correlation'):
            params+= ' -useCorrelation '
        if (self.ClusteringMethod=='classical'):
            params+= ' -classicalMultiref '

        launch_job.launch_job("xmipp_classify_CL2D",
                              params,
                              self.log,
                              True,
                              1,
                              self.thr,
                              self.SystemFlavour)

    def execute_core_analysis(self):
        import launch_job
        print '*********************************************************************'
        print '*  Executing CL2D core analysis :' 
        params= WorkingDir+'/class '+\
                str(self.thGoodClass)+' '+\
                str(self.thJunkZscore)+' '+\
                str(self.thPCAZscore)+' '+\
                str(self.thr)

        launch_job.launch_job("xmipp_classify_CL2D_core_analysis",
                              params,
                              self.log,
                              False,
                              1,
                              1,
                              self.SystemFlavour)
        
    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':
    CL2D=CL2D_class(InSelFile,
                    WorkingDir,
                    DoDeleteWorkingDir,
                    ProjectDir,
                    LogDir,
                    SamplingRate,
                    Highpass,
                    Lowpass,
                    DoCL2D,
                    NumberOfReferences,
                    NumberOfReferences0,
                    NumberOfIterations,
                    DoMirror,
                    DoFast,
                    ComparisonMethod,
                    ClusteringMethod,
                    AdditionalParameters,
                    DoCoreAnalysis,
                    thGoodClass,
                    thJunkZscore,
                    thPCAZscore,
                    NumberOfMpiProcesses,
                    SystemFlavour)
