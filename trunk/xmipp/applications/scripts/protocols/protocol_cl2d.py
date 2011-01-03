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
InSelFile='sort_junk.sel'
# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='CL2D/classes64'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/gpfs/fs1/home/bioinfo/coss/analu/Mcm467_cdt1_sinOli-7AporPix'
# {expert} Directory name for logfiles:
""" All logfiles will be stored here
"""
LogDir='Logs'

#------------------------------------------------------------------------------------------------
# {section} Class averages parameters
#------------------------------------------------------------------------------------------------
# Number of references (or classes) to be used:
NumberOfReferences=64

# {expert} Number of initial references
""" Initial number of initial models
"""
NumberOfReferences0=4

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=15

# {expert} Band pass filter
""" Apply a band pass filter before clustering """
DoFilter =True

# {expert} Highpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Highpass =0.02

# {expert} Lowpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Lowpass =0.4

# {expert}{list}|correntropy|correlation| Comparison method
""" Use correlation or correntropy """
ComparisonMethod='correlation'

# {expert}{list}|classical|robust| Clustering criterion
""" Use the classical clustering criterion or the robust clustering criterion """
ClusteringMethod='classical'

# {expert} Additional parameters for classify_CL2D
""" -verbose, -corr_split, ...
"""
AdditionalParameters=''

#------------------------------------------------------------------------------------------------
# {section} Core analysis
#------------------------------------------------------------------------------------------------
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
thJunkZscore=3

# PCA Zscore
""" Which is the PCA Z-score to be considered as junk. Typical values
    go from 1.5 to 3. For the Gaussian distribution 99.5% of the data is
    within a Z-score of 3. Lower Z-scores reject more images. Higher Z-scores
    accept more images.
    
    This Z-score is measured after projecting onto the PCA space."""
thPCAZscore=3

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of MPI processes to use:
NumberOfMpiProcesses=40

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

def getParameter(prm,filename):
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    for line in lines:
        tokens=line.split('=')
        if tokens[0]==prm:
            return tokens[1].strip()
    return ""

class CL2D_class:
    def saveAndCompareParameters(self, listOfParameters):
        fnOut=self.WorkingDir + "/protocolParameters.txt"
        linesNew=[];
        for prm in listOfParameters:
            eval("linesNew.append('"+prm +"='+str("+prm+")+'\\n')")
        if os.path.exists(fnOut):
            f = open(fnOut, 'r')
            linesOld=f.readlines()
            f.close()
            same=True;
            if len(linesOld)==len(linesNew):
                for i in range(len(linesNew)):
                    if not linesNew[i]==linesOld[i]:
                        same=False
                        break;
            else:
                same=False
            if not same:
                print("Deleting")
                self.log.info("Deleting working directory since it is run with different parameters")
                shutil.rmtree(self.WorkingDir)
                os.makedirs(self.WorkingDir)
        f = open(fnOut, 'w')
        f.writelines(linesNew)
        f.close()

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 ProjectDir,
                 LogDir,
                 NumberOfReferences,
                 NumberOfReferences0,
                 NumberOfIterations,
                 DoFilter,
                 Highpass,
                 Lowpass,
                 ComparisonMethod,
                 ClusteringMethod,
                 AdditionalParameters,
                 thGoodClass,
                 thJunkZscore,
                 thPCAZscore,
                 NumberOfMpiProcesses,
                 SystemFlavour):
	     
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.InSelFile=InSelFile
        self.NumberOfReferences=NumberOfReferences
        self.NumberOfReferences0=NumberOfReferences0
        self.NumberOfIterations=NumberOfIterations
        self.DoFilter=DoFilter
        self.Highpass=Highpass
        self.Lowpass=Lowpass
        self.ComparisonMethod=ComparisonMethod
        self.ClusteringMethod=ClusteringMethod
        self.AdditionalParameters=AdditionalParameters
        self.thGoodClass=thGoodClass
        self.thJunkZscore=thJunkZscore
        self.thPCAZscore=thPCAZscore
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour
   
        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Create directory if does not exist
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)
            self.doStep1=True
            self.doStep2=True
            self.doStep3=True
        else:
            fnParam=self.WorkingDir + "/protocolParameters.txt"
            oldInSelFile=getParameter("InSelFile",fnParam)
            oldDoFilter=getParameter("DoFilter",fnParam)
            oldHighpass=getParameter("Highpass",fnParam)
            oldLowpass=getParameter("Lowpass",fnParam)

            oldNumberOfReferences=getParameter("NumberOfReferences",fnParam)
            oldNumberOfReferences0=getParameter("NumberOfReferences0",fnParam)
            oldNumberOfIterations=getParameter("NumberOfIterations",fnParam)
            oldComparisonMethod=getParameter("ComparisonMethod",fnParam)
            oldClusteringMethod=getParameter("ClusteringMethod",fnParam)
            oldAdditionalParameters=getParameter("AdditionalParameters",fnParam)
            
            oldthGoodClass=getParameter("thGoodClass",fnParam)
            oldthJunkZscore=getParameter("thJunkZscore",fnParam)
            oldthPCAZscore=getParameter("thPCAZscore",fnParam)
            self.doStep1=False
            self.doStep2=False
            self.doStep3=False
            if oldInSelFile!=InSelFile or DoFilter!=oldDoFilter or \
               oldHighpass!=Highpass or oldLowpass!=Lowpass:
                self.doStep1=True
                self.doStep2=True
                self.doStep3=True
            elif oldNumberOfReferences!=NumberOfReferences or \
               oldNumberOfReferences0!=NumberOfReferences0 or \
               oldNumberOfIterations!=NumberOfIterations or \
               oldComparisonMethod!=ComparisonMethod or \
               oldClusteringMethod!=ClusteringMethod or \
               oldAdditionalParameters!=AdditionalParameters:
                self.doStep2=True
                self.doStep3=True
            elif oldthGoodClass!=thGoodClass or oldthJunkZscore!=thJunkZscore or \
               oldthPCAZscore!=thPCAZscore:
                self.doStep3=True                
                
        # Save parameters and compare to possible previous runs
        self.saveAndCompareParameters([
                 "InSelFile",
                 "NumberOfReferences",
                 "NumberOfReferences0",
                 "NumberOfIterations",
                 "DoFilter",
                 "Highpass",
                 "Lowpass",
                 "ComparisonMethod",
                 "ClusteringMethod",
                 "AdditionalParameters",
                 "thGoodClass",
                 "thJunkZscore",
                 "thPCAZscore"]);

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
            os.path.abspath(self.WorkingDir))

        # Run
        self.preprocess()
        self.execute_CLalign2D()
        self.execute_core_analysis()

    def preprocess(self):
        import launch_job
        if self.DoFilter:
            slope=max(self.Highpass/2,0.01)
            fnOut=self.WorkingDir+'/preprocessedImages.stk'
            self.selFileToUse=fnOut
            if self.doStep1:
                params= '-i '+str(self.InSelFile)+\
                        ' --fourier_mask raised_cosine '+str(slope)+\
                        ' -o '+fnOut
                if self.Highpass>0 and self.Lowpass>0:
                    params+=" --band_pass "+str(self.Highpass)+" "+str(self.Lowpass)
                elif self.Highpass>0:
                    params+=" --high_pass "+str(self.Highpass)
                elif self.Lowpass>0:
                    params+=" --low_pass "+str(self.Lowpass)
                launch_job.launch_job("xmipp_fourier_filter",
                                      params,
                                      self.log,
                                      False,
                                      1,
                                      1,
                                      self.SystemFlavour)
        else:
            self.selFileToUse=self.InSelFile

    def execute_CLalign2D(self):
        import launch_job
        if not self.doStep2:
            return
        params= '-i '+str(self.selFileToUse)+' -o '+WorkingDir+'/class '+\
                ' -codes '+str(self.NumberOfReferences)+\
                ' -codes0 '+str(self.NumberOfReferences0)+\
                ' -iter '+str(self.NumberOfIterations)
        params+=' '+self.AdditionalParameters
        if (self.ComparisonMethod=='correlation'):
            params+= ' -useCorrelation '
        if (self.ClusteringMethod=='classical'):
            params+= ' -classicalMultiref '

        launch_job.launch_job("xmipp_classify_CL2D",
                              params,
                              self.log,
                              True,
                              self.NumberOfMpiProcesses,
                              1,
                              self.SystemFlavour)
        
        # If the images have been filtered, translate the output of cl2d in
        # terms of the original images
        # Construct dictionary
        fnClasses=WorkingDir+'/class.sel'
        if os.path.exists(fnClasses) and self.DoFilter:
            import xmipp
            xmipp.substituteOriginalImages(fnClasses,self.InSelFile,
                                           self.WorkingDir+"/intermediate.sel",
                                           xmipp.MDL_IMAGE_ORIGINAL)
            os.system("mv -f "+self.WorkingDir+"/intermediate.sel "+fnClasses)

    def execute_core_analysis(self):
        import launch_job
        if not self.doStep2:
            return
        params= WorkingDir+'/class '+\
                str(self.thGoodClass)+' '+\
                str(self.thJunkZscore)+' '+\
                str(self.thPCAZscore)+' '+\
                str(self.NumberOfMpiProcesses)

        launch_job.launch_job("xmipp_classify_CL2D_core_analysis",
                              params,
                              self.log,
                              True,
                              self.NumberOfMpiProcesses,
                              1,
                              self.SystemFlavour)
        
# Preconditions
def preconditions(gui):
    retval=True
    # Check if there is workingdir
    if WorkingDir == "":
        message="No working directory given"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    # Check that there are any micrograph to process
    if not os.path.exists(InSelFile):
        message="The input selfile is not valid"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    # Check that the number of classes is correct
    if NumberOfReferences0<=0:
        message="The number of initial classes must be positive"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
        
    # Check that the number of classes is correct
    if NumberOfReferences0>NumberOfReferences:
        message="The number of initial classes cannot be larger than the number of final classes"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    # Check filter parameters
    if DoFilter and (Highpass<0 or Highpass>0.5 or Lowpass<0 or Lowpass>0.5):
        message="The filter frequencies must be between 0 and 0.5"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False

    # Check core parameters
    if thGoodClass<0 or thGoodClass>100:
        message="The good class threshold must be between 0 and 100"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False

    return retval

#		
# Main
#     
if __name__ == '__main__':
    CL2D=CL2D_class(
                 InSelFile,
                 WorkingDir,
                 ProjectDir,
                 LogDir,
                 NumberOfReferences,
                 NumberOfReferences0,
                 NumberOfIterations,
                 DoFilter,
                 Highpass,
                 Lowpass,
                 ComparisonMethod,
                 ClusteringMethod,
                 AdditionalParameters,
                 thGoodClass,
                 thJunkZscore,
                 thPCAZscore,
                 NumberOfMpiProcesses,
                 SystemFlavour)
