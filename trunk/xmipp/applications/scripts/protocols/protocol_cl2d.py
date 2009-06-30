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
InSelFile='imgs1000.sel'
# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='CL2D/CL2ref'
# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/gpfs/fs1/home/bioinfo/coss/Tilted_and_Untilted'
# {expert} Directory name for logfiles:
""" All logfiles will be stored here
"""
LogDir='Logs'
#------------------------------------------------------------------------------------------------
# {section} class_averages parameters
#------------------------------------------------------------------------------------------------
# Number of references (or classes) to be used:
NumberOfReferences=2

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=10

# Also include mirror transformation in the alignment?
"""  Including the mirror transformation is useful if your particles have a handedness and may fall either face-up or face-down on the grid.
 Note that when you want to use this CL2D run for later RCT reconstruction, you can NOT include the mirror transformation here.
"""
DoMirror=True

# {expert} Use the fast version of this algorithm?
DoFast=True

# {expert} Additional parameters for class_averages
""" -minsize <N>, -codes0 <N>
"""
AdditionalParameters=''

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of MPI processes to use:
NumberOfMpiProcesses=8
# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. 
"""
SystemFlavour=''

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
class CL2D_class:

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 NumberOfReferences,
                 NumberOfIterations,
                 DoMirror,
                 DoFast,
                 AdditionalParameters,
                 NumberOfMpiProcesses,
                 SystemFlavour):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,selfile

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.NumberOfReferences=NumberOfReferences
        self.NumberOfIterations=NumberOfIterations
        self.DoMirror=DoMirror
        self.DoFast=DoFast
        self.AdditionalParameters=AdditionalParameters
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
        self.InSelFile=os.path.abspath(self.WorkingDir+'/'+InSelFile)
        newsel.write(self.InSelFile)

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
            os.path.abspath(self.WorkingDir))

        # Execute protocol in the working directory
        self.execute_CLalign2D()
     
        # Finish
        self.close()

    def execute_CLalign2D(self):
        import os
        import launch_job
        print '*********************************************************************'
        print '*  Executing class_averages program :' 
        params= '-i '+str(self.InSelFile)+' -o '+WorkingDir+'/class '+\
                ' -codes '+str(self.NumberOfReferences)+\
                ' -iter '+str(self.NumberOfIterations)
        params+=' '+self.AdditionalParameters
        if (self.DoFast):
            params+= ' -fast '
        if (not self.DoMirror):
            params+= ' -no_mirror '

        program="xmipp_class_averages"
           
        launch_job.launch_job(program,
                              params,
                              self.log,
                              True,
                              self.NumberOfMpiProcesses,
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
                    NumberOfReferences,
                    NumberOfIterations,
                    DoMirror,
                    DoFast,
                    AdditionalParameters,
                    NumberOfMpiProcesses,
                    SystemFlavour)
