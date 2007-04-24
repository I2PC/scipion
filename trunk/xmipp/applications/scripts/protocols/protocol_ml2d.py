#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles, according to:
# {please cite} Scheres et al. (2005) J.Mol.Biol 348, 139-149 
# {please cite} Scheres et al. (2005) Bioinformatics, 21(suppl2), ii243-244 (fast version)
#
# Example use:
# ./protocol_ml2d.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Selfile with the input images (relative path from ProjectDir):
InSelFile="all_images.sel"
# Working subdirectory:
WorkingDir="ML3ref"
# Delete working subdirectory if it already exists?
""" The directory will not be deleted when only visualizing! 
"""
DoDeleteWorkingDir=True
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
""" All logfiles will be stored in $ProjectDir/$LogDir
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} ml_align2d parameters
#------------------------------------------------------------------------------------------------
# Perform 2D maximum-likelihood refinement?
DoML2D=False
# Number of references (or classes) to be used:
NumberOfReferences=3
# Also include mirror transformation in the alignment?
"""  Including the mirror transformation is useful if your particles have a handedness
     and may fall either face-up or face-down on the grid
"""
DoMirror=False
# Use the fast version of this algorithm?
""" See Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244
"""
DoFast=True
# {expert} Additional xmipp_ml_align2d parameters:
ExtraParamsMLalign2D=""
#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel? (see Expert options)
DoParallel=False
# Number of processors to use:
MyNumberOfCPUs=10
# A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, this file may be required. If not, just leave this entry blank.
    If your job submission system uses an environment variable, just type it here with the leading $
"""
MyMachineFile="/home2/bioinfo/scheres/machines.dat"
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_ml2d.py"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
class ML2D_class:

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 DoML2D,
                 NumberOfReferences,
                 DoMirror,
                 DoFast,
                 ExtraParamsMLalign2D,
                 DoParallel,
                 MyNumberOfCPUs,
                 MyMachineFile):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.InSelFile=self.ProjectDir+'/'+str(InSelFile)
        self.NumberOfReferences=NumberOfReferences
        self.DoMirror=DoMirror
        self.DoFast=DoFast
        self.ExtraParamsMLalign2D=ExtraParamsMLalign2D
        self.DoParallel=DoParallel
        self.MyNumberOfCPUs=MyNumberOfCPUs
        self.MyMachineFile=MyMachineFile

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Delete working directory if it exists, make a new one, and go there
        if (DoDeleteWorkingDir and DoML2D): 
            if os.path.exists(self.WorkingDir):
                shutil.rmtree(self.WorkingDir)
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Execute MLalign2D in the working directory
        os.chdir(self.WorkingDir)
        if (DoML2D):
            self.execute_MLalign2D()

        # Return to parent dir
        os.chdir(os.pardir)


    def execute_MLalign2D(self):
        import os
        import launch_parallel_job
        print '*********************************************************************'
        print '*  Executing ml_align2d program :' 
        params= ' -i '    + str(self.InSelFile) + \
                ' -o '    + str(self.WorkingDir) + \
                ' -nref ' + str(self.NumberOfReferences)
        if (self.DoFast):
            params+= ' -fast '
        if (self.DoMirror):
            params+= ' -mirror '
        # By default, do not write offsets for 2D alignments...
        # This will affect -restart, but that falls outside the scope of the protocol anyway
        params+=' '+self.ExtraParamsMLalign2D
        params+=' -dont_write_offsets'

        launch_parallel_job.launch_job(self.DoParallel,
                                       "xmipp_ml_align2d",
                                       "xmipp_mpi_ml_align2d",
                                       params,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)

    def visualize_ML2D(self):
        # Visualize class averages:
        import glob
        import os
        selfiles=glob.glob(str(self.WorkingDir)+'_it?????.sel')
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            lastselfile=selfiles[-1]
            command='xmipp_show -sel '+str(lastselfile)+ ' & '
            print '* ',command
            self.log.info(command)
            os.system(command)
            # print logfile to screen:
            logfiles=glob.glob('ml2d_it?????.log')
            lastlogfile=logfiles[-1]
            fh=open(lastlogfile,'r')
            loglines=fh.readlines()
            print "Logfile "+str(lastlogfile)+": "
            for line in loglines:
                print line[:-1]

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    # create ML2D_class object

    ML2D=ML2D_class(InSelFile,
                    WorkingDir,
                    DoDeleteWorkingDir,
                    ProjectDir,
                    LogDir,
                    DoML2D,
                    NumberOfReferences,
                    DoMirror,
                    DoFast,
                    ExtraParamsMLalign2D,
                    DoParallel,
                    MyNumberOfCPUs,
                    MyMachineFile)
    # close 
    ML2D.close()

