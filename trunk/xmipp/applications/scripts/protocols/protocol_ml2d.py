#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles, according to:
# {please cite} for ML2D:  Scheres et al. (2005) J.Mol.Biol 348, 139-149 
# {please cite} for MLF2D: Scheres et al. (2007) Structure 15, 1167-1177
#
# Example use:
# ./xmipp_protocol_ml2d.py
#
# Author: Sjors Scheres, January 2008
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {file} Selfile with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
InSelFile="all_images.sel"
# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir="ML2D/ML3ref"
# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir="/home/scheres/xmipp/applications/scripts/protocols"
# {expert} Directory name for logfiles:
""" All logfiles will be stored here
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} MLF-specific parameters
#------------------------------------------------------------------------------------------------
# Perform MLF2D instead of ML2D classification?
DoMlf=False
# {file} CTFdat file with the input images:
""" The names of both the images and the ctf-parameter files should be with absolute paths.
"""
InCtfDatFile="all_images.ctfdat"
# High-resolution limit (in Angstroms)
""" No frequencies higher than this limit will be taken into account. If zero is given, no limit is imposed
"""
HighResLimit=20
#------------------------------------------------------------------------------------------------
# {section} ml(f)_align2d parameters
#------------------------------------------------------------------------------------------------
# Perform 2D maximum-likelihood refinement?
DoML2D=True
# Number of references (or classes) to be used:
NumberOfReferences=3
# Also include mirror transformation in the alignment?
"""  Including the mirror transformation is useful if your particles have a handedness and may fall either face-up or face-down on the grid.
 Note that when you want to use this ML2D run for later RCT reconstruction, you can NOT include the mirror transformation here.
"""
DoMirror=False
# Use the fast version of this algorithm?
""" See Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244:
    http://dx.doi.org/10.1093/bioinformatics/bti1140
"""
DoFast=True
# Refine the normalization parameters for each image?
""" This variant of the algorithm deals with normalization errors. For more info see (and please cite) Scheres et. al. (2009) J. Struc. Biol., in press
"""
DoNorm=False
# {expert} Restart after iteration:
""" For previous runs that stopped before convergence, resume the calculations
    after the completely finished iteration. (Use zero to start from the beginning)
"""
RestartIter=0
# {expert} Additional xmipp_ml_align2d parameters:
""" For a complete description see the ml_align2d manual page at:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/MLalign2D
"""
ExtraParamsMLalign2D=""
#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads=1
# Use distributed-memory parallelization (MPI)?
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
DoParallel=False
# Number of MPI processes to use:
NumberOfMpiProcesses=5
# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. Or read: xxx
"""
SystemFlavour="SLURM-MPICH"
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
                 DoMlf,
                 InCtfDatFile,
                 HighResLimit,
                 DoML2D,
                 NumberOfReferences,
                 DoMirror,
                 DoFast,
                 DoNorm,
                 RestartIter,
                 ExtraParamsMLalign2D,
                 NumberOfThreads,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,selfile

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.DoMlf=DoMlf
        self.HighResLimit=HighResLimit
        self.NumberOfReferences=NumberOfReferences
        self.DoMirror=DoMirror
        self.DoFast=DoFast
        self.DoNorm=DoNorm
        self.ExtraParamsMLalign2D=ExtraParamsMLalign2D
        self.NumberOfThreads=NumberOfThreads
        self.DoParallel=DoParallel
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour
   
        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # This is not a restart
        if (RestartIter < 1):
            # Delete working directory if it exists, make a new one
            if (DoDeleteWorkingDir and DoML2D): 
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

            if (self.DoMlf):
                # Copy CTFdat to the workingdir as well
                shutil.copy(InCtfDatFile,self.WorkingDir+'/my.ctfdat')

            # Backup script
            log.make_backup_of_script_file(sys.argv[0],
                                               os.path.abspath(self.WorkingDir))
    
            # Execute protocol in the working directory
            os.chdir(self.WorkingDir)
            if (DoML2D):
                self.execute_MLalign2D()

        # Restarting a previous run...
        else:
            # Execute protocol in the working directory
            os.chdir(self.WorkingDir)
            self.restart_MLalign2D(RestartIter)
     
        # Return to parent dir
        os.chdir(os.pardir)

    def execute_MLalign2D(self):
        import os
        import launch_job
        print '*********************************************************************'
        print '*  Executing ml(f)_align2d program :' 
        params= ' -o ml2d -i '    + str(self.InSelFile) + \
                ' -nref ' + str(self.NumberOfReferences)
        params+=' '+self.ExtraParamsMLalign2D
        if (self.DoFast):
            params+= ' -fast '
        if (self.DoNorm):
            params+= ' -norm '
        if (self.DoMirror):
            params+= ' -mirror '
        if (self.NumberOfThreads > 1):
            params+= ' -thr ' + str(self.NumberOfThreads)
        if (self.DoMlf):
            params+= ' -ctfdat my.ctfdat'
            if (self.HighResLimit > 0):
                params += ' -high ' + str(self.HighResLimit)

        if (self.DoMlf):
            program="xmipp_mlf_alignd"
        else:
            program="xmipp_ml_align2d"
           
        launch_job.launch_job(program,
                              params,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)

    def restart_MLalign2D(self, iter):
        import os
        import launch_job, utils_xmipp
        print '*********************************************************************'
        print '*  Restarting ml(f)_align2d program :' 
        params= ' -restart ' + utils_xmipp.composeFileName('ml2d_it',iter,'log')

        if (self.DoMlf):
            program="xmipp_mlf_alignd"
        else:
            program="xmipp_ml_align2d"
           
        launch_job.launch_job(program,
                              params,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)

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
                    DoMlf,
                    InCtfDatFile,
                    HighResLimit,
                    DoML2D,
                    NumberOfReferences,
                    DoMirror,
                    DoFast,
                    DoNorm,
                    RestartIter,
                    ExtraParamsMLalign2D,
                    NumberOfThreads,
                    DoParallel,
                    NumberOfMpiProcesses,
                    SystemFlavour)
    # close 
    ML2D.close()

