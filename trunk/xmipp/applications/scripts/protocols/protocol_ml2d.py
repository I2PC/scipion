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
InSelFile="all_images.sel"
# Working subdirectory:
WorkingDir="ML2D/ML3ref"
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
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
# Use multiple processors in parallel?
DoParallel=False
# Number of processors to use:
MyNumberOfCPUs=10
# {file} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, this file may be required. If not, just leave this entry blank.
    If your job submission system uses an environment variable, just type it here with the leading $
"""
MyMachineFile="/home2/bioinfo/scheres/machines.dat"
# {expert} Control file
""" This is an ugly solution to have improved killing control over the mpi jobs.
    The script will create this file, and any mpi-job will be killed as soon as this file doesn't exist anymore. This is required with certain queueing systems.
"""
MyControlFile=""
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
                 RestartIter,
                 ExtraParamsMLalign2D,
                 DoParallel,
                 MyNumberOfCPUs,
                 MyMachineFile,
                 MyControlFile):
	     
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
        self.ExtraParamsMLalign2D=ExtraParamsMLalign2D
        self.DoParallel=DoParallel
        self.MyNumberOfCPUs=MyNumberOfCPUs
        if (MyMachineFile[0]=="$"):
            self.MyMachineFile=MyMachineFile
        else:
            self.MyMachineFile=os.path.abspath(MyMachineFile)
        if (MyControlFile==""):
            self.DoControl=False
        else:
            self.DoControl=True
   
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

            # Create a CONTROL file for improved killing control
            if (self.DoControl):
                self.MyControlFile=os.path.abspath(self.WorkingDir+
                                                   '/'+MyControlFile)
                FILE = open(self.MyControlFile,"w")
                FILE.write("Delete this file to kill all current processes\n")
                FILE.close()

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
            # Create a CONTROL file for improved killing control
            if (self.DoControl):
                self.MyControlFile=os.path.abspath(self.WorkingDir+
                                                   '/'+MyControlFile)
                FILE = open(self.MyControlFile,"w")
                FILE.write("Delete this file to kill current processes\n")
                FILE.close()

            # Execute protocol in the working directory
            os.chdir(self.WorkingDir)
            self.restart_MLalign2D(RestartIter)
     
        # Return to parent dir
        os.chdir(os.pardir)

        # Delete the CONTROL file for improved killing control
        if (self.DoControl):
            os.remove(self.MyControlFile)

    def execute_MLalign2D(self):
        import os
        import launch_parallel_job
        print '*********************************************************************'
        print '*  Executing ml(f)_align2d program :' 
        params= ' -o ml2d -i '    + str(self.InSelFile) + \
                ' -nref ' + str(self.NumberOfReferences)
        params+=' '+self.ExtraParamsMLalign2D
        if (self.DoFast):
            params+= ' -fast '
        if (self.DoMirror):
            params+= ' -mirror '
        if (self.DoMlf):
            params+= ' -ctfdat my.ctfdat'
            if (self.HighResLimit > 0):
                params += ' -high ' + str(self.HighResLimit)
        if (self.DoControl):
            params+=' -control ' + self.MyControlFile

        if (self.DoMlf):
            program="xmipp_ml_alignd"
            mpiprogram="xmipp_mpi_mlf_align2d"
        else:
            program="xmipp_ml_align2d"
            mpiprogram="xmipp_mpi_ml_align2d"
           
        launch_parallel_job.launch_job(self.DoParallel,
                                       program,
                                       mpiprogram,
                                       params,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)

    def restart_MLalign2D(self, iter):
        import os
        import launch_parallel_job
        print '*********************************************************************'
        print '*  Restarting ml(f)_align2d program :' 
        params= ' -restart ml2d_it'    + str(iter).zfill(5) + '.log'

        if (self.DoMlf):
            program="xmipp_ml_alignd"
            mpiprogram="xmipp_mpi_mlf_align2d"
        else:
            program="xmipp_ml_align2d"
            mpiprogram="xmipp_mpi_ml_align2d"
           
        launch_parallel_job.launch_job(self.DoParallel,
                                       program,
                                       mpiprogram,
                                       params,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)

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
                    RestartIter,
                    ExtraParamsMLalign2D,
                    DoParallel,
                    MyNumberOfCPUs,
                    MyMachineFile,
                    MyControlFile)
    # close 
    ML2D.close()

