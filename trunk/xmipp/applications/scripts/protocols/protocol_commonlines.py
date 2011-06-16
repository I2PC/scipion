#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for building initial references by common lines
# This protocol is based on EMAN 1
#
# Example use:
# ./xmipp_protocol_commonlines.py
#
# Author:Carlos Oscar Sorzano, January 2011
#
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# {file} Selfile or stack with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName='CL2D/classes8/class_level_00.stk'
# Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='CommonLines/test1'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/media/usbdisk/Experiments/TestProtocols/CommonLines/Tercera_tanda_LTag_RPA'
# {hidden} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Common lines parameters
#------------------------------------------------------------------------------------------------
# Particle radius
Radius=25
""" Maximum radius of the particle
"""
# Symmetry
""" c1=No symmetry; c2, c3, ...
"""
Symmetry='c1'
#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use distributed-memory parallelization (MPI)?
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
DoParallel=False
# Number of MPI processes to use:
NumberOfMpiProcesses=1
# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. 
"""
SystemFlavour=''

#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results
AnalysisScript='visualize_commonlines.py'
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
import os,shutil,sys,time
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path

class commonline_class:
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
                SelFileName,
                WorkingDir,
                ProjectDir,
                Radius,
                Symmetry,
                DoParallel,
                NumberOfMpiProcesses,
                SystemFlavour):
       import log

       self.SelFileName=SelFileName
       self.WorkingDir=WorkingDir
       self.ProjectDir=ProjectDir
       self.Radius=Radius
       self.Symmetry=Symmetry
       self.DoParallel=DoParallel
       self.NumberOfMpiProcesses=NumberOfMpiProcesses
       self.SystemFlavour=SystemFlavour
       self.log=log.init_log_system(ProjectDir,
                                      LogDir,
                                      sys.argv[0],
                                      WorkingDir)
       
       # Create directory if does not exist
       if not os.path.exists(self.WorkingDir):
           os.makedirs(self.WorkingDir)

       # Save parameters and compare to possible previous runs
       self.saveAndCompareParameters([
                 "SelFileName",
                 "Radius",
                 "Symmetry"]);

       #made backup of this script
       log.make_backup_of_script_file(sys.argv[0],self.WorkingDir)

       # Update status
       fh=open(self.WorkingDir + "/status.txt", "a")
       fh.write("Step 0: Process started at " + time.asctime() + "\n")
       fh.close()

       # Effectively run
       self.execute_commonlines()
  
       fh=open(self.WorkingDir + "/status.txt", "a")
       fh.write("Step F: Process finished at " + time.asctime() + "\n")
       fh.close()
       
   #------------------------------------------------------------------------
   # Common lines
   #------------------------------------------------------------------------
   def execute_commonlines(self):
      import launch_job

      print '*********************************************************************'
      print '* Computing the common lines ...'
      params=' -i '+self.SelFileName+' -o '+self.WorkingDir+"/inputImages.hed"
      launchJob("xmipp_convert_image",
                            params,
                            self.log,
                            False,1,1,'')
      currentDir=os.getcwd()
      os.chdir(self.WorkingDir)
      params="inputImages.hed mask="+str(self.Radius)
      launchJob("cenalignint",
                            params,
                            self.log,
                            False,1,1,'')
      os.system("rm -f avg.* inputImages.*")
      params="ali.hed mask="+str(self.Radius)+" rounds=5"
      if self.DoParallel:
          params+=" proc="+str(self.NumberOfMpiProcesses)
      if self.Symmetry!="c1":
          params+=" sym "+self.Symmetry
      launchJob("startAny",
                            params,
                            self.log,
                            False,1,1,'')
      os.chdir(currentDir)
      
      if os.path.exists(self.WorkingDir+"/threed.0a.mrc"):
          fh=open(self.WorkingDir + "/status.txt", "a")
          fh.write("Step 1: Common line volume finished at " + time.asctime() + "\n")
          fh.close()

# Preconditions
def checkErrors():
    errors = []
    # Check if there is workingdir
    if WorkingDir == "":
        errors.append("No working directory given")
    # Check that there are any micrograph to process
    if not os.path.exists(SelFileName):
        errors.append("The input selfile is not valid")
    # Check that Eman is accessible
    import which
    startAny=which.which('startAny')
    if startAny=='':
        errors.append("EMAN is not accesible")
        
    return errors

#
# main
#     
if __name__ == '__main__':
   commonline_class(SelFileName,
                    WorkingDir,
                    ProjectDir,
                    Radius,
                    Symmetry,
                    DoParallel,
                    NumberOfMpiProcesses,
                    SystemFlavour)
