#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for building initial references by common lines
# This protocol is based on EMAN 1
#
# Example use:
# ./xmipp_protocol_commonlines.py
#
# Author:Carlos Oscar Sorzano, January 2011
#
import os,shutil,time
from protlib_base import *
from config_protocols import protDict

class ProtCommonLines(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.commonlines.name, scriptname, project)
        self.Import = 'from protocol_commonlines import *'
    
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
#   def __init__(self,
#                SelFileName,
#                WorkingDir,
#                ProjectDir,
#                Radius,
#                Symmetry,
#                DoParallel,
#                NumberOfMpi,
#                SystemFlavour):
#       import log
#
#       self.SelFileName=SelFileName
#       self.WorkingDir=WorkingDir
#       self.ProjectDir=ProjectDir
#       self.Radius=Radius
#       self.Symmetry=Symmetry
#       self.DoParallel=DoParallel
#       self.NumberOfMpi=NumberOfMpi
#       self.SystemFlavour=SystemFlavour
#       self.log=log.init_log_system(ProjectDir,
#                                      LogDir,
#                                      sys.argv[0],
#                                      WorkingDir)
#       
#       # Create directory if does not exist
#       if not os.path.exists(self.WorkingDir):
#           os.makedirs(self.WorkingDir)
#
#       # Save parameters and compare to possible previous runs
#       self.saveAndCompareParameters([
#                 "SelFileName",
#                 "Radius",
#                 "Symmetry"]);
#
#       #made backup of this script
#       log.make_backup_of_script_file(sys.argv[0],self.WorkingDir)
#
#       # Update status
#       fh=open(self.WorkingDir + "/status.txt", "a")
#       fh.write("Step 0: Process started at " + time.asctime() + "\n")
#       fh.close()
#
#       # Effectively run
#       self.execute_commonlines()
#  
#       fh=open(self.WorkingDir + "/status.txt", "a")
#       fh.write("Step F: Process finished at " + time.asctime() + "\n")
#       fh.close()
       
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
          params+=" proc="+str(self.NumberOfMpi)
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
    def validate():
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

