#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for rotational spectra classification
# using self-organizing maps
#
# Example use:
# ./xmipp_rotational_spectra.py
#
# Author:Roberto Marabini, March 2007
#        Carlos Oscar Sorzano, January 2011
#
# required files: log.py, spider_header.py, Sel_Files.py
# by defaults searchs for python files in ../python directory
#
import os, shutil, time
from protlib_base import *

def stepPerformed(step,filename):
    import re
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    expr = re.compile(step)
    return len(filter(expr.search,lines))>0

class ProtRotSpectra(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.rotspectra.name, scriptname, project)
        self.Import = 'from protocol_rotspectra import *'    
    
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
                self.mylog.info("Deleting working directory since it is run with different parameters")
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
#                HowCenter,
#                SpectraInnerRadius,
#                SpectraOuterRadius,
#                SpectraLowHarmonic,
#                SpectraHighHarmonic,
#                SomXdim,
#                SomYdim,
#                SomReg0,
#                SomReg1,
#                SomSteps,                
#                KerdensomExtraCommand):
#       import log
#
#       self.SelFileName=SelFileName
#       self.WorkingDir=WorkingDir
#       self.ProjectDir=ProjectDir
#       self.HowCenter=HowCenter
#       self.SpectraInnerRadius=SpectraInnerRadius
#       self.SpectraOuterRadius=SpectraOuterRadius
#       self.SpectraLowHarmonic=SpectraLowHarmonic
#       self.SpectraHighHarmonic=SpectraHighHarmonic
#       self.SomXdim=SomXdim
#       self.SomYdim=SomYdim
#       self.SomReg0=SomReg0
#       self.SomReg1=SomReg1
#       self.SomSteps=SomSteps              
#       self.KerdensomExtraCommand=KerdensomExtraCommand
#       self.mylog=log.init_log_system(ProjectDir,
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
#                 "HowCenter",
#                 "SpectraInnerRadius",
#                 "SpectraOuterRadius",
#                 "SpectraLowHarmonic",
#                 "SpectraHighHarmonic",
#                 "SomXdim",
#                 "SomYdim",
#                 "SomReg0",
#                 "SomReg1",
#                 "SomSteps",
#                 "KerdensomExtraCommand"]);
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
#       self.execute_findcenter()
#       self.execute_spectra()
#       self.execute_KerDenSOM()
#  
#       fh=open(self.WorkingDir + "/status.txt", "a")
#       fh.write("Step F: Process finished at " + time.asctime() + "\n")
#       fh.close()
       
   #------------------------------------------------------------------------
   #execute_findcenter
   #------------------------------------------------------------------------
    
    def execute_findcenter(self):
      import os, string, xmipp
      print '*********************************************************************'
      print '* Looking for the position of the center of symmetry ...'
      fn=xmipp.FileName(self.SelFileName)
      if fn.isMetaData():
          MD=xmipp.MetaData(self.SelFileName)
          R=xmipp.SingleImgSize(MD.getValue(xmipp.MDL_IMAGE))[1]/2
      else:
          R=xmipp.SingleImgSize(self.SelFileName)[1]/2.0
      if self.HowCenter=='Minimize first harmonic':
          command='xmipp_find_center2d'+ \
                  ' -i ' + self.SelFileName + \
                  ' --oroot '+self.WorkingDir+"/center2d" + \
                  ' --r1 ' + str(100.0*self.SpectraInnerRadius/R) + \
                  ' --r2 ' + str(100.0*self.SpectraOuterRadius/R) + \
                  ' --r3 ' + str(100.0*(self.SpectraOuterRadius+2)/R) + \
                  ' --r4 ' + str(100.0*(self.SpectraOuterRadius+5)/R)
          print '* ',command
          self.mylog.info(command)
          os.system(command)
          if os.path.exists(self.WorkingDir+"/center2d_center.xmd"):
              MD=xmipp.MetaData(self.WorkingDir+"/center2d_center.xmd")
              self.xOffset=MD.getValue(xmipp.MDL_X)
              self.yOffset=MD.getValue(xmipp.MDL_Y)
              fh=open(self.WorkingDir + "/status.txt", "a")
              fh.write("Step 1: Center found at " + time.asctime() + "\n")
              fh.close()
      elif self.HowCenter=='Use the middle of the image':
          if fn.isMetaData():
              MD=xmipp.MetaData(fn)
              fnAux=MD.getValue(xmipp.MDL_IMAGE)
              dimensions=xmipp.SingleImgSize(fnAux)
          else:
              dimensions=xmipp.SingleImgSize(fn)
          self.xOffset=dimensions[0]/2.0
          self.yOffset=dimensions[1]/2.0
          fh=open(self.WorkingDir + "/status.txt", "a")
          fh.write("Step 1: Center defined at " + time.asctime() + "\n")
          fh.close()

   #------------------------------------------------------------------------
   #execute_spectra
   #------------------------------------------------------------------------
    def execute_spectra(self):
      import os, string
      import launch_job
      if not stepPerformed("Step 1",self.WorkingDir + "/status.txt"):
         return
      print '*********************************************************************'
      print '* Computing rotational power spectra'
      command=' -i ' + self.SelFileName + \
              ' -o ' + self.WorkingDir+"/rotSpectra.txt" + \
              ' --x0 '  + str(self.xOffset) + \
              ' --y0 '  + str(self.yOffset) + \
              ' --r1 '  + str(self.SpectraInnerRadius) + \
              ' --r2 '   + str(self.SpectraOuterRadius) + \
              ' --low ' + str(self.SpectraLowHarmonic) + \
              ' --high ' + str(self.SpectraHighHarmonic)
      launchJob("xmipp_make_spectra",
                            command,
                            self.mylog,
                            False,1,1,'')
      if os.path.exists(self.WorkingDir+"/rotSpectra.txt"):
          fh=open(self.WorkingDir + "/status.txt", "a")
          fh.write("Step 2: Spectra created at " + time.asctime() + "\n")
          fh.close()
  
   #------------------------------------------------------------------------
   #execute_KerDenSOM
   #------------------------------------------------------------------------
    def execute_KerDenSOM(self):
      import launch_job

      if not stepPerformed("Step 2",self.WorkingDir + "/status.txt"):
         return
      print '*********************************************************************'
      print '* Computing kerdensom ...'
      command=' -v 1 -i '  + self.WorkingDir+"/rotSpectra.txt" + \
              ' -o '    + self.WorkingDir+"/som"  + \
              ' --xdim ' + str(self.SomXdim) + \
              ' --ydim ' + str(self.SomYdim) + \
              ' --reg0 ' + str(self.SomReg0) + \
              ' --reg1 ' + str(self.SomReg1) + \
              ' --steps ' + str(self.SomSteps) + \
              ' '  + str(self.KerdensomExtraCommand)
      launchJob("xmipp_classify_kerdensom",
                            command,
                            self.mylog,
                            False,1,1,'')
      os.system("rm -f "+self.WorkingDir+"/som_*")
      if os.path.exists(self.WorkingDir+"/som.cod"):
          fh=open(self.WorkingDir + "/status.txt", "a")
          fh.write("Step 3: KerDenSOM finished at " + time.asctime() + "\n")
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
        
        return errors
