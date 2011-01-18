#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for image classification using self-organizing maps
#
# Example use:
# ./xmipp_protocol_kerdensom.py
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
SelFileName='CL2D/classes4/class_core_classes/class_core_00000.stk'
# Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='KerDenSOM/test1'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/media/usbdisk/Experiments/TestProtocols/KerDenSOM/Small_TAG_normalized'
# {hidden} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Mask parameters
#------------------------------------------------------------------------------------------------
# Design your mask graphically?
DoXmask=True
""" This will launch a graphical program to design your own mask.
    Be careful NOT to submit your job via a queueing system!
"""
# {file} OR provide an already existing mask:
MaskFileName=''
#-----------------------------------------------------------------------------
# {section} Classification: classify_kerdensom 
#-----------------------------------------------------------------------------
# X-dimension of the self-organizing map:
SomXdim=7
# Y-dimension of the self-organizing map:
SomYdim=7
# {expert} Initial regularization factor:
""" The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    If the output map is too smooth, lower the regularization factors
    If the output map is not organized, higher the regularization factors
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
SomReg0=1000
# {expert} Final regularization factor:
SomReg1=200
# {expert} Number of steps to lower the regularization factor:
SomSteps=5
# {expert} Additional kerdenSOM parameters:
""" For a complete description 
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraCommand=''
#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results
AnalysisScript='visualize_kerdensom.py'
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
import os,shutil,sys,time
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path

def stepPerformed(step,filename):
    import re
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    expr = re.compile(step)
    return len(filter(expr.search,lines))>0

class kerdensom_class:
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
                DoXmask,
                MaskFileName,
                SomXdim,
                SomYdim,
                SomReg0,
                SomReg1,
                SomSteps,                
                KerdensomExtraCommand):
       import log

       self.SelFileName=SelFileName
       self.WorkingDir=WorkingDir
       self.ProjectDir=ProjectDir
       self.DoXmask=DoXmask
       self.MaskFileName=MaskFileName
       self.SomXdim=SomXdim
       self.SomYdim=SomYdim
       self.SomReg0=SomReg0
       self.SomReg1=SomReg1
       self.SomSteps=SomSteps              
       self.KerdensomExtraCommand=KerdensomExtraCommand
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
                 "DoXmask",
                 "MaskFileName",
                 "SomXdim",
                 "SomYdim",
                 "SomReg0",
                 "SomReg1",
                 "SomSteps",
                 "KerdensomExtraCommand"]);

       #made backup of this script
       log.make_backup_of_script_file(sys.argv[0],self.WorkingDir)

       # Update status
       fh=open(self.WorkingDir + "/status.txt", "a")
       fh.write("Step 0: Process started at " + time.asctime() + "\n")
       fh.close()

       # Effectively run
       if (self.DoXmask):
            self.execute_xmask()
       self.execute_img2data()
       self.execute_kerdensom()
       self.execute_data2img()
  
       fh=open(self.WorkingDir + "/status.txt", "a")
       fh.write("Step F: Process finished at " + time.asctime() + "\n")
       fh.close()
       
   #------------------------------------------------------------------------
   #Auxiliary functions
   #------------------------------------------------------------------------
   def execute_xmask(self):
        import os
        import launch_job
        self.MaskFileName=self.WorkingDir+'/mask_design.msk'
        command=' -sel '+self.SelFileName + ' -save_as '+self.MaskFileName
        launch_job.launch_job("xmipp_mask_design",
                              command,
                              self.log,
                              False,1,1,'')

   def execute_img2data(self):
        import os
        import launch_job
        command=' -i '+ self.SelFileName + \
                ' -o ' + self.WorkingDir+"/data.txt"
        if self.MaskFileName=='':
            command+=' -nomask '
        else:
            command+=' -mask '+ self.MaskFileName
        launch_job.launch_job("xmipp_convert_img2data",
                              command,
                              self.log,
                              False,1,1,'')

   def execute_data2img(self):
        import os
        import launch_job
        command=' -i '+ self.WorkingDir+"/som.cod" + \
                ' -o '+ self.WorkingDir+"/som.stk" + \
                ' -rows ' + str(self.SomXdim) + \
                ' -cols ' + str(self.SomYdim)
        if self.MaskFileName=='':
            command+=' -nomask '
        else:
            command+=' -mask '+ self.MaskFileName
        launch_job.launch_job("xmipp_convert_data2img",
                              command,
                              self.log,
                              False,1,1,'')

   #------------------------------------------------------------------------
   #execute_KerDenSOM
   #------------------------------------------------------------------------
   def execute_kerdensom(self):
      import launch_job

      if stepPerformed("Step 1",self.WorkingDir + "/status.txt"):
         return
      print '*********************************************************************'
      print '* Computing kerdensom ...'
      command=' -v 1 -i '  + self.WorkingDir+"/data.txt" + \
              ' -o '    + self.WorkingDir+"/som"  + \
              ' --xdim ' + str(self.SomXdim) + \
              ' --ydim ' + str(self.SomYdim) + \
              ' --reg0 ' + str(self.SomReg0) + \
              ' --reg1 ' + str(self.SomReg1) + \
              ' --steps ' + str(self.SomSteps) + \
              ' '  + str(self.KerdensomExtraCommand)
      launch_job.launch_job("xmipp_classify_kerdensom",
                            command,
                            self.log,
                            False,1,1,'')
      os.system("rm -f "+self.WorkingDir+"/som_*")
      if os.path.exists(self.WorkingDir+"/som.cod"):
          fh=open(self.WorkingDir + "/status.txt", "a")
          fh.write("Step 1: KerDenSOM finished at " + time.asctime() + "\n")
          fh.close()

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
    if not os.path.exists(SelFileName):
        message="The input selfile is not valid"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    return retval

#
# main
#     
if __name__ == '__main__':
   mykerdensom=kerdensom_class(SelFileName,
                                      WorkingDir,
                                      ProjectDir,
                                      DoXmask,
                                      MaskFileName,
                                      SomXdim,
                                      SomYdim,
                                      SomReg0,
                                      SomReg1,
                                      SomSteps,                
                                      KerdensomExtraCommand)
