#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based manual particle picking
#
# For each micrograph in the PreprocessingDir, this program will launch
# the xmipp_mark program 
# A graphical interface exists to identify micrographs that have been finished
#
# Author: Sjors Scheres, March 2007
# Author: Carlos Oscar Sorzano, June, 2011
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {begin_of_header} 
#------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------
# Comment
""" Specify the particularities of this run """
Comment=''

# Working subdirectory:
""" Working directory for this protocol 
"""
RunName = "particles_001"

# {dir} Preprocessing dir
""" Directory with the preprocessing (output of the Preprocessing Micrographs protocol)
"""
PreprocessingDir = "Preprocessing/micrographs_001"

# Perform automatic particle picking
""" Perform automatic particle picking """
AutomaticPicking = False

#------------------------------------------------------------------------------------------
# {section}{condition}(AutomaticPicking=True) Parallelization issues for automatic particle picking
#------------------------------------------------------------------------------------------
# Number of MPI processes to use:
NumberOfMpiProcesses = 1

#Submmit to queue
"""Submmit to queue"""
SubmmitToQueue=False
#------------------------------------------------------------------------------------------------
# {section}{expert}{condition}(SubmmitToQueue=True) Queue 
#------------------------------------------------------------------------------------------------

# Queue name
"""Name of the queue to submit the job"""
QueueName="default"

# Queue hours
"""This establish a maximum number of hours the job will
be running, after that time it will be killed by the 
queue system"""
QueueHours=72
#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# Do not repeat already taken steps
IsIter=False
ContinueAtIteration=1

# Debugging 
Verify=True                 # Check that some output files are created.
PrintWrapperCommand=True    # Print wrapper name
PrintWrapperParameters=True # Print wrapper parameters
ViewVerifyedFiles=True      # Show file verification 

import os,shutil,sys,threading,time
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path
import xmipp
from protlib_gui import *
from protlib_base import *
from protlib_filesystem import *

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(XmippProtocol):
    def __init__(self, scriptname,project=None):
        super(ProtParticlePicking,self).__init__(protDict.particle_pick.key, scriptname, RunName, project)
        self.Import = 'from xmipp_protocol_particle_pick import *'

    def defineActions(self):
        # Create link to input micrographs
        micrographSelfile=os.path.join(PreprocessingDir,"micrographs.sel")
        self.Db.insertAction('createLinkToMicrographs', [micrographSelfile], MicrographSelfile=micrographSelfile, WorkingDir=self.WorkingDir)       

        # Launch GUI
        self.Db.insertAction('launchParticlePickingGUI', AutomaticPicking=AutomaticPicking, ProjectDir=self.projectDir,
                             PreprocessingDir=PreprocessingDir, WorkingDir=self.WorkingDir)       

    def summary(self):
        pass
    
    def validate(self):
        errors = []

        # Check that there is a valid list of micrographs
        MicrographSelfile=PreprocessingDir+"/micrographs.sel"
        if not os.path.exists(MicrographSelfile)>0:
            errors.append("Cannot find "+MicrographSelfile)
        
        # Check that all micrographs exist
        mD = xmipp.MetaData()
        xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, mD)
        isPairList = mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                     not xmipp.FileName(MicrographSelfile).isStar1()
        NnotFound=0
        message=""
        for id in mD:
             micrograph = mD.getValue(xmipp.MDL_IMAGE,id)
             if not os.path.exists(micrograph):
                message+="  "+micrograph
                NnotFound=NnotFound+1
             if isPairList:
                 micrograph = mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
                 if not os.path.exists(micrograph):
                     message+="  "+micrograph
                     NnotFound=NnotFound+1
        
        if NnotFound>0:
            errors.append("Cannot find the following micrographs: "+message)
                
        # Check that automatic particle picking is not for tilted
        if isPairList and AutomaticPicking:
            errors.append("Automatic particle picking cannot be done on tilt pairs")
        
        return errors

# GUI For marking micrographs
class ProtParticlePickingGUI(BasicGUI):
    def __init__(self,log,AutomaticPicking,ProjectDir,PreprocessingDir,WorkingDir):
        self.log=log
        self.AutomaticPicking=AutomaticPicking
        self.ProjectDir=ProjectDir
        self.PreprocessingDir=PreprocessingDir
        self.WorkingDir=WorkingDir

        # Read micrographs
        self.mD=xmipp.MetaData();
        self.MicrographSelfile=PreprocessingDir+"/micrographs.sel"
        xmipp.readMetaDataWithTwoPossibleImages(self.MicrographSelfile, self.mD)
        self.IsPairList = self.mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                          not xmipp.FileName(self.MicrographSelfile).isStar1()
    
    def createGUI(self):
        self.createBasicGUI()
        self.createScrollableCanvas()
 
    def fillGUI(self):
        # Init GUI variables
        self.whichmark=StringVar()
        self.row={}
        if self.AutomaticPicking:
            self.selectedForAutomaticPickingAuto=[]
            self.selectedForAutomaticPickingName=[]
            self.selectedForAutomaticPickingMark=[]
            titleSpan=5
        else:
            titleSpan=3
        if self.IsPairList:
            self.whichtilted={}

        # Window title
        self.master.title("GUI for Xmipp particle picking")
        
        # Frame header
        self.addLabel("Project directory: "+self.ProjectDir,0,0,titleSpan,fgColor="medium blue",sticky=W)
        self.addLabel("Preprocessing directory: "+self.PreprocessingDir,1,0,titleSpan,fgColor="medium blue",sticky=W)
        if (AutomaticPicking):
            Label(self.frame, text="Manual").grid(row=1,column=3)
            Label(self.frame, text="Auto").grid(row=1,column=4)
 
        # Add all micrographs
        containsEnable=self.mD.containsLabel(xmipp.MDL_ENABLED)
        for id in self.mD:
            micrograph=self.mD.getValue(xmipp.MDL_IMAGE,id)
            if containsEnable:
                if self.mD.getValue(xmipp.MDL_ENABLED,id)==0:
                    continue
            if self.IsPairList:
                self.whichtilted[micrograph]=self.mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
            self.GuiAddMarkEntry(micrograph)

        # Add blue line surrounded by two empty rows 
        self.addLabel("")
        self.addLine("medium blue",0,6)
        self.addLabel("")

        # Add some buttons
        self.buttonrow=(self.frame.grid_size()[1]+1)
        self.addButton("Show preprocessing",self.buttonrow,0,self.ShowPreprocessing,sticky=E)
        if (AutomaticPicking):
            self.addButton("Invert Selection",self.buttonrow,1,self.InvertSelection)
            nextColumn=2
        else:
            nextColumn=1
        self.addButton("Update Total Count",self.buttonrow,nextColumn,self.GuiUpdateCount,binding='<Control_L><U>')
        nextColumn+=1
        if (AutomaticPicking):
            self.addButton("AutoSelect",self.buttonrow+1,1,self.AutomaticallyDetect)

        # Update counts
        self.GuiUpdateCount()
        
    def GuiAddMarkEntry(self,micrograph):
        row=self.frame.grid_size()[1]
        self.row[micrograph]=row

        # Add Micrograph name
        label=micrograph.split("/")[-2]
        if self.IsPairList:
            tiltedName=self.whichtilted[micrograph].split("/")[-2]
            label+=' : '+tiltedName
        self.addLabel(label, row, 0, sticky=E)
        
        # If automatic particle picking, add checkbox
        if AutomaticPicking:
            self.selectedForAutomaticPickingAuto.append(self.addCheckButton("Auto", row, 1, 0, self.AutoSelectionChanged, NW))
            self.selectedForAutomaticPickingName.append(micrograph)
            nextColumn=2
        else:
            nextColumn=1
        
        # Add Mark button
        radio=self.addRadioButton("Mark", row, nextColumn, self.whichmark, micrograph, self.LaunchSingleMark, N)
        if AutomaticPicking:
            self.selectedForAutomaticPickingMark.append(radio)

    def InvertSelection(self):
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            self.selectedForAutomaticPickingAuto[i].set(1-
                self.selectedForAutomaticPickingAuto[i].get());
        self.AutoSelectionChanged();

    def AutoSelectionChanged(self):
        for i in range(0,len(self.selectedForAutomaticPickingTrain)):
            if (self.selectedForAutomaticPickingAuto[i].get()):
                self.selectedForAutomaticPickingMark[i].config(state=DISABLED)
            else:
                self.selectedForAutomaticPickingMark[i].config(state=NORMAL)

    def AutomaticallyDetect(self):
        command_file = open(WorkingDir+"/pick.sh", "w")
        directoryPreprocessing,dummy=os.path.split(self.MicrographSelfile)
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            if (self.selectedForAutomaticPickingAuto[i].get()):
               directory,micrograph=os.path.split(
                  self.selectedForAutomaticPickingName[i])
               filename=WorkingDir+"/"+micrograph+".Common.auto.pos"
               command_file.write(
                  "( xmipp_micrograph_mark -i "+directoryPreprocessing+"/"+\
                  self.selectedForAutomaticPickingName[i]+\
                  " --auto "+WorkingDir+"/Common.auto --autoSelect "+\
                  " --outputRoot "+WorkingDir+"/"+micrograph+\
                  "; if [ -e " + filename + ' ]; then ' + \
                  'echo "Step A: "'+micrograph+' automatically marked on `date` >> ' + \
                  WorkingDir + "/status.txt; fi )\n")
                 
        command_file.write("MPI_BARRIER\n")
        command_file.write('echo "Step F: " finished marking on `date` >> ' + \
                  WorkingDir + "/status.txt \n")
        command_file.close();

        import tkMessageBox
        answer=tkMessageBox._show("Execute autodetection",
                                  "Use a job queueing system?",
                                  tkMessageBox.QUESTION, 
                                  tkMessageBox.YESNO)
        import launch_job
        command=""
        if DoParallel:
            command=launchJob("xmipp_run",
                                 "-i "+WorkingDir+"/pick.sh",
                                 self.log,
                                 True,
                                 NumberOfMpiProcesses,
                                 1,
                                 SystemFlavour,
                                 answer=="yes" or answer==True)
        else:
            os.system("chmod 755 "+WorkingDir+"/pick.sh");
            command=launchJob(WorkingDir+"/pick.sh",
                                 "",
                                 self.log,
                                 False,
                                 NumberOfMpiProcesses,
                                 1,
                                 SystemFlavour,
                                 answer=="yes" or answer==True)
        if (answer=="yes" or answer==True):
            import xmipp_protocol_gui
            python_file = open(WorkingDir+"/pick.py", "w")
            python_file.write("WorkingDir='"+WorkingDir+"'\n")
            python_file.write("DoParallel="+str(DoParallel)+"\n")
            python_file.write("NumberOfMpiProcesses="+str(NumberOfMpiProcesses)+"\n")
            python_file.write("# {end-of-header}\n")
            python_file.write("import os\n")
            python_file.write('os.system("'+command+'")\n');
            python_file.close();

            d = xmipp_protocol_gui.MyQueueLaunch(self.master,"python "+WorkingDir+"/pick.py")
        else:
            os.system(command);

    def GuiUpdateCount(self):
        total_manual=0
        total_auto=0
        if AutomaticPicking:
            familyColumn=2
        else:
            familyColumn=3

        # Count all micrographs
        for micrograph,row in self.row.items():
            manual=self.CountPicked(micrograph,"Common")
            total_manual+=manual
            self.addLabel(str(manual).zfill(5),row,familyColumn, sticky=N+S+W+E)
            if AutomaticPicking:
                auto=self.CountPicked(micrograph,"Common.auto")
                total_auto+=auto
                self.addLabel(str(auto).zfill(5),row,familyColumn+1, sticky=N+S+W+E)

        # Add summaries
        self.addLabel(str(total_manual).zfill(5),self.buttonrow,familyColumn, sticky=N+S+W+E)
        if AutomaticPicking:
            self.addLabel(str(total_auto).zfill(5),self.buttonrow,familyColumn+1, sticky=N+S+W+E)

    def CountPicked(self,micrograph,label):
        micrographName=micrograph.split("/")[-2]
        posfile=self.WorkingDir+"/"+str(micrographName)+'.'+label+'.pos'
        if os.path.exists(posfile):
            mD=xmipp.MetaData(posfile);
            return mD.size()
        return 0
    
    def ShowPreprocessing(self):
        self.GuiUpdateCount()
        command='xmipp_visualize_preprocessing_micrographj -i '+self.MicrographSelfile+" &"
        os.system(command)

    def LaunchSingleMark(self):
        self.GuiUpdateCount()
        micrograph=self.whichmark.get()
        if micrograph=='':
            return
        micrographName=micrograph.split("/")[-2]
        
        arguments=' -i '+micrograph+" --outputRoot "+self.WorkingDir+"/"+micrographName
        if (AutomaticPicking):
            arguments+=' --auto '+self.WorkingDir+"/Common.auto"
            verifyFile=[self.WorkingDir+"/"+micrographName+".Common.auto.pos"]
        else:
            verifyFile=[self.WorkingDir+"/"+micrographName+".Common.pos"]
        print arguments
        RunJobThread(self.log,"xmipp_micrograph_mark",arguments,verifyFile)

    def LaunchPairMark(self):
        self.GuiUpdateCount()
        untilted=self.whichmark.get()
        tilted=self.whichtilted[untilted]
        uname=untilted.split("/")[-2]
        tname=tilted.split("/")[-2]
        arguments=' -i '+untilted+" --tilted "+tilted+" --outputRoot "+self.WorkingDir+"/"+uname
        verifyFile=[self.WorkingDir+"/"+uname+".Common.pos"]
        RunJobThread(self.log,"xmipp_micrograph_mark",arguments,verifyFile).start()

# RunJobThread
class RunJobThread(threading.Thread):
    def __init__(self,log,program,arguments,verifyFiles):
        self.log=log
        self.program=program
        self.arguments=arguments
        self.verifyFiles=verifyFiles
    
    def run(self):
        runJob(self.log,self.program,self.arguments,parameters,False,1,1,"")
        OK=True
        for file in verifyFiles:
            if not os.path.exists(verifyFile):
                OK=False

# Actions
def createLinkToMicrographs(log,WorkingDir,MicrographSelfile):
    inputSelfile=os.path.join(WorkingDir,"micrographs.sel")
    if not os.path.exists(inputSelfile):
        os.system("ln -s "+os.path.relpath(os.path.abspath(MicrographSelfile),WorkingDir)+" "+inputSelfile)

# Execute protocol in the working directory
def launchParticlePickingGUI(log,AutomaticPicking,ProjectDir,PreprocessingDir,WorkingDir):
    gui=ProtParticlePickingGUI(log,AutomaticPicking,ProjectDir,PreprocessingDir,WorkingDir)
    gui.createGUI()
    gui.fillGUI()
    gui.launchGUI()

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
