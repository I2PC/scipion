#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based manual particle picking
#
# For each micrograph in the self.PreprocessingRun, this program will launch
# the xmipp_mark program 
# A graphical interface exists to identify micrographs that have been finished
#
# Author: Carlos Oscar Sorzano, June, 2011
#
#------------------------------------------------------------------------------------------------
# This program is interactive and cannot use MPI or queues
NumberOfMpiProcesses = 1
SubmmitToQueue=False
QueueName="default"
QueueHours=72

# Do not repeat already taken steps
IsIter=False
ContinueAtIteration=1

# Debugging 
Verify=True                 # Check that some output files are created.
PrintWrapperCommand=True    # Print wrapper name
PrintWrapperParameters=True # Print wrapper parameters
ViewVerifyedFiles=True      # Show file verification 

import os,shutil,sys,time
import xmipp
from protlib_gui import *
from protlib_base import *
from protlib_filesystem import *
from Tkinter import IntVar

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick.key, scriptname, project)
        self.Import = 'from xmipp_protocol_particle_pick import *'

    def defineSteps(self):
        self.mD = xmipp.MetaData()
        self.MicrographSelfile = os.path.join(self.self.PreprocessingRun, "micrographs.sel")
        xmipp.readMetaDataWithTwoPossibleImages(self.MicrographSelfile, self.mD)
        self.isPairList = self.mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                          not xmipp.FileName(self.MicrographSelfile).isStar1()
    
        # Create link to input micrographs
        micrographSelfile=os.path.join(self.PreprocessingRun,"micrographs.sel")
        self.Db.insertStep('createLinkToMicrographs', [micrographSelfile], MicrographSelfile=micrographSelfile, WorkingDir=self.WorkingDir)       

        # Create a dictionary of relevant micrograph information
        micrographDict={}
        containsEnable=self.mD.containsLabel(xmipp.MDL_ENABLED)
        for id in self.mD:
            micrograph = self.mD.getValue(xmipp.MDL_IMAGE,id)
            if containsEnable:
                if self.mD.getValue(xmipp.MDL_ENABLED,id)==0:
                    continue
            micrographName=micrograph.split("/")[-2]
            micrographDict[micrographName]=[micrograph]
            if self.isPairList:
                tilted=self.mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
                micrographDict[micrographName].append(tilted)
                micrographDict[micrographName].append(tilted.split("/")[-2])

        # Launch GUI
        self.Db.insertStep('launchParticlePickingGUI', passDb=True, ProjectDir=self.projectDir, WorkingDir=self.WorkingDir,
                             micrographDict=micrographDict)       

    def summary(self):
        summary = []
        #FIXME:
        return summary
        MicrographSelfile = os.path.join(self.PreprocessingRun,"micrographs.sel")
        if self.isPairList:
            summary=["Input: "+MicrographSelfile+" (Tilt pairs)"]
        else:
            summary=["Input: "+MicrographSelfile]
        
        total_manual = 0
        total_auto = 0
        N_manual = 0
        N_auto = 0
        mD = MetaData(MicrographSelfile)
        for id in mD:
             micrograph = mD.getValue(xmipp.MDL_IMAGE,id)
             manual=CountPicked(self.WorkingDir,micrograph,"Common")
             if manual>0:
                 total_manual+=manual
                 N_manual+=1
             if AutomaticPicking:
                 auto=CountPicked(self.WorkingDir,micrograph,"Common.auto")
                 if auto>0:
                     total_auto+=auto
                     N_auto+=1
        summary.append("# Manually picked: %d (from %d micrographs)" % (total_manual, N_manual))
        if AutomaticPicking:
            summary.append("# Automatically picked: %d (from %d micrographs) " % (total_auto, N_auto))
        return summary
    
    def validate(self):
        errors = []
        
        #Fixme:
        return errors

        # Check that there is a valid list of micrographs
        if not os.path.exists(self.MicrographSelfile)>0:
            errors.append("Cannot find "+self.MicrographSelfile)
        
        # Check that all micrographs exist
        NnotFound=0
        message=""
        for id in self.mD:
             micrograph = self.mD.getValue(xmipp.MDL_IMAGE,id)
             if not os.path.exists(micrograph):
                message+="  "+micrograph
                NnotFound=NnotFound+1
             if self.isPairList:
                 micrograph = self.mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
                 if not os.path.exists(micrograph):
                     message+="  "+micrograph
                     NnotFound=NnotFound+1
        
        if NnotFound>0:
            errors.append("Cannot find the following micrographs: "+message)
                
        # Check that automatic particle picking is not for tilted
        if self.isPairList and AutomaticPicking:
            errors.append("Automatic particle picking cannot be done on tilt pairs")
        
        return errors

# GUI For marking micrographs
class ProtParticlePickingGUI(BasicGUI):
    def __init__(self,log,Db,ProjectDir,WorkingDir,micrographDict):
        self.log=log
        self.Db=Db
        self.ProjectDir=ProjectDir
        self.WorkingDir=WorkingDir
        self.micrographDict=micrographDict

        # Read micrographs
        self.mD=xmipp.MetaData();
        self.MicrographSelfile=self.PreprocessingRun+"/micrographs.sel"
        xmipp.readMetaDataWithTwoPossibleImages(self.MicrographSelfile, self.mD)
        self.IsPairList = self.mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                          not xmipp.FileName(self.MicrographSelfile).isStar1()
    
    def createGUI(self):
        self.createBasicGUI()
        self.createScrollableCanvas()
 
    def fillGUI(self):
        # Init GUI variables
        self.whichmark=StringVar() # Which mark button has been pressed
        if AutomaticPicking:
            self.selectedForAutomaticPickingAuto=[] # Checkboxes for automatic particle picking
            self.selectedForAutomaticPickingMicrograph=[] # List of micrograph names associated to each button
            self.selectedForAutomaticPickingMark=[] # Mark buttons in automatic particle picking must be disabled
            titleSpan=5
        else:
            titleSpan=3

        # Window title
        self.master.title("GUI for Xmipp particle picking")
        
        # Frame header
        self.addLabel("Project directory: "+self.ProjectDir,0,0,titleSpan,fgColor="medium blue",sticky=W)
        self.addLabel("Preprocessing directory: "+self.PreprocessingRun,1,0,titleSpan,fgColor="medium blue",sticky=W)
 
        sectionFrame, sectionLabel, sectionContent = createSection(self.frame, 'Micrographs')
        
        # Add all micrographs
        containsEnable=self.mD.containsLabel(xmipp.MDL_ENABLED)
        for micrograph in self.micrographDict.keys():
            self.GuiAddMarkEntry(sectionFrame, micrograph)
        
        sectionFrame.grid(row=2, column=0)
        
        # Add blue line surrounded by two empty rows 
        self.addLabel("")
        self.addLine("medium blue",0,6)
        if AutomaticPicking:
            row=self.currentRow()
            self.addLabel("Manual",row=row,column=3)
            self.addLabel("Auto",row=row,column=4)
        else:
            self.addLabel("")

        # Add some buttons
        self.buttonrow=(self.frame.grid_size()[1]+1)
        self.addButton("Show preprocessing",self.buttonrow,0,self.ShowPreprocessing,sticky=E)
        if AutomaticPicking:
            self.addButton("Invert Selection",self.buttonrow,1,self.InvertSelection)
            nextColumn=2
        else:
            nextColumn=1
        self.addButton("Update Count",self.buttonrow,nextColumn,self.GuiUpdateCount,binding='<Control_L><U>')
        nextColumn+=1
        if AutomaticPicking:
            self.addButton("AutoSelect",self.buttonrow+1,1,self.AutomaticallyDetect)

        # Update counts
        self.GuiUpdateCount()
        
    def GuiAddMarkEntry(self, parent, micrograph):
        return
        row=parent.grid_size()[1]
        self.micrographDict[micrograph].append(row)

        # Add Micrograph name
        label=micrograph
        if self.IsPairList:
            label+=' : '+self.micrographDict[micrograph][2]
        self.addLabel(label, row, 0, sticky=E, parent=parent)
        
        # If automatic particle picking, add checkbox
        if AutomaticPicking:
            self.selectedForAutomaticPickingAuto.append(IntVar())
            self.selectedForAutomaticPickingMicrograph.append(micrograph)
            self.addCheckButton("Auto", row, 1, self.selectedForAutomaticPickingAuto[-1], 0, self.AutoSelectionChanged, "", parent)
            nextColumn=2
        else:
            nextColumn=1
        
        # Add Mark button
        radio=self.addRadioButton("Mark", row, nextColumn, self.whichmark, micrograph, self.LaunchMarkFromGUI, N, parent)
        if AutomaticPicking:
            self.selectedForAutomaticPickingMark.append(radio)

    def InvertSelection(self):
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            self.selectedForAutomaticPickingAuto[i].set(1-
                self.selectedForAutomaticPickingAuto[i].get());
        self.AutoSelectionChanged();

    def AutoSelectionChanged(self):
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            if (self.selectedForAutomaticPickingAuto[i].get()):
                self.selectedForAutomaticPickingMark[i].config(state=DISABLED)
            else:
                self.selectedForAutomaticPickingMark[i].config(state=NORMAL)

    def AutomaticallyDetect(self):
        mDauto=xmipp.MetaData();
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            if (self.selectedForAutomaticPickingAuto[i].get()):
                id=mDauto.addObject()
                mDauto.setValue(xmipp.MDL_IMAGE,self.micrographDict[self.selectedForAutomaticPickingMicrograph[i]][0],id)
        if mDauto.size()>0:
            mDauto.write(os.path.join(self.WorkingDir, "forAutomaticPicking.xmd"))

            # Create a project
            project=XmippProject(os.getcwd())
            project.load()
            
            # Get run parameters
            run=project.newOrLoadProtocol(protDict.particle_pick_auto.key, RunName+"_auto")
            run["comment"]="Automatic particle picking from ParticlePicking/"+RunName
            
            gui = ProtocolGUI()
            gui.createGUI(project, run, Toplevel())
            gui.fillGUI()
            gui.launchGUI()
            
    def GuiUpdateCount(self):
        total_manual=0
        total_auto=0
        if AutomaticPicking:
            familyColumn=3
        else:
            familyColumn=2

        # Count all micrographs
        for micrograph in self.micrographDict.keys():
            row=self.micrographDict[micrograph][-1]
            manual=CountPicked(self.WorkingDir,micrograph,"Common")
            total_manual+=manual
            self.addLabel(str(manual).zfill(5),row,familyColumn, sticky=N+S+W+E)
            if AutomaticPicking:
                auto=CountPicked(self.WorkingDir,micrograph,"Common.auto")
                total_auto+=auto
                self.addLabel(str(auto).zfill(5),row,familyColumn+1, sticky=N+S+W+E)

        # Add summaries
        self.addLabel(str(total_manual).zfill(5),self.buttonrow,familyColumn, sticky=N+S+W+E)
        if AutomaticPicking:
            self.addLabel(str(total_auto).zfill(5),self.buttonrow,familyColumn+1, sticky=N+S+W+E)

    def ShowPreprocessing(self):
        self.GuiUpdateCount()
        command='xmipp_visualize_preprocessing_micrographj -i '+self.MicrographSelfile+" &"
        os.system(command)

    def LaunchMarkFromGUI(self):
        self.GuiUpdateCount()
        micrograph=self.whichmark.get()
        if micrograph=='':
            return
        arguments=' -i '+self.micrographDict[micrograph][0]+" --outputRoot "+self.WorkingDir+"/"+micrograph
        if AutomaticPicking:
            arguments+=' --auto '+self.WorkingDir+"/Common.auto"
        elif tilted!="":
            arguments+=" --tilted "+self.micrographDict[micrograph][1]
        arguments+=" --thr "+str(NumberOfThreads)
        runJob(self.log,'xmipp_micrograph_mark',arguments,False,1,1,"",True)
    
def CountPicked(WorkingDir,micrograph,label):
    posfile=WorkingDir+"/"+str(micrograph)+'.'+label+'.pos'
    if os.path.exists(posfile):
        mD=xmipp.MetaData(posfile);
        return mD.size()
    return 0
    
# Actions
def createLinkToMicrographs(log,WorkingDir,MicrographSelfile):
    inputSelfile=os.path.join(WorkingDir,"micrographs.sel")
    if not os.path.exists(inputSelfile):
        os.system("ln -s "+os.path.relpath(os.path.abspath(MicrographSelfile),WorkingDir)+" "+inputSelfile)

# Execute protocol in the working directory
def launchParticlePickingGUI(log,Db,ProjectDir,WorkingDir,micrographDict):
    gui=ProtParticlePickingGUI(log,Db,ProjectDir,WorkingDir,micrographDict)
    gui.createGUI()
    gui.fillGUI()
    gui.launchGUI()

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
