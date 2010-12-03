#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script (part A) for Xmipp-based manual particle picking
#
# For each micrograph in the MicrographSelfile, this program will launch
# the xmipp_mark program 
# A graphical interface exists to identify micrographs that have been finished
#
# Example use:
# python xmipp_particle_pick.py &
#
# Author: Sjors Scheres, March 2007
# Author: Carlos Oscar Sorzano, November, 2010
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {dir} Working subdirectory:
""" Do not use the same directory where you executed xmipp_protocol_preprocess_micrographs.py
"""
WorkingDir='ParticlePicking'

# {file} Selfile with all micrographs to pick particles from:
MicrographSelfile='Preprocessing/all_micrographs.sel'

# Is this selfile a list of untilted-tilted pairs?
""" True for RCT-processing. In that case, provide a 2-column selfile as follows:
    untilted_pair1.raw tilted_pair1.raw
    untilted_pair2.raw tilted_pair2.raw
    etc...
    Use relative paths from the Preprocessing directory!
"""
IsPairList=False

# Name of the position files (or family name)
""" This is specified inside the micrograph_mark program (Common by default)
"""
PosName='Common'

# Perform automatic particle picking
""" Perform automatic particle picking """
AutomaticPicking=False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/media/usbdisk/Experiments/TestProtocols'

# {expert} Directory name for logfiles:
LogDir='Logs'

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues for automatic particle picking
#------------------------------------------------------------------------------------------------
# Use distributed-memory parallelization (MPI)?
""" This option provides parallelization on clusters with distributed memory architecture.
    It requires mpi to be installed.
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
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
from Tkinter import *
import os,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path
import xmipp

# Create a GUI automatically from a selfile of micrographs
class particle_pick_class:
    def __init__(self,
                 WorkingDir,
                 MicrographSelfile,
                 IsPairList,
                 PosName,
                 AutomaticPicking,
                 ProjectDir,
                 LogDir,
                 DoParallel,
		         NumberOfMpiProcesses,
                 SystemFlavour
                 ):

        self.SYSTEMSCRIPTDIR=scriptdir
        import log

        # get color definition from protocol_gui
        import xmipp_protocol_gui
        self.LabelBackgroundColour=xmipp_protocol_gui.LabelBackgroundColour
        self.ButtonBackgroundColour=xmipp_protocol_gui.ButtonBackgroundColour
        self.ButtonActiveBackgroundColour=xmipp_protocol_gui.ButtonActiveBackgroundColour
        self.HighlightBackgroundColour=xmipp_protocol_gui.HighlightBackgroundColour
        self.BooleanSelectColour=xmipp_protocol_gui.BooleanSelectColour

        self.WorkingDir=WorkingDir
        self.MicrographSelfile=os.path.abspath(MicrographSelfile)
        self.IsPairList=IsPairList
        self.PosName=PosName
        self.AutomaticPicking=AutomaticPicking
        self.ProjectDir=ProjectDir
        self.LogDir=LogDir
        self.DoParallel=DoParallel
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour

        # Delete working directory if exists, make a new one
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Make working directory if it does not exist yet
        import time
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)
            fh=open(self.WorkingDir + "/status.txt", "w")
            fh.write("Process started at " + time.asctime() + "\n")
            fh.write("Step 0: Directory created at " + time.asctime() + "\n")
            fh.close()
        else:
            fh=open(self.WorkingDir + "/status.txt", "a")
            fh.write("Process started at " + time.asctime() + "\n")
            fh.close()

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
        			       os.path.abspath(self.WorkingDir))
    
        # Read micrographs
        self.mD=xmipp.MetaData();
        self.mD.read(MicrographSelfile)

        # Execute protocol in the working directory
        self.print_warning()
        self.MakeGui()

    def print_warning(self):
        import os, sys
        print '*********************************************************************'
        print '*  Perform manual particle picking for micrographs in: '+os.path.basename(self.MicrographSelfile)
        print '*'
        print '* DONT FORGET TO SAVE YOUR COORDINATES REGULARLY, AND ALWAYS BEFORE CLOSING!'
        if (self.IsPairList):
            print '* AND ALSO SAVE THE ANGLES IN THE UNTILTED MICROGRAPHS!'
        print '*'

    def MakeGui(self):
        import xmipp_protocol_gui

        self.master=Tk()
        self.total_count=0
        self.total_count_auto=0
        self.whichmark=StringVar()
        self.whichtilted={}
        self.row={}

        # Create the Canvas with Scrollbars
        self.canvas,self.frame=xmipp_protocol_gui.PrepareCanvas(self.master)

        # Fill the GUI
        self.FillMarkGui()

        # Launch the window
        xmipp_protocol_gui.LaunchCanvas(self.master,self.canvas,self.frame)
        xmipp_protocol_gui.GuiResize(self.master,self.frame)

        # Enter main loop
        self.master.mainloop()
 
    def FillMarkGui(self):
        import os

        # Script title
        self.master.title("GUI for particle picking")
        headertext='GUI for Xmipp particle picking\n'
        headertext+="Executed in directory: "+str(os.getcwd())
        l1=Label(self.frame, text=headertext, fg="medium blue")
        if (self.AutomaticPicking):
            l1.grid(row=0, column=0,columnspan=5,sticky=EW)
            Label(self.frame, text="Manual").grid(row=1,column=3)
            Label(self.frame, text="Auto").grid(row=1,column=4)
        else:
            l1.grid(row=0, column=0,columnspan=3,sticky=EW)
 
        total=0
        total_auto=0
        self.selectedForAutomaticPickingTrain=[]
        self.selectedForAutomaticPickingAuto=[]
        self.selectedForAutomaticPickingName=[]
        self.selectedForAutomaticPickingMark=[]
        directory,dummy=os.path.split(self.MicrographSelfile)
        if not self.IsPairList:
            for id in self.mD:
                micrograph=self.mD.getValue(xmipp.MDL_IMAGE)
                if not os.path.exists(directory+"/"+micrograph):
                    c=0
                    cauto=0
                    row=self.GuiAddSingleMarkEntry(micrograph,c,cauto)
                    self.row[micrograph]=row
                else:
                    c=self.CountPicked(micrograph,str(self.PosName))
                    cauto=self.CountPicked(micrograph,str(self.PosName+'.auto'))
                    total+=c
                    total_auto+=cauto
                    row=self.GuiAddSingleMarkEntry(micrograph,c,cauto)
                    self.row[micrograph]=row
        else:
            for id in self.mD:
                micrograph=self.mD.getValue(xmipp.MDL_IMAGE)
                tilted=self.mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1)
                c=self.CountPicked(micrograph,str(self.PosName))
                total+=c
                self.whichtilted[micrograph]=tilted
                row=self.GuiAddPairMarkEntry(micrograph,c)
                self.row[micrograph]=row

        row=(self.frame.grid_size()[1]+1)
        Label(self.frame,text="", 
              bg=self.LabelBackgroundColour).grid(row=row)
        l2=Frame(self.frame, height=2, bd=1, bg="medium blue",relief=RIDGE)
        l2.grid(row=row+1, column=0,columnspan=6,sticky=EW)
        Label(self.frame,text="", 
              bg=self.LabelBackgroundColour).grid(row=row+2)
        self.buttonrow=(self.frame.grid_size()[1]+1)

        b = Button(self.frame, text="Close", 
                   command=self.GuiClose,underline=0, 
                   bg=self.ButtonBackgroundColour, 
                   activebackground=self.ButtonActiveBackgroundColour)
        b.grid(row=self.buttonrow,column=0,sticky=W)
        self.master.bind('<Control_L><c>', self.GuiClose)

        if (self.AutomaticPicking):
            b = Button(self.frame, text="Invert Selection", 
                       command=self.InvertSelection, 
                       bg=self.ButtonBackgroundColour, 
                       activebackground=self.ButtonActiveBackgroundColour)
            b.grid(row=self.buttonrow,column=1,sticky=N+S+W+E)
            nextColumn=2
        else:
            nextColumn=1

        b = Button(self.frame, text="Update Total Count:", 
                   command=self.GuiUpdateCount, underline=0,
                   bg=self.ButtonBackgroundColour, 
                   activebackground=self.ButtonActiveBackgroundColour)
        b.grid(row=self.buttonrow,column=nextColumn)
        self.master.bind('<Control_L><U>', self.GuiUpdateCount)
        nextColumn+=1

        label=str(total).zfill(5)
        l = Label(self.frame, text=label, 
                  bg=self.LabelBackgroundColour)
        l.grid(row=self.buttonrow,column=nextColumn)
        nextColumn+=1
        
        if (self.AutomaticPicking):
            label=str(total_auto).zfill(5)
            l = Label(self.frame, text=label, 
                      bg=self.LabelBackgroundColour)
            l.grid(row=self.buttonrow,column=nextColumn)

            b = Button(self.frame, text="AutoSelect",
                       command=self.AutomaticallyDetect,underline=0, 
                       bg=self.ButtonBackgroundColour, 
                       activebackground=self.ButtonActiveBackgroundColour)
            b.grid(row=self.buttonrow+1,column=1,sticky=N+S+W+E)
        
    def GuiAddSingleMarkEntry(self,micrograph,count,count_auto):
        row=self.frame.grid_size()[1]

        label=os.path.basename(micrograph)
        l=Label(self.frame, text=label, 
                bg=self.LabelBackgroundColour)
        l.grid(row=row, column=0, sticky=E)

        if (self.AutomaticPicking):
            controlTrain = IntVar()
            c=Checkbutton(self.frame, text="Train", variable=controlTrain,
                          command=self.TrainSelectionChanged, 
                          selectcolor=self.BooleanSelectColour)
            c.grid(row=row, column=1, sticky=N+W)
            controlTrain.set(1)
            self.selectedForAutomaticPickingTrain.append(controlTrain)

            controlAuto = IntVar()
            c=Checkbutton(self.frame, text="Auto", variable=controlAuto,
                          command=self.AutoSelectionChanged, 
                          selectcolor=self.BooleanSelectColour)
            c.grid(row=row, column=1, sticky=S+E)
            controlAuto.set(0)
            self.selectedForAutomaticPickingAuto.append(controlAuto)

            self.selectedForAutomaticPickingName.append(micrograph)
            nextColumn=2
        else:
            nextColumn=1

        r=Radiobutton(self.frame,text="Mark",variable=self.whichmark,
                      value=micrograph,indicatoron=0,
                      command=self.LaunchSingleMark, 
                      bg=self.ButtonBackgroundColour, 
                      activebackground=self.ButtonActiveBackgroundColour,
                      highlightbackground=self.HighlightBackgroundColour, 
                      selectcolor=self.ButtonActiveBackgroundColour)
        r.grid(row=row, column=nextColumn,sticky=N)
        self.selectedForAutomaticPickingMark.append(r)
        nextColumn+=1

        label=str(count).zfill(5)
        l=Label(self.frame, text=label, 
                bg=self.LabelBackgroundColour)
        l.grid(row=row, column=nextColumn, sticky=N+S+W+E)
        nextColumn+=1

        if (self.AutomaticPicking):
            label=str(count_auto).zfill(5)
            l=Label(self.frame, text=label, 
                    bg=self.LabelBackgroundColour)
            l.grid(row=row, column=nextColumn, sticky=N+S+W+E)
        return row

    def GuiAddPairMarkEntry(self,micrograph,count):
        row=self.frame.grid_size()[1]
        label=os.path.basename(micrograph)+' : '+os.path.basename(self.whichtilted[micrograph])
        l=Label(self.frame, text=label, 
                bg=self.LabelBackgroundColour)
        l.grid(row=row, column=0, sticky=E)
        label=str(count).zfill(5)
        l=Label(self.frame, text=label, 
                bg=self.LabelBackgroundColour)
        l.grid(row=row, column=2)
        r=Radiobutton(self.frame,text="Mark",variable=self.whichmark,
                           value=micrograph, indicatoron=0,
                           command=self.LaunchPairMark)
        r.grid(row=row, column=1,sticky=N)
        return row

    def InvertSelection(self):
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            self.selectedForAutomaticPickingAuto[i].set(1-
                self.selectedForAutomaticPickingAuto[i].get());
        self.AutoSelectionChanged();

    def TrainSelectionChanged(self):
        for i in range(0,len(self.selectedForAutomaticPickingTrain)):
            if (self.selectedForAutomaticPickingTrain[i].get()):
                self.selectedForAutomaticPickingMark[i].config(state=NORMAL)
                self.selectedForAutomaticPickingAuto[i].set(0)
            else:
                self.selectedForAutomaticPickingMark[i].config(state=DISABLED)
                self.selectedForAutomaticPickingAuto[i].set(1)

    def AutoSelectionChanged(self):
        for i in range(0,len(self.selectedForAutomaticPickingTrain)):
            if (self.selectedForAutomaticPickingAuto[i].get()):
                self.selectedForAutomaticPickingMark[i].config(state=DISABLED)
                self.selectedForAutomaticPickingTrain[i].set(0)
            else:
                self.selectedForAutomaticPickingMark[i].config(state=NORMAL)
                self.selectedForAutomaticPickingTrain[i].set(1)

    def AutomaticallyDetect(self):
        command_file = open("pick.sh", "w")
        directoryPreprocessing,dummy=os.path.split(self.MicrographSelfile)
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            if (self.selectedForAutomaticPickingAuto[i].get()):
               directory,micrograph=os.path.split(
                  self.selectedForAutomaticPickingName[i])
               filename=self.WorkingDir+"/"+micrograph+"."+self.PosName+".auto.pos"
               command_file.write(
                  "( xmipp_micrograph_mark -i "+directoryPreprocessing+"/"+\
                  self.selectedForAutomaticPickingName[i]+\
                  " --auto "+self.WorkingDir+"/"+self.PosName+".auto --autoSelect "+\
                  " --outputRoot "+self.WorkingDir+"/"+micrograph+\
                  "; if [ -e " + filename + ' ]; then ' + \
                  'echo "Step A: "'+micrograph+' automatically marked on `date` >> ' + \
                  self.WorkingDir + "/status.txt; fi )\n")
                 
        command_file.write("MPI_Barrier\n")
        command_file.write('echo "Step F: " finished marking on `date` >> ' + \
                  self.WorkingDir + "/status.txt \n")
        command_file.write("rm pick.sh pick.py\n")
        command_file.close();

        import tkMessageBox
        answer=tkMessageBox._show("Execute autodetection",
                                  "Use a job queueing system?",
                                  tkMessageBox.QUESTION, 
                                  tkMessageBox.YESNOCANCEL)
        import launch_job
        command=""
        if self.DoParallel:
            command=launch_job.launch_job("xmipp_run",
                                 "-i pick.sh",
                                 self.log,
                                 True,
                                 self.NumberOfMpiProcesses,
                                 1,
                                 self.SystemFlavour,
                                 answer=="yes" or answer==True)
        else:
            os.system("chmod 755 pick.sh");
            command=launch_job.launch_job("./pick.sh",
                                 "",
                                 self.log,
                                 False,
                                 self.NumberOfMpiProcesses,
                                 1,
                                 self.SystemFlavour,
                                 answer=="yes" or answer==True)
        if (answer=="yes" or answer==True):
            import xmipp_protocol_gui
            python_file = open("pick.py", "w")
            python_file.write("WorkingDir='"+self.WorkingDir+"'\n")
            python_file.write("DoParallel="+str(self.DoParallel)+"\n")
            python_file.write("NumberOfMpiProcesses="+str(self.NumberOfMpiProcesses)+"\n")
            python_file.write("# {end-of-header}\n")
            python_file.write("import os\n")
            python_file.write('os.system("'+command+'")\n');
            python_file.close();

            d = xmipp_protocol_gui.MyQueueLaunch(self.master,"python pick.py")
        else:
            os.system(command);

    def CountPicked(self,micrograph,label):
        directory,fnMicrograph=os.path.split(micrograph);
        posfile=self.WorkingDir+"/"+str(fnMicrograph)+'.'+label+'.pos'
        if os.path.exists(posfile):
            import xmipp
            mD=xmipp.MetaData(posfile);
            return mD.size()
        return 0
    
    def CountAll(self, family=''):
        total=0
        for mic,row in self.row.items():
            if family=='':
                c=self.CountPicked(mic,str(self.PosName))
            else:
                c=self.CountPicked(mic,family)
            total=total+c
            label=str(c).zfill(5)
            l=Label(self.frame, text=label, 
                    bg=self.LabelBackgroundColour)
            if (self.AutomaticPicking):
                if family=='':
                    l.grid(row=row, column=3)
                else:
                    l.grid(row=row, column=4)
            else:
                l.grid(row=row, column=2)
        return total
        
    def LaunchSingleMark(self):
        import launch_job

        self.GuiUpdateCount()
        print "* Marking ... "
        name=self.whichmark.get()
        if name=='':
            return
        basedirectory,dummy=os.path.split(self.MicrographSelfile)
        directory,micrograph=os.path.split(name)
        command='( xmipp_micrograph_mark -i '+basedirectory+"/"+directory+"/"+micrograph+\
                  " --outputRoot "+self.WorkingDir+"/"+micrograph
        if (self.AutomaticPicking):
            command+=' --auto '+self.WorkingDir+"/"+self.PosName+'.auto'
            filename=self.WorkingDir+"/"+micrograph+"."+self.PosName+".auto.pos"
            command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step L: "'+micrograph+' used for learning on `date` >> ' + \
                 self.WorkingDir + "/status.txt; fi "
        filename=self.WorkingDir+"/"+micrograph+"."+self.PosName+".pos"
        command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step M: "'+micrograph+' manually marked on `date` >> ' + \
                 self.WorkingDir + "/status.txt; fi "
        command+=") &"
        self.log.info(command)
        os.system(command)

    def LaunchPairMark(self):
        import launch_job
        self.GuiUpdateCount()
        print "* Marking pair ... "
        
        untilted=self.whichmark.get()
        tilted=self.whichtilted[self.whichmark.get()]
        directory,uname=os.path.split(untilted)
        if (len(directory)>0):
            os.chdir(directory)
        tname='../'+tilted
        command=' -i '+uname+' --tilted '+tname + ' &'
        launch_job.launch_job("xmipp_micrograph_mark",
                              command,
                              self.log,
                              False,1,1,'')
        os.chdir(os.pardir)

    def GuiUpdateCount(self):
        print "* Updating count ..."
        total=self.CountAll()
        label=str(total).zfill(5)
        l=Label(self.frame, text=label, 
                bg=self.LabelBackgroundColour)
        if (self.AutomaticPicking):
            l.grid(row=self.buttonrow, column=3)
            totalAuto=self.CountAll(self.PosName+'.auto')
            label=str(totalAuto).zfill(5)
            l=Label(self.frame, text=label, 
                    bg=self.LabelBackgroundColour)
            l.grid(row=self.buttonrow, column=4)
        else:
            l.grid(row=self.buttonrow, column=2)

    def GuiClose(self):
        import sys
        self.master.quit()
        self.master.destroy()
        sys.exit(0)
        
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
    
    # Check that automatic particle picking is not for tilted
    if IsPairList and AutomaticPicking:
        message="No working directory given"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
        
    # Check that there is a valid list of micrographs
    if not os.path.exists(MicrographSelfile)>0:
        message="Cannot find "+MicrographSelfile
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    # Check that all micrographs exist
    import xmipp
    mD = xmipp.MetaData()
    mD.read(MicrographSelfile)
    message="Cannot find the following micrographs:\n"
    NnotFound=0
    for id in mD:
         micrograph = mD.getValue(xmipp.MDL_IMAGE)
         if not os.path.exists(micrograph):
            message+=micrograph+"\n"
            NnotFound=NnotFound+1
    
    if not NnotFound>0:
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
            
    return retval

#		
# Main
#     
if __name__ == '__main__':
	particle_pick=particle_pick_class(WorkingDir,
                                      MicrographSelfile,
                                      IsPairList,
                                      PosName,
                                      AutomaticPicking,
                                      ProjectDir,
                                      LogDir,
                                      DoParallel,
                                      NumberOfMpiProcesses,
                                      SystemFlavour
                                      )
