#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script (part A) for Xmipp-based manual particle picking
#
# For each micrograph in the MicrographSelfile, this program will launch
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
# {dir} Working subdirectory:
""" Do not use the same directory where you executed xmipp_protocol_preprocess_micrographs.py
"""
WorkingDir = "ParticlePicking"

# {file} Selfile with all micrographs to pick particles from:
MicrographSelfile = "Preprocessing/micrographs.sel"

# Perform automatic particle picking
""" Perform automatic particle picking """
AutomaticPicking = False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir = "/media/usbdisk/Experiments/TestProtocols"

#------------------------------------------------------------------------------------------
# {section}{condition}(AutomaticPicking=True) Parallelization issues for automatic particle picking
#------------------------------------------------------------------------------------------
# Use distributed-memory parallelization (MPI)?
""" This option provides parallelization on clusters with distributed memory architecture.
    It requires mpi to be installed.
"""
DoParallel = False

# Number of MPI processes to use:
NumberOfMpiProcesses = 1

# MPI system Flavour
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. 
"""
SystemFlavour = ""

#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
from Tkinter import *
import tkFont
import os,shutil,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path
import xmipp
from protlib_gui import *

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePicking(BasicGUI):
    def saveAndCompareParameters(self, listOfParameters):
        fnOut=WorkingDir + "/protocolParameters.txt"
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
                shutil.rmtree(WorkingDir)
                os.makedirs(WorkingDir)
        f = open(fnOut, 'w')
        f.writelines(linesNew)
        f.close()

    def __init__(self):
        self.SYSTEMSCRIPTDIR=scriptdir
        import log

        # get color definition from protocol_gui
        self.PosName="Common"

        # Setup logging
        self.log=log.init_log_system(ProjectDir,
                                     "Log",
                                     sys.argv[0],
                                     WorkingDir)
                
        # Read micrographs
        self.mD=xmipp.MetaData();
        xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, self.mD)
        self.IsPairList = self.mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                          not xmipp.FileName(MicrographSelfile).isStar1()

        # Store parameters
        if os.path.exists(WorkingDir):
            self.saveAndCompareParameters(["MicrographSelfile",
                                           "self.IsPairList"]);
            
        # Make working directory if it does not exist yet
        import time, protlib_filesystem
        if not os.path.exists(WorkingDir):
            os.makedirs(WorkingDir)
            fh=open(WorkingDir + "/status.txt", "w")
            fh.write("Process started at " + time.asctime() + "\n")
            fh.write("Step 0: Directory created at " + time.asctime() + "\n")
            fh.close()
            os.system("ln -s "+\
                      os.path.relpath(os.path.abspath(MicrographSelfile),WorkingDir)+" "+\
                      WorkingDir+"/micrographs.sel")
            self.saveAndCompareParameters(["MicrographSelfile",
                                           "IsPairList"]);
        else:
            fh=open(WorkingDir + "/status.txt", "a")
            fh.write("Process started at " + time.asctime() + "\n")
            fh.close()

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
        			       os.path.abspath(WorkingDir))
    
        # Execute protocol in the working directory
        self.print_warning()
        self.createGUI()
        self.fillGUI()
        self.launchGUI()

    def print_warning(self):
        import os, sys
        print '*********************************************************************'
        print '*  Perform manual particle picking for micrographs in: '+os.path.basename(MicrographSelfile)
        print '*'
        print '* DONT FORGET TO SAVE YOUR COORDINATES REGULARLY, AND ALWAYS BEFORE CLOSING!'
        if (self.IsPairList):
            print '* AND ALSO SAVE THE ANGLES IN THE UNTILTED MICROGRAPHS!'
        print '*'

    def createGUI(self):
        self.createBasicGUI()
        self.createScrollableCanvas()
        self.total_count=0
        self.total_count_auto=0
        self.whichmark=StringVar()
        self.whichtilted={}
        self.row={}
        self.LabelBackgroundColour=self.style.LabelBgColor
        self.ButtonBackgroundColour=self.style.ButtonBgColor
        self.ButtonActiveBackgroundColour=self.style.ButtonActiveBgColor
        self.HighlightBackgroundColour=self.style.HighlightBgColor
        self.BooleanSelectColour=self.style.BooleanSelectColor
 
    def fillGUI(self):
        import os

        # Script title
        self.master.title("GUI for Xmipp particle picking")
        headertext='GUI for Xmipp particle picking\n'
        headertext+="Executed in directory: "+str(os.getcwd())
        l1=Label(self.frame, text=headertext, fg="medium blue")
        if (AutomaticPicking):
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
        directory,dummy=os.path.split(MicrographSelfile)
        containsEnable=self.mD.containsLabel(xmipp.MDL_ENABLED)
        if not self.IsPairList:
            for id in self.mD:
                micrograph=self.mD.getValue(xmipp.MDL_IMAGE,id)
                if containsEnable:
                    if self.mD.getValue(xmipp.MDL_ENABLED,id)==0:
                        continue
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
                micrograph=self.mD.getValue(xmipp.MDL_IMAGE,id)
                tilted=self.mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
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

        if (AutomaticPicking):
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
        
        if (AutomaticPicking):
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

        label=micrograph.split("/")[-2]
        l=Label(self.frame, text=label, 
                bg=self.LabelBackgroundColour)
        l.grid(row=row, column=0, sticky=E)

        if (AutomaticPicking):
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

        if (AutomaticPicking):
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
        command_file = open(WorkingDir+"/pick.sh", "w")
        directoryPreprocessing,dummy=os.path.split(MicrographSelfile)
        for i in range(0,len(self.selectedForAutomaticPickingAuto)):
            if (self.selectedForAutomaticPickingAuto[i].get()):
               directory,micrograph=os.path.split(
                  self.selectedForAutomaticPickingName[i])
               filename=WorkingDir+"/"+micrograph+"."+self.PosName+".auto.pos"
               command_file.write(
                  "( xmipp_micrograph_mark -i "+directoryPreprocessing+"/"+\
                  self.selectedForAutomaticPickingName[i]+\
                  " --auto "+WorkingDir+"/"+self.PosName+".auto --autoSelect "+\
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

    def CountPicked(self,micrograph,label):
        directory,fnMicrograph=os.path.split(micrograph);
        posfile=WorkingDir+"/"+str(fnMicrograph)+'.'+label+'.pos'
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
            if (AutomaticPicking):
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
        basedirectory,dummy=os.path.split(MicrographSelfile)
        micrograph=name.split("/")[-2]
        command='( xmipp_micrograph_mark -i '+name+\
                  " --outputRoot "+WorkingDir+"/"+micrograph
        if (AutomaticPicking):
            command+=' --auto '+WorkingDir+"/"+self.PosName+'.auto'
            filename=WorkingDir+"/"+micrograph+"."+self.PosName+".auto.pos"
            command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step L: "'+micrograph+' used for learning on `date` >> ' + \
                 WorkingDir + "/status.txt; fi "
        filename=WorkingDir+"/"+micrograph+"."+self.PosName+".pos"
        command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step M: "'+micrograph+' manually marked on `date` >> ' + \
                 WorkingDir + "/status.txt; fi "
        command+=") &"
        self.log.info(command)
        os.system(command)

    def LaunchPairMark(self):
        import launch_job
        self.GuiUpdateCount()
        print "* Marking pair ... "
        
        untilted=self.whichmark.get()
        tilted=self.whichtilted[self.whichmark.get()]
        uname=untilted.split("/")[-2]
        tname=tilted.split("/")[-2]
        filename=WorkingDir+"/"+uname+"."+self.PosName+".pos"
        command='( xmipp_micrograph_mark -i '+untilted+\
                ' --tilted '+tilted+\
                " --outputRoot "+WorkingDir+"/"+uname
        command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step M: "'+uname+' manually marked pair on `date` >> ' + \
                 WorkingDir + "/status.txt; fi "
        command+=") &"
        self.log.info(command)
        os.system(command)

    def GuiUpdateCount(self):
        print "* Updating count ..."
        total=self.CountAll()
        label=str(total).zfill(5)
        l=Label(self.frame, text=label, 
                bg=self.LabelBackgroundColour)
        if (AutomaticPicking):
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
def checkErrors():
    errors = []
    # Check if there is workingdir
    if WorkingDir == "":
        errors.append("No working directory given")
    # Check that there is a valid list of micrographs
    if not os.path.exists(MicrographSelfile)>0:
        errors.append("Cannot find ")+MicrographSelfile
    # Check that all micrographs exist
    import xmipp
    mD = xmipp.MetaData()
    xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, mD)
    isPairList = mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                 not xmipp.FileName(MicrographSelfile).isStar1()
    NnotFound=0
    message=[]
    for id in mD:
         micrograph = mD.getValue(xmipp.MDL_IMAGE,id)
         if not os.path.exists(micrograph):
            message.append("  "+micrograph)
            NnotFound=NnotFound+1
         if isPairList:
             micrograph = mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1,id)
             if not os.path.exists(micrograph):
                 message.append("  "+micrograph)
                 NnotFound=NnotFound+1
    
    if NnotFound>0:
        errors.append("Cannot find the following micrographs:")
        errors.append(message)
            
    # Check that automatic particle picking is not for tilted
    if isPairList and AutomaticPicking:
        errors.append("Automatic particle picking cannot be done on tilt pairs")
        
    return errors

#		
# Main
#     
if __name__ == '__main__':
	particle_pick=ProtParticlePicking()
