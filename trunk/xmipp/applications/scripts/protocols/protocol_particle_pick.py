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

# Perform automatic particle picking
""" Perform automatic particle picking """
AutomaticPicking=False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/media/usbdisk/Experiments/TestProtocols'

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
import os,shutil,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path
import xmipp

# Taken from Python 2.6
import posixpath
def relpath(path, start=posixpath.curdir):
    """Return a relative version of a path"""
    if not path:
        raise ValueError("no path specified")
    start_list=posixpath.abspath(start).split(posixpath.sep)
    path_list=posixpath.abspath(path).split(posixpath.sep)
    # Work out how much of the filepath is shared by start and path.
    i=len(posixpath.commonprefix([start_list, path_list]))
    rel_list=[posixpath.pardir] * (len(start_list) - i) + path_list[i:]
    if not rel_list:
        return curdir
    return posixpath.join(*rel_list)

# Create a GUI automatically from a selfile of micrographs
class particle_pick_class:
    def saveAndCompareParameters(self, listOfParameters):
        fnOut=self.WorkingDir + "/protocolParameters.txt"
        linesNew=[];
        for prm in listOfParameters:
            eval("linesNew.append('"+prm +"='+str(self."+prm+")+'\\n')")
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

    def __init__(self,
                 WorkingDir,
                 MicrographSelfile,
                 AutomaticPicking,
                 ProjectDir,
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
        self.MicrographSelfile=MicrographSelfile
        self.PosName="Common"
        self.AutomaticPicking=AutomaticPicking
        self.ProjectDir=ProjectDir
        self.LogDir="Logs"
        self.DoParallel=DoParallel
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     self.LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Read micrographs
        self.mD=xmipp.MetaData();
        xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, self.mD)
        self.IsPairList = self.mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                          not xmipp.FileName(MicrographSelfile).isStar1()

        # Store parameters
        if os.path.exists(self.WorkingDir):
            self.saveAndCompareParameters(["MicrographSelfile",
                                           "IsPairList"]);
            
        # Make working directory if it does not exist yet
        import time
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)
            fh=open(self.WorkingDir + "/status.txt", "w")
            fh.write("Process started at " + time.asctime() + "\n")
            fh.write("Step 0: Directory created at " + time.asctime() + "\n")
            fh.close()
            os.system("ln -s "+\
                      relpath(os.path.abspath(self.MicrographSelfile),self.WorkingDir)+" "+\
                    self.WorkingDir+"/micrographs.sel")
            self.saveAndCompareParameters(["MicrographSelfile",
                                           "IsPairList"]);
        else:
            fh=open(self.WorkingDir + "/status.txt", "a")
            fh.write("Process started at " + time.asctime() + "\n")
            fh.close()

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
        			       os.path.abspath(self.WorkingDir))
    
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
        command_file = open(self.WorkingDir+"/pick.sh", "w")
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
        command_file.close();

        import tkMessageBox
        answer=tkMessageBox._show("Execute autodetection",
                                  "Use a job queueing system?",
                                  tkMessageBox.QUESTION, 
                                  tkMessageBox.YESNO)
        import launch_job
        command=""
        if self.DoParallel:
            command=launchJob("xmipp_run",
                                 "-i "+self.WorkingDir+"/pick.sh",
                                 self.log,
                                 True,
                                 self.NumberOfMpiProcesses,
                                 1,
                                 self.SystemFlavour,
                                 answer=="yes" or answer==True)
        else:
            os.system("chmod 755 "+self.WorkingDir+"/pick.sh");
            command=launchJob(self.WorkingDir+"/pick.sh",
                                 "",
                                 self.log,
                                 False,
                                 self.NumberOfMpiProcesses,
                                 1,
                                 self.SystemFlavour,
                                 answer=="yes" or answer==True)
        if (answer=="yes" or answer==True):
            import xmipp_protocol_gui
            python_file = open(self.WorkingDir+"/pick.py", "w")
            python_file.write("WorkingDir='"+self.WorkingDir+"'\n")
            python_file.write("DoParallel="+str(self.DoParallel)+"\n")
            python_file.write("NumberOfMpiProcesses="+str(self.NumberOfMpiProcesses)+"\n")
            python_file.write("# {end-of-header}\n")
            python_file.write("import os\n")
            python_file.write('os.system("'+command+'")\n');
            python_file.close();

            d = xmipp_protocol_gui.MyQueueLaunch(self.master,"python "+self.WorkingDir+"/pick.py")
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
        basedirectory,dummy=os.path.split(self.MicrographSelfile)
        directoryUname,uname=os.path.split(untilted)
        directoryTname,tname=os.path.split(tilted)
        filename=self.WorkingDir+"/"+uname+"."+self.PosName+".pos"
        command='( xmipp_micrograph_mark -i '+basedirectory+"/"+directoryUname+"/"+uname+\
                ' --tilted '+basedirectory+"/"+directoryTname+"/"+tname+\
                " --outputRoot "+self.WorkingDir+"/"+uname
        command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step M: "'+uname+' manually marked pair on `date` >> ' + \
                 self.WorkingDir + "/status.txt; fi "
        command+=") &"
        self.log.info(command)
        os.system(command)

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
    errors.append("Cannot find the following micrographs:\n")
    NnotFound=0
    for id in mD:
         micrograph = mD.getValue(xmipp.MDL_IMAGE)
         if not os.path.exists(micrograph):
            message+=micrograph+"\n"
            NnotFound=NnotFound+1
         if isPairList:
             micrograph = mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1)
             if not os.path.exists(micrograph):
                 message+=micrograph+"\n"
                 NnotFound=NnotFound+1
    
    if not NnotFound>0:
        errors.append(message)
            
    # Check that automatic particle picking is not for tilted
    if isPairList and AutomaticPicking:
        errors.append("Automatic particle picking cannot be done on tilt pairs")
        
    return errors

#		
# Main
#     
if __name__ == '__main__':
	particle_pick=particle_pick_class(WorkingDir,
                                      MicrographSelfile,
                                      AutomaticPicking,
                                      ProjectDir,
                                      DoParallel,
                                      NumberOfMpiProcesses,
                                      SystemFlavour
                                      )
