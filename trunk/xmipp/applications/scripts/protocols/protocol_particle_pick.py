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
# {dir} Working subdirectory:
""" Do not use the same directory where you executed xmipp_protocol_preprocess_micrographs.py
"""
WorkingDir = "ParticlePicking"

# {dir}{expert} Directory with the preprocessing
PreprocessingDir = "Preprocessing"

# Perform automatic particle picking
""" Perform automatic particle picking """
AutomaticPicking = False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir = "/media/usbdisk/Experiments/TestProtocols/Marcaje_automatico"

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
import os,shutil,sys,time
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

        # Setup logging
        self.log=log.init_log_system(ProjectDir,
                                     "Log",
                                     sys.argv[0],
                                     WorkingDir)
                
        # Read micrographs
        self.mD=xmipp.MetaData();
        self.MicrographSelfile=PreprocessingDir+"/micrographs.sel"
        xmipp.readMetaDataWithTwoPossibleImages(self.MicrographSelfile, self.mD)
        self.IsPairList = self.mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1) and \
                          not xmipp.FileName(self.MicrographSelfile).isStar1()

        # Store parameters
        if os.path.exists(WorkingDir):
            self.saveAndCompareParameters(["PreprocessingDir",
                                           "self.IsPairList"]);
            
        # Make working directory if it does not exist yet
        if not os.path.exists(WorkingDir):
            os.makedirs(WorkingDir)
            fh=open(WorkingDir + "/status.txt", "w")
            fh.write("Process started at " + time.asctime() + "\n")
            fh.write("Step 0: Directory created at " + time.asctime() + "\n")
            fh.close()
            os.system("ln -s "+\
                      os.path.relpath(os.path.abspath(self.MicrographSelfile),WorkingDir)+" "+\
                      WorkingDir+"/micrographs.sel")
            self.saveAndCompareParameters(["PreprocessingDir",
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
        print '*********************************************************************'
        print "* DON'T FORGET TO SAVE YOUR COORDINATES REGULARLY, AND ALWAYS BEFORE CLOSING!"
        if (self.IsPairList):
            print '* AND ALSO SAVE THE ANGLES IN THE UNTILTED MICROGRAPHS!'

    def createGUI(self):
        self.createBasicGUI()
        self.createScrollableCanvas()
 
    def addLabel(self,_text,_row=-1,_column=0,_columnspan=1,_sticky=EW,_bgColor="",_fgColor=""):
        if _fgColor=="":
            _fgColor=self.style.LabelTextColor
        if _bgColor=="":
            _bgColor=self.style.LabelBgColor
        if _row==-1:
            _row=self.frame.grid_size()[1]+1
        label=Label(self.frame, text=_text, bg=_bgColor, fg=_fgColor)
        label.grid(row=_row,column=_column,sticky=_sticky)
        return label

    def addCheckButton(self,_text,_row,_column,_default,_command,_sticky):
        controlVar = IntVar()
        controlVar.set(_default)
        check=Checkbutton(self.frame, _text, variable=controlVar,
                      command=_command, 
                      selectcolor=self.style.BooleanSelectColor)
        check.grid(row=_row, column=_column, sticky=_sticky)
        return check

    def addRadioButton(self,_text,_row,_column,_variable,_value,_command,_sticky=""):
        radio=Radiobutton(self.frame,text=_text,variable=_variable,
                          value=_value,indicatoron=0,
                          command=_command, 
                          bg=self.style.ButtonBgColor, 
                          activebackground=self.style.ButtonActiveBgColor,
                          highlightbackground=self.style.HighlightBgColor, 
                          selectcolor=self.style.ButtonActiveBgColor)
        radio.grid(row=_row, column=_column,sticky=_sticky)
        return radio

    def addButton(self,_text,_row,_column,_command,_sticky="",_binding=""):
        button = Button(self.frame, text=_text, 
                        command=_command, 
                        bg=self.style.ButtonBgColor, 
                        activebackground=self.style.ButtonActiveBgColor)
        button.grid(row=_row,column=_column,sticky=_sticky)
        if _binding!="":
            self.master.bind(_binding, _command)
        return button

    def addLine(self,_color,_column,_columnspan):
        line=Frame(self.frame, height=2, bd=1, bg=_color,relief=RIDGE)
        line.grid(row=self.frame.grid_size()[1]+1, column=_column,columnspan=_columnspan,sticky=EW)
        return line

    def fillGUI(self):
        # Init GUI variables
        self.whichmark=StringVar()
        self.row={}
        if AutomaticPicking:
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
        self.addLabel("Project directory: "+ProjectDir,0,0,titleSpan,_fgColor="medium blue",_sticky=W)
        self.addLabel("Preprocessing directory: "+PreprocessingDir,1,0,titleSpan,_fgColor="medium blue",_sticky=W)
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
        self.addButton("Show preprocessing",self.buttonrow,0,self.ShowPreprocessing,_sticky=E)
        if (AutomaticPicking):
            self.addButton("Invert Selection",self.buttonrow,1,self.InvertSelection)
            nextColumn=2
        else:
            nextColumn=1
        self.addButton("Update Total Count",self.buttonrow,nextColumn,self.GuiUpdateCount,_binding='<Control_L><U>')
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
        self.addLabel(label, row, 0, _sticky=E)
        
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
            self.addLabel(str(manual).zfill(5),row,familyColumn, _sticky=N+S+W+E)
            if AutomaticPicking:
                auto=self.CountPicked(micrograph,"Common.auto")
                total_auto+=auto
                self.addLabel(str(auto).zfill(5),row,familyColumn+1, _sticky=N+S+W+E)

        # Add summaries
        self.addLabel(str(total_manual).zfill(5),self.buttonrow,familyColumn, _sticky=N+S+W+E)
        if AutomaticPicking:
            self.addLabel(str(total_auto).zfill(5),self.buttonrow,familyColumn+1, _sticky=N+S+W+E)

    def CountPicked(self,micrograph,label):
        micrographName=micrograph.split("/")[-2]
        posfile=WorkingDir+"/"+str(micrographName)+'.'+label+'.pos'
        if os.path.exists(posfile):
            mD=xmipp.MetaData(posfile);
            return mD.size()
        return 0
    
    def ShowPreprocessing(self):
        self.GuiUpdateCount()
        command='xmipp_visualize_preprocessing_micrographj -i '+self.MicrographSelfile+" &"
        self.log.info(command)
        os.system(command)

    def LaunchSingleMark(self):
        self.GuiUpdateCount()
        micrograph=self.whichmark.get()
        if micrograph=='':
            return
        micrographName=micrograph.split("/")[-2]
        command='( xmipp_micrograph_mark -i '+micrograph+\
                  " --outputRoot "+WorkingDir+"/"+micrographName
        if (AutomaticPicking):
            command+=' --auto '+WorkingDir+"/Common.auto"
            filename=WorkingDir+"/"+micrographName+".Common.auto.pos"
            command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step L: "'+micrographName+' used for learning on `date` >> ' + \
                 WorkingDir + "/status.txt; fi "
        else:
            filename=WorkingDir+"/"+micrographName+".Common.pos"
            command+="; if [ -e " + filename + ' ]; then ' + \
                     'echo "Step M: "'+micrograph+' manually marked on `date` >> ' + \
                     WorkingDir + "/status.txt; fi "
        command+=") &"
        self.log.info(command)
        os.system(command)

    def LaunchPairMark(self):
        self.GuiUpdateCount()
        untilted=self.whichmark.get()
        tilted=self.whichtilted[untilted]
        uname=untilted.split("/")[-2]
        tname=tilted.split("/")[-2]
        filename=WorkingDir+"/"+uname+".Common.pos"
        command='( xmipp_micrograph_mark -i '+untilted+\
                ' --tilted '+tilted+\
                " --outputRoot "+WorkingDir+"/"+uname
        command+="; if [ -e " + filename + ' ]; then ' + \
                 'echo "Step M: "'+uname+' manually marked pair on `date` >> ' + \
                 WorkingDir + "/status.txt; fi "
        command+=") &"
        self.log.info(command)
        os.system(command)

# Preconditions
def checkErrors():
    errors = []
    # Check if there is workingdir
    if WorkingDir == "":
        errors.append("No working directory given")
    # Check that there is a valid list of micrographs
    MicrographSelfile=PreprocessingDir+"/micrographs.sel"
    if not os.path.exists(MicrographSelfile)>0:
        errors.append("Cannot find ")+MicrographSelfile
    
    # Check that all micrographs exist
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
