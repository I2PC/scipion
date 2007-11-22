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
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {dir} Working subdirectory:
""" Use the same directory where you executed protocol_preprocess_micrographs.py
"""
WorkingDir="Preprocessing"
# {file} Selfile with all micrographs to pick particles from:
MicrographSelfile='Preprocessing/all_micrographs.sel'
# Is this selfile a list of untilted-tilted pairs?
""" True for RCT-processing. In that case, provide a 3-column selfile as follows:
    untilted_pair1.raw tilted_pair1.raw 1
    untilted_pair2.raw tilted_pair2.raw 1
    etc...
    Where 1 in the third column means active pair, and -1 means inactive pair
    Use relative paths from the Preprocessing directory!
"""
IsPairList=False
# Name of the position files (or family name)
""" This is specified inside the micrograph_mark program (Common by default)
"""
PosName='Common'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/home/scheres/work/protocols'
# {expert} Directory name for logfiles:
LogDir='Logs'
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
from Tkinter import *
# Create a GUI automatically from a selfile of micrographs
class particle_pick_class:

    def __init__(self,
                 WorkingDir,
                 MicrographSelfile,
                 IsPairList,
                 PosName,
                 ProjectDir,
                 LogDir):

        import os,sys
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        self.SYSTEMSCRIPTDIR=scriptdir
        import log

        self.WorkingDir=WorkingDir
        self.MicrographSelfile=os.path.abspath(MicrographSelfile)
        self.IsPairList=IsPairList
        self.PosName=PosName
        self.ProjectDir=ProjectDir
        self.LogDir=LogDir

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Make working directory if it does not exist yet
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
    
        # Execute protocol in the working directory
        os.chdir(self.WorkingDir)
        
        # Execute protocol in the working directory
        self.print_warning()
        self.sellines=self.ReadSelfile()
        self.MakeGui()

        # Return to parent dir
        os.chdir(os.pardir)

    def MakeGui(self):
        import xmipp_protocol_gui

        self.master=Tk()
        self.total_count=0
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

    def ReadSelfile(self):
        import os
        newlines=[]
        if not os.path.exists(self.MicrographSelfile):
            message='Error: '+self.MicrographSelfile+' does not exist'
            print '*',message
            self.log.error(message)
            sys.exit()
        fh=open(self.MicrographSelfile,'r')
        lines=fh.readlines()
        fh.close()
        words=lines[0].split()
        if ((not self.IsPairList) and (not len(words)==2)):
            message='Error: '+self.MicrographSelfile+' is not a valid selection file'
            print '*',message
            self.log.error(message)
            sys.exit()
        if ((self.IsPairList) and (not len(words)==3)):
            message='Error: '+self.MicrographSelfile+' is not a valid pairlist file'
            print '*',message
            self.log.error(message)
            sys.exit()
        for line in lines:
            words=line.split()
            newlines.append(words)
        return newlines
 
    def FillMarkGui(self):
        import os

        # Script title
        self.master.title("GUI for particle picking")
        headertext='GUI for Xmipp particle picking\n'
        headertext+="Executed in directory: "+str(os.getcwd())
        l1=Label(self.frame, text=headertext, fg="medium blue")
        l1.grid(row=0, column=0,columnspan=3,sticky=EW)
 
        total=0
        if not self.IsPairList:
            for micrograph,state in self.sellines:
                if (state.find('-1') == -1):
                    c=self.CountPicked(micrograph)
                    total+=c
                    row=self.GuiAddSingleMarkEntry(micrograph,c)
                    self.row[micrograph]=row
        else:
            for micrograph,tilted,state in self.sellines:
                if (state.find('-1') == -1):
                    c=self.CountPicked(micrograph)
                    total+=c
                    self.whichtilted[micrograph]=tilted
                    row=self.GuiAddPairMarkEntry(micrograph,c)
                    self.row[micrograph]=row

        row=(self.frame.grid_size()[1]+1)
        Label(self.frame,text="").grid(row=row)
        l2=Frame(self.frame, height=2, bd=1, bg="medium blue",relief=RIDGE)
        l2.grid(row=row+1, column=0,columnspan=6,sticky=EW)
        Label(self.frame,text="").grid(row=row+2)
        self.buttonrow=(self.frame.grid_size()[1]+1)
        b = Button(self.frame, text="Close", command=self.GuiClose,underline=0)
        b.grid(row=self.buttonrow,column=0,sticky=W)
        self.master.bind('<Control_L><c>', self.GuiClose)
        b = Button(self.frame, text="Update Total Count:", command=self.GuiUpdateCount)
        b.grid(row=self.buttonrow,column=1)
        label=str(total).zfill(5)
        l = Label(self.frame, text=label)
        l.grid(row=self.buttonrow,column=2)
    
    def GuiAddSingleMarkEntry(self,micrograph,count):
        import os
        row=self.frame.grid_size()[1]
        label=os.path.basename(micrograph)
        l=Label(self.frame, text=label)
        l.grid(row=row, column=0, sticky=E)
        label=str(count).zfill(5)
        l=Label(self.frame, text=label)
        l.grid(row=row, column=2)
        r=Radiobutton(self.frame,text="Mark!",variable=self.whichmark,
                           value=micrograph,indicatoron=0,command=self.LaunchSingleMark)
        r.grid(row=row, column=1,sticky=N)
        return row

    def GuiAddPairMarkEntry(self,micrograph,count):
        import os
        row=self.frame.grid_size()[1]
        label=os.path.basename(micrograph)+' : '+os.path.basename(self.whichtilted[micrograph])
        l=Label(self.frame, text=label)
        l.grid(row=row, column=0, sticky=E)
        label=str(count).zfill(5)
        l=Label(self.frame, text=label)
        l.grid(row=row, column=2)
        r=Radiobutton(self.frame,text="Mark!",variable=self.whichmark,
                           value=micrograph, indicatoron=0, command=self.LaunchPairMark)
        r.grid(row=row, column=1,sticky=N)
        return row

    def CountPicked(self,micrograph):
        import os
        posfile=str(micrograph)+'.'+str(self.PosName)+'.pos'
        if os.path.exists(posfile):
            fh=open(posfile,'r')
            lines=fh.readlines()
            fh.close()
            picked=len(lines)-1
            if picked>0:
                return picked
        return 0
    
    def CountAll(self):
        total=0
        for mic,row in self.row.items():
            c=self.CountPicked(mic)
            total=total+c
            label=str(c).zfill(5)
            l=Label(self.frame, text=label)
            l.grid(row=row, column=2)
        return total
        
    def LaunchSingleMark(self):
        import os
        self.GuiUpdateCount()
        print "* Marking... "
        self.perform_picking(self.whichmark.get())

    def LaunchPairMark(self):
        import os
        self.GuiUpdateCount()
        print "* Marking... "
        untilted=self.whichmark.get()
        self.perform_picking_pair(self.whichmark.get(),self.whichtilted[self.whichmark.get()])

    def GuiUpdateCount(self):
        print "* Updating count..."
        total=self.CountAll()
        label=str(total).zfill(5)
        l=Label(self.frame, text=label)
        l.grid(row=self.buttonrow, column=2)

    def GuiClose(self):
        import sys
        self.master.quit()
        self.master.destroy()
        sys.exit(0)
        
    def print_warning(self):
        import os
        print '*********************************************************************'
        print '*  Perform manual particle picking for micrographs in: '+os.path.basename(self.MicrographSelfile)
        print '*'
        print '* DONT FORGET TO SAVE YOUR COORDINATES REGULARLY, AND ALWAYS BEFORE CLOSING!'
        if (self.IsPairList):
            print '* AND ALSO SAVE THE ANGLES IN THE UNTILTED MICROGRAPHS!'
        print '*'

    def perform_picking(self,name):
        import os
        directory,micrograph=os.path.split(name)
        if (len(directory)>0):
            os.chdir(directory)
        command='xmipp_micrograph_mark -i '+micrograph+' &'
        print '* ',command
        self.log.info(command)
        os.system(command)
        os.chdir(os.pardir)

    def perform_picking_pair(self,untilted,tilted):
        import os
        directory,uname=os.path.split(untilted)
        if (len(directory)>0):
            os.chdir(directory)
        tname='../'+tilted
        command='xmipp_micrograph_mark -i '+uname+' -tilted '+tname+' &'
        print '* ',command
        self.log.info(command)
        os.system(command)
        os.chdir(os.pardir)


    def close(self):
        message=" Exiting ... "
        print '* ',message
        print '*********************************************************************'
        self.log.info(message)
#		
# Main
#     
if __name__ == '__main__':

   	# create preprocess_A_class object

	particle_pick=particle_pick_class(WorkingDir,
                                          MicrographSelfile,
                                          IsPairList,
                                          PosName,
                                          ProjectDir,
                                          LogDir)

	# close 
	particle_pick.close()
