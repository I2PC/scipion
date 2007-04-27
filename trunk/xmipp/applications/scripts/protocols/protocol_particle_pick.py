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
# {file} Selfile with all micrographs to pick particles from:
MicrographSelfile="/home2/bioinfo/scheres/work/protocols/Preprocessing/all_micrographs.sel"
# Is this selfile a list of untilted-tilted pairs?
""" True for RCT-processing. In that case, provide a 3-column selfile as follows:
    untilted_pair1.raw tilted_pair1.raw 1
    untilted_pair2.raw tilted_pair2.raw 1
    etc...
    Where 1 in the third column means active pair, and -1 means inactive pair
"""
IsPairList=False
# Name of the position files (or family name)
""" This is specified inside the micrograph_mark program (raw.Common.pos by default)
"""
PosName="raw.Common.pos"
# Use GUI to display list of finished micrographs?
DoUseGui=True
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
LogDir="Logs"
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
                 MicrographSelfile,
                 IsPairList,
                 PosName,
                 DoUseGui,
                 ProjectDir,
                 LogDir):

        import os,sys
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        self.SYSTEMSCRIPTDIR=scriptdir
        import log

        self.MicrographSelfile=MicrographSelfile
        self.IsPairList=IsPairList
        self.PosName=PosName
        self.ProjectDir=ProjectDir
        self.LogDir=LogDir
        self.DoUseGui=DoUseGui

        self.log=log.init_log_system(self.ProjectDir,
                                     self.LogDir,
                                     sys.argv[0],
                                     '.')

        # Execute protocol:
        self.print_warning()
        self.sellines=self.ReadSelfile()
        if (DoUseGui):
            self.MakeGui()
        else:
            self.process_all_without_gui()
            

    def MakeGui(self):
        import gui_lib
        
        self.master=Tk()
        self.total_count=0
        self.guiwidth=600
        self.guiheight=500

        # Create the Canvas with Scrollbars
        self.canvas,self.frame=gui_lib.PrepareCanvas(self.master,self.guiwidth,self.guiheight)

        # Fill the GUI
        self.FillMarkGui()

        # Launch the window
        gui_lib.LaunchCanvas(self.canvas,self.frame)

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
 
    def process_all_without_gui(self):
        if not self.IsPairList:
            for micrograph,state in self.sellines:
                if (state.find('-1') == -1):
                    self.perform_picking(micrograph)
        else:
            for untilted,tilted,state in sellines:
                if (state.find('-1') == -1):
                    self.perform_picking_pair(untilted,tilted)

                for line in self.pairlines:
                    words=line.split()
                    untilted=words[0]
                    tilted=words[1]
                    state=words[2]
                    if (state.find('-1') == -1):
                        self.update_have_picked()
                        if not (self.have_already_picked(untilted)):
                            self.perform_picking_pair(untilted,tilted)

    def FillMarkGui(self):
        import os
        import gui_lib

        self.whichmark=StringVar()
        self.row={}

        # Script title
        self.master.title("GUI for particle picking")
        headertext='GUI for Xmipp particle picking\n'
        headertext+="Executed in directory: "+str(os.getcwd())
        self.l1=Label(self.frame, text=headertext, fg="medium blue")
        self.l1.grid(row=0, column=0,columnspan=3,sticky=EW)
 
        total=0
        if not self.IsPairList:
            for micrograph,state in self.sellines:
                if (state.find('-1') == -1):
                    c=self.CountPicked(micrograph)
                    total+=c
                    row=self.GuiAddSingleEntry(micrograph,c)
                    self.row[micrograph]=row
        else:
            for micrograph,tilted,state in sellines:
                if (state.find('-1') == -1):
                    self.GuiAddPairEntry(micrograph)
 
        row=(self.frame.grid_size()[1]+1)
        gui_lib.AddSeparator(self.frame,row,6)
        self.buttonrow=(self.frame.grid_size()[1]+1)
        b = Button(self.frame, text="Close", command=self.GuiClose)
        b.grid(row=self.buttonrow,column=0,sticky=W)
        b = Button(self.frame, text="Update Total Count:", command=self.GuiUpdateCount)
        b.grid(row=self.buttonrow,column=1)
        label=str(total).zfill(5)
        l = Label(self.frame, text=label)
        l.grid(row=self.buttonrow,column=2)
    
    def GuiAddSingleEntry(self,micrograph,count):
        import os
        row=(self.frame.grid_size()[1])
        label='Micrograph '+os.path.basename(micrograph)
        l=Label(self.frame, text=label)
        l.grid(row=row, column=0, sticky=E)
        label=str(count).zfill(5)
        l=Label(self.frame, text=label)
        l.grid(row=row, column=2)
        r=Radiobutton(self.frame,text="Mark!",variable=self.whichmark,
                           value=micrograph,indicatoron=0, command=self.LaunchSingleMark)
        r.grid(row=row, column=1,sticky=N)
        return row

    def CountPicked(self,micrograph):
        import os
        posfile=str(micrograph).replace('.raw','')+'.'+str(self.PosName)
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
        print "* Marking... "
        self.perform_picking(self.whichmark.get())
        self.GuiUpdateCount()

    def GuiUpdateCount(self):
        print "* Updating count..."
        total=self.CountAll()
        label=str(total).zfill(5)
        l=Label(self.frame, text=label)
        l.grid(row=self.buttonrow, column=2)

    def GuiClose(self):
        import sys
        print "* Closing..."
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
        os.chdir(directory)
        command='xmipp_micrograph_mark -i '+micrograph
        print '* ',command
        self.log.info(command)
        os.system(command)
        os.chdir(os.pardir)

    def perform_picking_pair(self,untilted,tilted):
        import os
        directory,uname=os.path.split(untilted)
        tname='../'+tilted
        os.chdir(directory)
        command='xmipp_micrograph_mark -i '+uname+' -tilted '+tname
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

	particle_pick=particle_pick_class(MicrographSelfile,
                                          IsPairList,
                                          PosName,
                                          DoUseGui,
                                          ProjectDir,
                                          LogDir)

	# close 
	particle_pick.close()

if __name__ == '__main__':

    import sys
    pick=particle_pick_class(MicrographSelfile,
                             IsPairList,
                             PosName,
                             DoUseGui,
                             ProjectDir,
                             LogDir)
