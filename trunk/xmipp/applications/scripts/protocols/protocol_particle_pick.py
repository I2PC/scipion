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
MicrographSelfile='/home/scheres/work/protocols/all.sel'
# Is this selfile a list of untilted-tilted pairs?
""" True for RCT-processing. In that case, provide a 3-column selfile as follows:
    untilted_pair1.raw tilted_pair1.raw 1
    untilted_pair2.raw tilted_pair2.raw 1
    etc...
    Where 1 in the third column means active pair, and -1 means inactive pair
"""
IsPairList=True
# Name of the position files (or family name)
""" This is specified inside the micrograph_mark program (raw.Common.pos by default)
"""
PosName='raw.Common.pos'
# Use GUI to display list of finished micrographs?
DoUseGui=True
# {expert} Root directory name for this project:
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
# A scrollbar that hides itself if it's not needed.
class AutoScrollbar(Scrollbar):
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)

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

    def MakeGui(self):
        self.master=Tk()
        self.total_count=0
        self.whichmark=StringVar()
        self.whichtilted={}
        self.row={}

        # Create the Canvas with Scrollbars
        self.canvas,self.frame=self.PrepareCanvas(self.master)

        # Fill the GUI
        self.FillMarkGui()

        # Launch the window
        self.LaunchCanvas(self.master,self.canvas,self.frame)
        self.GuiResize(self.master,self.frame)

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
 
    def PrepareCanvas(self,master):

        # Stuff to make the scrollbars work
        vscrollbar = AutoScrollbar(master)
        vscrollbar.grid(row=0, column=1, sticky=N+S)
        hscrollbar = AutoScrollbar(master, orient=HORIZONTAL)
        hscrollbar.grid(row=1, column=0, sticky=E+W)
        canvas = Canvas(master,
                        yscrollcommand=vscrollbar.set,
                        xscrollcommand=hscrollbar.set)
        canvas.grid(row=0, column=0, sticky=N+S+E+W)
        vscrollbar.config(command=canvas.yview)
        hscrollbar.config(command=canvas.xview)
        master.grid_rowconfigure(0, weight=1)
        master.grid_columnconfigure(0, weight=1)
        frame = Frame(canvas)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        return canvas,frame

    def LaunchCanvas(self,master,canvas,frame):
        # Launch the window
        canvas.create_window(0, 0, anchor=NW, window=frame)
        frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))
        
    def GuiResize(self,master,frame):
        height=frame.winfo_reqheight()+25
        width=frame.winfo_reqwidth()+25
        if (height>600):
            height=600
        if (width>800):
            width=800
        master.geometry("%dx%d%+d%+d" % (width,height,0,0))

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
                    row=self.GuiAddPairMarkEntry(micrograph,c)
                    self.row[micrograph]=row
                    self.whichtilted[micrograph]=tilted

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
        label=os.path.basename(micrograph)
        l=Label(self.frame, text=label)
        l.grid(row=row, column=0, sticky=E)
        label=str(count).zfill(5)
        l=Label(self.frame, text=label)
        l.grid(row=row, column=2)
        print "mic=",micrograph
        r=Radiobutton(self.frame,text="Mark!",variable=self.whichmark,
                           value=micrograph, indicatoron=0, command=self.LaunchPairMark)
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

    def LaunchPairMark(self):
        import os
        print "* Marking... "
        untilted=self.whichmark.get()
        print "untilted=",untilted
        print "tilted=",self.whichtilted[untilted]
        self.perform_picking_pair(self.whichmark.get(),self.whichtilted[self.whichmark.get()])
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
