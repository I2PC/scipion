#!/usr/bin/env python
import sys
import os
import string
from Tkinter import *
 
# Create a GUI automatically from a script
class automated_gui_class:

    def __init__(self,
                 scriptname):

        self.SYSTEMSCRIPTDIR=os.path.expanduser('~')+'/scripts/'

        self.scriptname=scriptname
        self.master=Tk()
        self.expert_mode=False
        self.is_setupgui=False
        
        self.ScriptRead()
        self.ScriptParseVariables()
        self.SetupGuiParameters()
        self.GuiFill()
        
    def ScriptRead(self):
        fh=open(self.scriptname,'r')
        self.script_header_lines=[]
        self.script_body_lines=[]
        isheader=False
        while (isheader!=True):
            line=fh.readline()
            if line=="":
                print "Error, this file does not have a {end-of-header} label"
                sys.exit()
            if "{setup-" in line:
                self.is_setupgui=True
            self.script_header_lines.append(line)
            if "{end-of-header}" in line:
                isheader=True

        self.script_body_lines=fh.readlines()
        fh.close()

    def ScriptWrite(self):
        fh=open(self.scriptname,'w')
        linesb=[]
        for i in range(len(self.script_header_lines)):
            if (self.script_header_lines[i][0]!="#" and
                self.script_header_lines[i][0]!="\"" and
                self.script_header_lines[i][0]!="\n" and
                self.script_header_lines[i][0]!=" " and
                self.script_header_lines[i][0]!="\t"):
                args=self.script_header_lines[i].split("=")
                lineb=str(args[0])+'='+str(self.variables[args[0]][2].get())+"\n"
            else:
                lineb=self.script_header_lines[i]
            linesb.append(lineb)

        fh.writelines(linesb)
        fh.writelines(self.script_body_lines)
        fh.close()        

    def ScriptParseComments(self,i):
        import string
        found_comment=False
        morehelp=[]
        j=i
        # comment will be the first line above the i^th line that begins with a "#"
        while not found_comment:
            j-=1
            if (self.script_header_lines[j][0]=='#'):
                found_comment=True
                comment=self.script_header_lines[j][1:]
                if "{expert}" in comment:
                    isexpert="expert"
                    comment=comment.replace ('{expert}', '' )
                elif "{setup-pre}" in comment:
                    isexpert="setup-pre"
                    comment=comment.replace ('{setup-pre}', '' )
                elif "{setup-2d}" in comment:
                    isexpert="setup-2d"
                    comment=comment.replace ('{setup-2d}', '' )
                elif "{setup-3d}" in comment:
                    isexpert="setup-3d"
                    comment=comment.replace ('{setup-3d}', '' )
                else:
                    isexpert="normal"
        if (i-j>1):
            while (j<i-1):
                j=j+1
                morehelp+=self.script_header_lines[j]

        morehelp=string.join(morehelp,'')
        morehelp=morehelp.replace('\"\"\"','')
        return comment,isexpert,morehelp
        

    def ScriptParseVariables(self):
        self.variables={}
        self.vfields=[]
        self.have_analyse_results=False

        for i in range(len(self.script_header_lines)):
            # Get section headers
            if ("{section}" in self.script_header_lines[i]):
                section=self.script_header_lines[i][1:].replace('{section}','')
                args=self.script_header_lines[i].split("}")
                self.variables[section]=[section[:-1],]
                self.vfields.append(section)
                self.variables[section].append("Section")
            # Get a variable
            elif (self.script_header_lines[i][0]!="#"
                  and self.script_header_lines[i][0]!="\n"
                  and self.script_header_lines[i][0]!="\""
                  and self.script_header_lines[i][0]!="\t"
                  and self.script_header_lines[i][0]!=" "):
                args=self.script_header_lines[i].split("=")
                if (args[0]=="AnalysisScript"):
                    self.have_analyse_results=True;
                self.variables[args[0]]=[args[1][:-1],]
                self.vfields.append(args[0])
                if ("True" in args[1]) or ("False" in args[1]):
                    self.variables[args[0]].append("Boolean")
                    newvar=BooleanVar()                    
                    self.variables[args[0]].append(newvar)
                else:
                    self.variables[args[0]].append("String")
                    newvar=StringVar()
                    self.variables[args[0]].append(newvar)
            
                comment,isexpert,morehelp=self.ScriptParseComments(i)
                self.variables[args[0]].append(comment)
                self.variables[args[0]].append(isexpert)
                self.variables[args[0]].append(morehelp)


    def SetupGuiParameters(self):
        if (self.is_setupgui):
            self.guiwidth=680
            self.guiheight=350
            self.columnspantextlabel=2
            self.columntextentry=2
            self.column_pre=0
            self.column_2d=1
            self.column_3d=2
        else:
            self.guiwidth=700
            self.guiheight=650
            self.columnspantextlabel=3
            self.columntextentry=3

    def GuiFill(self):
                       
        # Stuff to make the scrollbars work
        vscrollbar = AutoScrollbar(self.master)
        vscrollbar.grid(row=0, column=1, sticky=N+S)
        hscrollbar = AutoScrollbar(self.master, orient=HORIZONTAL)
        hscrollbar.grid(row=1, column=0, sticky=E+W)
        canvas = Canvas(self.master,
                        width=self.guiwidth,
                        height=self.guiheight,
                        yscrollcommand=vscrollbar.set,
                        xscrollcommand=hscrollbar.set)
        canvas.grid(row=0, column=0, sticky=N+S+E+W)
        vscrollbar.config(command=canvas.yview)
        hscrollbar.config(command=canvas.xview)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_columnconfigure(0, weight=1)
        self.frame = Frame(canvas)
        self.frame.rowconfigure(1, weight=1)
        self.frame.columnconfigure(1, weight=1)

        # Fill the entire GUI
        if (self.is_setupgui):
            self.FillSetupGui()
        else:
            self.FillProtocolGui()
        
        # Launch the window
        canvas.create_window(0, 0, anchor=NW, window=self.frame)
        self.frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))

        # Remove expert options in normal mode
        if (self.expert_mode==False):
            for w in self.widgetexpertlist:
                w.grid_remove()  

        self.master.mainloop()


    def FillProtocolGui(self):

        self.morehelp=StringVar()

        # Script title
        headertext="GUI for Xmipp "
        programname=sys.argv[1]
        headertext+=programname.replace('.py','')
        self.GuiAddSection(headertext,0)

        headertext="Executed in directory: "+str(os.getcwd())
        self.GuiAddSection(headertext,1)

        # Add all the variables in the script header
        self.widgetexpertlist=[]
        for var in self.vfields:
            if (self.variables[var][1]=="Section"):
                self.GuiAddSection(self.variables[var][0],2)
            elif (self.variables[var][1]=="String"):
                self.GuiAddTextEntry(self.variables[var][3],
                                     self.variables[var][0],
                                     self.variables[var][2],
                                     self.variables[var][4],
                                     self.variables[var][5])
            elif (self.variables[var][1]=="Boolean"):
                newvar=BooleanVar()
                self.variables[var].append(newvar)
                self.GuiAddBooleanEntry(self.variables[var][3],
                                        self.variables[var][0],
                                        self.variables[var][2],
                                        self.variables[var][4],
                                        self.variables[var][5])
            else:
                print "ERROR",self.variables[var][1]," variable type not recognized"
                sys.exit()

        # Add bottom row buttons
        self.buttonrow=(self.frame.grid_size()[1]+1)
        self.GuiAddRestProtocolButtons()


    def FillSetupGui(self):

        self.which_setup=StringVar()
        self.morehelp=StringVar()
      
        # Add labels for different protocol categories
        self.GuiAddSection("Which Xmipp protocol do you want to run?",0)
        row=(self.frame.grid_size()[1]+1)
        self.row_pre=row
        self.row_2d=row
        self.row_3d=row
        self.l=Label(self.frame, text="Preprocessing", fg="blue", width=30)
        self.l.grid(row=row, column=self.column_pre,columnspan=1, sticky=E)
        self.l=Label(self.frame, text="2D analysis", fg="blue", width=30)
        self.l.grid(row=row, column=self.column_2d,columnspan=1, sticky=E)
        self.l=Label(self.frame, text="3D analysis", fg="blue", width=30)
        self.l.grid(row=row, column=self.column_3d,columnspan=1, sticky=E)

        # Add all the variables in the script header
        self.widgetexpertlist=[]
        for var in self.vfields:
            if (self.variables[var][1]=="Section"):
                self.GuiAddSection(self.variables[var][0],2)
            elif (self.variables[var][1]=="String"):
                self.GuiAddTextEntry(self.variables[var][3],
                                     self.variables[var][0],
                                     self.variables[var][2],
                                     self.variables[var][4],
                                     self.variables[var][5])
            elif (self.variables[var][1]=="Boolean"):
                self.GuiAddLaunchButton(self.variables[var][3],
                                        var,
                                        self.variables[var][4])
            else:
                print "ERROR",self.variables[var][1]," variable type not recognized"
                sys.exit()

        # Add bottom row buttons
        self.buttonrow=(self.frame.grid_size()[1]+1)
        self.GuiAddRestSetupButtons()

    def GuiAddSection(self,label,size):
        row=(self.frame.grid_size()[1]+1)
        line="-----------------------------------------------------------"
        if size==0:
            self.l1=Label(self.frame, text=label, font=("Helvetica", 18), fg="blue")
            self.l2=Label(self.frame, text=line, font=("Helvetica", 18), fg="blue")
            self.l1.grid(row=row, column=0,columnspan=5,sticky=EW)
            self.l2.grid(row=row+1, column=0,columnspan=5,sticky=EW)
        if size==1:
            self.l1=Label(self.frame, text=label, font=("Helvetica", 14), fg="blue")
            self.l2=Label(self.frame, text="", font=("Helvetica", 14), fg="blue")
            self.l1.grid(row=row, column=0,columnspan=5,sticky=EW)
            self.l2.grid(row=row+1, column=0,columnspan=5,sticky=EW)
        if size==2:
            self.l1=Label(self.frame, text=label, fg="blue")
            self.l2=Label(self.frame, text=line, fg="blue")
            self.l1.grid(row=row, column=0,columnspan=self.columnspantextlabel,sticky=E)
            self.l2.grid(row=row+1, column=0,columnspan=self.columnspantextlabel,sticky=E)

    def GuiAddLaunchButton(self,label,value,expert):
        if (expert=="setup-pre"):
            column=self.column_pre
            self.row_pre+=1
            row=self.row_pre
        elif (expert=="setup-2d"):
            column=self.column_2d
            self.row_2d+=1
            row=self.row_2d
        elif (expert=="setup-3d"):
            column=self.column_3d
            self.row_3d+=1
            row=self.row_3d

        self.bGet = Radiobutton(self.frame, text=label, variable=self.which_setup, width=30,
                                value=value,  indicatoron=0, command=self.GuiLanchSetup)
        self.bGet.grid(row=row,column=column)

    def GuiPositionLabel(self,label,default,variable,expert,morehelp):
        row=(self.frame.grid_size()[1]+1)
        if (expert=="expert"):
            bg="yellow"
        else:
            bg="white"
        self.l=Label(self.frame, text=label, bg=bg)
        self.l.configure(wraplength=350)
        self.l.grid(row=row, column=0,columnspan=self.columnspantextlabel, sticky=E)
        self.r=Radiobutton(self.frame,text="More help",variable=self.morehelp, bg=bg,
                           value=morehelp,indicatoron=0, command=self.GuiShowMoreHelp )
        if (morehelp!=""):
            self.r.grid(row=row, column=self.columntextentry+2, sticky=W)
            
        return row,self.l,self.r

    def GuiAddBooleanEntry(self,label,default,variable,expert,morehelp):
        row,self.l,self.r=self.GuiPositionLabel(label,default,variable,expert,morehelp)
        self.r1 = Radiobutton(self.frame, text="Yes", variable=variable, value=True)
        self.r1.grid(row=row, column=self.columntextentry)
        self.r2 = Radiobutton(self.frame, text="No", variable=variable, value=False)
        self.r2.grid(row=row, column=self.columntextentry+1)
        if (default=="True"):
            self.r1.select()
        else:
            self.r2.select()
        if (expert=="expert"):
            self.widgetexpertlist.append(self.l)
            self.widgetexpertlist.append(self.r1)
            self.widgetexpertlist.append(self.r2)
            self.widgetexpertlist.append(self.r)

    def GuiAddTextEntry(self,label,default,variable,expert,morehelp):
        row,self.l,self.r=self.GuiPositionLabel(label,default,variable,expert,morehelp)
        self.e = Entry(self.frame, text=label, textvariable=variable)
        self.e.delete(0, END) 
        self.e.insert(0,default)
        self.e.grid(row=row, column=self.columntextentry,columnspan=2,sticky=W+E)
        if (expert=="expert"):
            self.widgetexpertlist.append(self.l)
            self.widgetexpertlist.append(self.e)
            self.widgetexpertlist.append(self.r)

    def GuiBrowseWindow(self):
        import tkFileDialog
        fileformats = [('All Files ','*.*')]
        self.FILENAME = tkFileDialog.askopenfilename(title='Choose File',
                                                     filetypes=fileformats)
        # This somehow doesn't work yet...
        self.e.delete(0, END) 
        self.e.insert(0,self.FILENAME)

    def GuiAddFileNameEntry(self,label,default,variable,expert):
        row,self.l=self.GuiPositionLabel(label,default,variable,expert)
        self.e = Entry(self.frame, text=label, textvariable=variable)
        self.e.delete(0, END) 
        self.e.insert(0,default)
        self.e.grid(row=row, column=3)
        self.FILENAME=StringVar()
        self.b = Button(self.frame, text="Browse", command=self.GuiBrowseWindow)
        self.b.grid(row=row, column=4)
        if (expert=="expert"):
            self.widgetexpertlist.append(self.l)
            self.widgetexpertlist.append(self.e)
            self.widgetexpertlist.append(self.b)

    def GuiAddRestProtocolButtons(self):
        if (self.expert_mode==True):
            text2=" Hide expert options "
        else:
            text2="Show expert options"
        self.bGet = Button(self.frame, text=text2, command=self.GuiTockleExpertMode)
        self.bGet.grid(row=self.buttonrow,column=0)
        self.button = Button(self.frame, text="QUIT", command=self.master.quit)
        self.button.grid(row=self.buttonrow,column=1)

        if (self.have_analyse_results):
            self.bGet = Button(self.frame, text="Analyse Results", command=self.AnalyseResults)
            self.bGet.grid(row=self.buttonrow,column=2)
        self.bGet = Button(self.frame, text="Save", command=self.GuiSave)
        self.bGet.grid(row=self.buttonrow,column=3)
        self.bGet = Button(self.frame, text="Save & Execute", command=self.GuiSaveExecute)
        self.bGet.grid(row=self.buttonrow,column=4)

    def GuiAddRestSetupButtons(self):
        if (self.expert_mode==True):
            text2=" Hide expert options "
        else:
            text2="Show expert options"
        self.bGet = Button(self.frame, text=text2, command=self.GuiTockleExpertMode)
        self.bGet.grid(row=self.buttonrow,column=0)
        self.button = Button(self.frame, text="QUIT", command=self.master.quit)
        self.button.grid(row=self.buttonrow,column=2)

    def GuiTockleExpertMode(self):
        if (self.expert_mode==True):
            for w in self.widgetexpertlist:
                w.grid_remove()
            self.expert_mode=False
        else:
            self.expert_mode=True
            for w in self.widgetexpertlist:
                w.grid()  
        if (self.is_setupgui):
            self.GuiAddRestSetupButtons()
        else:
            self.GuiAddRestProtocolButtons()

    def AnalyseResults(self):
        print "* Analyzing..."
        command='python '+str(self.SYSTEMSCRIPTDIR)+'/protocol_gui.py '+\
                 self.variables["WorkingDir"][0]+'/'+self.variables["AnalysisScript"][0]+' &'
        print command
        os.system(command)
         
    def GuiSaveExecute(self):
        self.GuiSave()
        self.GuiExecute()
        
    def GuiExecute(self):
        print "* Executing..."
        os.system('python '+self.scriptname+' &')
        
    def GuiSave(self):
        print "* Saving..."
        self.ScriptWrite()

    def GuiLanchSetup(self):
        import protocol_setup
        print "* Saving setup ..."
        self.ScriptWrite()
        print "* Launching protocol ..."
        command='python '+str(self.scriptname)+' '+str(self.which_setup.get())+' &'
        os.system(command)
       
    def GuiShowMoreHelp(self):
        import tkMessageBox
        message=str(self.morehelp.get())
        print message
        tkMessageBox.showinfo('More Help',message)
        
       
# A scrollbar that hides itself if it's not needed.
class AutoScrollbar(Scrollbar):
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)


if __name__ == '__main__':

    import sys
    args=sys.argv[1]
    automated_gui=automated_gui_class(args)
