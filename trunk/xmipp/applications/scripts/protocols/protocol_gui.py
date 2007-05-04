#!/usr/bin/env python
import sys
import os
import string
from Tkinter import *
""" Guidelines for python script header formatting:

The idea of the GUI class is that it provides a graphical editor of
only the variables contained in the header of a python script, and
that it can also be used to launch the script.

Usage:
  python protocol_gui.py script.py

Where the header of script.py should be organized as follows:

Obligatory:
    * Include a {end-of-header} label at the end of the header
    * Variable declaration (visible in GUI): Variablename=XXX
          o If XXX is True or False, the variable is considered as a
            Boolean (yes/no button in GUI)
          o If XXX starts with \" or \', the variable is considered as
            a string (text field in GUI)
          o If XXX is a number, the variable is considered as a number
            (text field in GUI) 
    * The first single comment line above each variable (starting with a #)
      will be displayed in the GUI, as a label left from the variable entry field

    * More detailed help for each variable can be placed between the comment
      line and the variable declaration line using \""" ...\""".
      This text will become visible in the GUI by pressing a -what's this?-
      button, at the right side of the variable entry field.
      !!!NOTE that the first character in newlines within these comments
      should start with spaces or a tab!!!
    * An {expert} label on the comment line (#) marks the option as -expert-
      (by default not shown in the GUI, unless you press the -show expert options-
      button, at the left side of the variable entry field. Then, they will be
      shown in yellow.
    * A {file} or {dir} label on the comment line (#) marks the option as a
      Filename or Directory, and will add a corresponding Browse-button to the GUI
      Note that this button return absolute paths

Optional:

    * A {please cite} label on a comment line (starting with a #) will display
      a message at the top of the protocol stating -If you publish results obtained with
      this protocol, please cite-, followed by the text on rest of the comment line.
      If more than one citation lines are present, they will all be displayed.
      DONT use very long citations, as this will results in an ugly gui.
    * Include a {section} label on a comment line (starting with a #) to separate
      variables by a blue line + the corresponding title in the GUI 

"""
# Create a GUI automatically from a script
class automated_gui_class:

    def __init__(self,
                 scriptname):

        self.SYSTEMSCRIPTDIR=os.path.expanduser('~')+'/scripts/'

        self.scriptname=scriptname
        self.master=Tk()
        self.expert_mode=False
        self.is_setupgui=False
        self.is_markgui=False
        self.have_publication=False

    def MakeGui(self):
        self.ScriptRead()
        self.ScriptParseVariables()
        self.SetupGuiParameters()
        self.GuiFill()
        
    def SetupGuiParameters(self):
        if (self.is_setupgui):
            self.columnspantextlabel=2
            self.columntextentry=2
            self.column_pre=0
            self.column_2d=1
            self.column_3d=2
        else:
            self.columnspantextlabel=3
            self.columntextentry=3

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
                # Add "" for strings
                if (self.variables[args[0]][1]=="String" or
                    self.variables[args[0]][1]=="File" or
                    self.variables[args[0]][1]=="Directory"):
                    lineb=str(args[0])+'=\''+str(self.variables[args[0]][2].get())+"\'\n"
                else:
                    lineb=str(args[0])+'='+str(self.variables[args[0]][2].get())+"\n"
            else:
                lineb=self.script_header_lines[i]
            linesb.append(lineb)

        fh.writelines(linesb)
        fh.writelines(self.script_body_lines)
        fh.close() 
        os.chmod(self.scriptname,0755)

    def ScriptParseComments(self,i):
        import string
        morehelp=[]
        found_comment=False
        j=i
        # comment will be the first line above the i^th line that begins with a "#"
        while not found_comment:
            j-=1
            if (self.script_header_lines[j][0]=='#'):
                found_comment=True
                comment=self.script_header_lines[j][1:]
                # Check setup or expert options
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
                # Check browse options
                if "{file}" in comment:
                    browse="file"
                    comment=comment.replace ('{file}', '' )
                elif "{dir}" in comment:
                    browse="dir"
                    comment=comment.replace ('{dir}', '' )
                else:
                    browse="none"
        # Checkout more help
        if (i-j>1):
            while (j<i-1):
                j=j+1
                morehelp+=self.script_header_lines[j]

        morehelp=string.join(morehelp,'')
        morehelp=morehelp.replace('\"\"\"','')
        return comment,isexpert,morehelp,browse
        
    def Is_a_number(self,string):
        try:
            test=float(string)
        except ValueError:
            return False
        return True

    def ScriptParseVariables(self):
        self.variables={}
        self.vfields=[]
        self.publications=[]
        self.have_analyse_results=False

        for i in range(len(self.script_header_lines)):
            # Get section headers
            if ("{section}" in self.script_header_lines[i]):
                section=self.script_header_lines[i][1:].replace('{section}','')
                args=self.script_header_lines[i].split("}")
                self.variables[section]=[section[:-1],]
                self.vfields.append(section)
                self.variables[section].append("Section")
            # Get corresponding publication
            if ("{please cite}" in self.script_header_lines[i]):
                self.have_publication=True
                self.publications.append(self.script_header_lines[i][1:].replace('{please cite}',''))
            
            # Get a variable
            elif (self.script_header_lines[i][0]!="#"
                  and self.script_header_lines[i][0]!="\n"
                  and self.script_header_lines[i][0]!="\""
                  and self.script_header_lines[i][0]!="\t"
                  and self.script_header_lines[i][0]!=" "):
                args=self.script_header_lines[i].split("=")
                value=args[1][:-1]
                value=value.strip('\'')
                value=value.strip('\"')
                if (args[0]=="AnalysisScript"):
                    self.have_analyse_results=True;
                self.vfields.append(args[0])
                if ("True" in args[1]) or ("False" in args[1]):
                    # boolean
                    self.variables[args[0]]=[value,]
                    self.variables[args[0]].append("Boolean")
                    newvar=BooleanVar()                    
                    self.variables[args[0]].append(newvar)
                elif (self.Is_a_number(value)):
                    # number
                    self.variables[args[0]]=[value,]
                    self.variables[args[0]].append("Number")
                    newvar=StringVar()                    
                    self.variables[args[0]].append(newvar)
                else:
                    # string
                    self.variables[args[0]]=[value,]
                    self.variables[args[0]].append("String")
                    newvar=StringVar()
                    self.variables[args[0]].append(newvar)
            
                comment,isexpert,morehelp,browse=self.ScriptParseComments(i)
                if (browse=="file"):
                    self.variables[args[0]][1]="File"
                elif (browse=="dir"):
                    self.variables[args[0]][1]="Directory"
                self.variables[args[0]].append(comment)
                self.variables[args[0]].append(isexpert)
                self.variables[args[0]].append(morehelp)
                row=0
                self.variables[args[0]].append(row)

        if self.is_setupgui:
            # Set PWD as default for ProjectDir
            self.variables['ProjectDir'][0]=str(os.getcwd())

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

    def GuiResize(self,master,frame):
        height=frame.winfo_reqheight()+25
        width=frame.winfo_reqwidth()+25
        if (height>600):
            height=600
        if (width>800):
            width=800
        master.geometry("%dx%d%+d%+d" % (width,height,0,0))

    def LaunchCanvas(self,master,canvas,frame):
        # Launch the window
        canvas.create_window(0, 0, anchor=NW, window=frame)
        frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))
        
    def GuiFill(self):
                       
        # Create the Canvas with Scrollbars
        self.canvas,self.frame=self.PrepareCanvas(self.master)

        # Fill the entire GUI
        if (self.is_setupgui):
            self.FillSetupGui()
        else:
            self.FillProtocolGui()
        
        # Launch the window
        self.LaunchCanvas(self.master,self.canvas,self.frame)

        # Remove expert options in normal mode
        if (self.expert_mode==False):
            for w in self.widgetexpertlist:
                w.grid_remove()  

        # Resize window
        self.GuiResize(self.master,self.frame)

        # Enter main loop
        self.master.mainloop()

    def FillProtocolGui(self):

        import os,sys
        self.morehelp=StringVar()
        self.whichfile=StringVar()

        # Script title
        programname=(os.path.basename(sys.argv[1])).replace('.py','')
        self.master.title(programname)
        headertext='GUI for Xmipp '+programname+'\n'
        headertext+="Executed in directory: "+str(os.getcwd())
        self.l1=Label(self.frame, text=headertext, fg="medium blue")
        self.l1.grid(row=0, column=0,columnspan=6,sticky=EW)
        if (self.have_publication):
            headertext="If you publish results obtained with this protocol, please cite:"
            for pub in self.publications:
                headertext+='\n'+pub.replace('\n','')
            self.l2=Label(self.frame, text=headertext, fg="dark green")
            self.l2.grid(row=1, column=0,columnspan=5,sticky=EW)
            self.AddSeparator(2)
        else:
            self.AddSeparator(1)

        # Add all the variables in the script header
        self.widgetexpertlist=[]
        for var in self.vfields:
            if (self.variables[var][1]=="Section"):
                self.GuiAddSection(self.variables[var][0])
            elif (self.variables[var][1]=="File" or
                  self.variables[var][1]=="Directory"):                 
                if (self.variables[var][1]=="File"):
                    is_a_file=True
                else:
                    is_a_file=False
                self.GuiAddBrowseEntry(var,
                                       self.variables[var][3],
                                       self.variables[var][0],
                                       self.variables[var][2],
                                       self.variables[var][4],
                                       self.variables[var][5],
                                       is_a_file)
            elif (self.variables[var][1]=="String" or
                  self.variables[var][1]=="Number"):
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
        self.buttonrow=(self.frame.grid_size()[1])
        self.GuiAddRestProtocolButtons()

    def FillSetupGui(self):
        self.which_setup=StringVar()
        self.morehelp=StringVar()
      
        # Script title
        self.master.title('Xmipp protocols')
        headertext="Which Xmipp protocol do you want to run?"
        self.l1=Label(self.frame, text=headertext, fg="medium blue")
        self.l1.grid(row=0, column=0,columnspan=5,sticky=EW)
        self.AddSeparator(2)

        # Add labels for different protocol categories
        row=(self.frame.grid_size()[1])
        self.row_pre=row
        self.row_2d=row
        self.row_3d=row
        self.l=Label(self.frame, text="Preprocessing", fg="medium blue", width=30)
        self.l.grid(row=row, column=self.column_pre,columnspan=1, sticky=E)
        self.l=Label(self.frame, text="2D analysis", fg="medium blue", width=30)
        self.l.grid(row=row, column=self.column_2d,columnspan=1, sticky=E)
        self.l=Label(self.frame, text="3D analysis", fg="medium blue", width=30)
        self.l.grid(row=row, column=self.column_3d,columnspan=1, sticky=E)

        # Add all the variables in the script header
        self.widgetexpertlist=[]
        for var in self.vfields:
            if (self.variables[var][1]=="Section"):
                self.GuiAddSection(self.variables[var][0])
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
        self.buttonrow=(self.frame.grid_size()[1])
        self.GuiAddRestSetupButtons()

    def GuiAddSection(self,label):
        row=(self.frame.grid_size()[1])
        line="-----------------------------------------------------------"
        self.l1=Label(self.frame, text=label, fg="medium blue")
        self.l2=Label(self.frame, text=line, fg="medium blue")
        self.l1.grid(row=row, column=0,columnspan=self.columnspantextlabel,sticky=E)
        self.l2.grid(row=row+1, column=0,columnspan=self.columnspantextlabel,sticky=E)

    def AddSeparator(self,row):
        self.l1=Label(self.frame,text="")
        self.l1.grid(row=row)
        self.l2=Frame(self.frame, height=2, bd=1, bg="medium blue",relief=RIDGE)
        self.l2.grid(row=row+1, column=0,columnspan=self.columnspantextlabel+3,sticky=EW)
        self.l3=Label(self.frame,text="")
        self.l3.grid(row=row+2)

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

    def GuiPositionLabel(self,label,default,variable,expert,morehelp,has_browse=False):
        row=(self.frame.grid_size()[1])
        if (expert=="expert"):
            self.l=Label(self.frame, text=label, bg="yellow")
        else: 
            self.l=Label(self.frame, text=label)
        self.l.configure(wraplength=350)
        self.l.grid(row=row, column=0,columnspan=self.columnspantextlabel, sticky=E)
        self.r=Radiobutton(self.frame,text="What's this?",variable=self.morehelp,
                           value=morehelp,indicatoron=0, command=self.GuiShowMoreHelp )
        if (has_browse):
            column=self.columntextentry+3
        else:
            column=self.columntextentry+2
        if (morehelp!=""):
            self.r.grid(row=row, column=column, sticky=W)
            
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
            if (morehelp!=""):
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
            if (morehelp!=""):
                self.widgetexpertlist.append(self.r)

    def GuiAddBrowseEntry(self,varname,label,default,variable,expert,morehelp,is_file):
        row,self.l,self.r=self.GuiPositionLabel(label,default,variable,expert,morehelp,True)
        self.variables[varname][6]=row
        self.e = Entry(self.frame, text=label, textvariable=variable)
        self.e.delete(0, END) 
        self.e.insert(0,default)
        self.e.grid(row=row, column=self.columntextentry,columnspan=2,sticky=W+E)
        if (is_file):
            self.b = Radiobutton(self.frame, text="Browse",indicatoron=0,variable=self.whichfile,
                                 value=varname, command=self.GuiBrowseFile)
        else:
            self.b = Radiobutton(self.frame, text="Browse",indicatoron=0,variable=self.whichfile,
                                 value=varname, command=self.GuiBrowseDirectory)
        self.b.grid(row=row, column=self.columntextentry+2,sticky=W)

        if (expert=="expert"):
            self.widgetexpertlist.append(self.l)
            self.widgetexpertlist.append(self.e)
            self.widgetexpertlist.append(self.b)
            if (morehelp!=""):
                self.widgetexpertlist.append(self.r)

    def GuiAddRestProtocolButtons(self):
        self.AddSeparator(self.buttonrow)
        self.button = Button(self.frame, text="Close", command=self.GuiClose)
        self.button.grid(row=self.buttonrow+3,column=0, sticky=W)

        if (self.expert_mode==True):
            text2=" Hide expert options "
        else:
            text2="Show expert options"
        self.bGet = Button(self.frame, text=text2, command=self.GuiTockleExpertMode)
        self.bGet.grid(row=self.buttonrow+3,column=1)

        self.bGet = Button(self.frame, text="Load", command=self.GuiLoad)
        self.bGet.grid(row=self.buttonrow+3,column=2)
        self.bGet = Button(self.frame, text="Save", command=self.GuiSave)
        self.bGet.grid(row=self.buttonrow+3,column=3)
        self.bGet = Button(self.frame, text="Save & Execute", command=self.GuiSaveExecute)
        self.bGet.grid(row=self.buttonrow+3,column=4)
        if (self.have_analyse_results):
            self.bGet = Button(self.frame, text="Analyse Results", command=self.AnalyseResults)
            self.bGet.grid(row=self.buttonrow+3,column=5)

    def GuiAddRestSetupButtons(self):
        self.button = Button(self.frame, text="Close", command=self.master.quit)
        self.button.grid(row=self.buttonrow,column=0, sticky=W)
        if (self.expert_mode==True):
            text2=" Hide expert options "
        else:
            text2="Show expert options"
        self.bGet = Button(self.frame, text=text2, command=self.GuiTockleExpertMode)
        self.bGet.grid(row=self.buttonrow,column=1)
        if (self.have_publication):
            headertext="Read more about Xmipp protocols in:"
            for pub in self.publications:
                headertext+='\n'+pub.replace('\n','')
            self.l2=Label(self.frame, text=headertext, fg="dark green")
            self.l2.grid(row=self.buttonrow+1, column=0,columnspan=5,sticky=EW)

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

    def relpath(self,target, base=os.curdir):
        import os
        if not os.path.exists(target):
            raise OSError, 'Target does not exist: '+target
        if not os.path.isdir(base):
            raise OSError, 'Base is not a directory or does not exist: '+base
        base_list = (os.path.abspath(base)).split(os.sep)
        target_list = (os.path.abspath(target)).split(os.sep)
        for i in range(min(len(base_list), len(target_list))):
            if base_list[i] <> target_list[i]: break
        else:
            i+=1
        rel_list = [os.pardir] * (len(base_list)-i) + target_list[i:]
        return os.path.join(*rel_list)

    def GuiBrowseFile(self):
        import tkFileDialog
        import os
        fileformats = [('All Files ','*.*')]
        fname = tkFileDialog.askopenfilename(title='Choose File',
                                             filetypes=fileformats)
        if (len(fname)>0):
#            fname=self.relpath(fname,os.curdir)
            self.e = Entry(self.frame,textvariable=self.variables[self.whichfile.get()][2])
            self.e.delete(0, END) 
            self.e.insert(0,fname)
            self.e.grid(row=self.variables[self.whichfile.get()][6],
                        column=self.columntextentry,
                        columnspan=2,sticky=W+E)
            print "var=",self.whichfile.get(),self.variables[self.whichfile.get()][2].get()

    def GuiBrowseDirectory(self):
        import tkFileDialog
        import os
        fname = tkFileDialog.askdirectory()       
        if (len(fname)>0):
#            fname=self.relpath(fname,os.curdir)
            self.e = Entry(self.frame,textvariable=self.variables[self.whichfile.get()][2])
            self.e.delete(0, END) 
            self.e.insert(0,fname)
            self.e.grid(row=self.variables[self.whichfile.get()][6],
                        column=self.columntextentry,
                        columnspan=2,sticky=W+E)

    def AnalyseResults(self):
        self.GuiSave()
        print "* Analyzing..."
        command='python '+str(self.SYSTEMSCRIPTDIR)+'/protocol_gui.py '+\
                 self.variables["AnalysisScript"][0]+' &'
        print command
        os.system(command)
         
    def GuiClose(self):
        print "* Closing..."
        self.master.destroy()
        
    def GuiSave(self):
        print "* Saving..."
        self.ScriptWrite()

    def GuiSaveExecute(self):
        import tkMessageBox
        self.GuiSave()
        command="python "+self.scriptname+' &'
        answer=tkMessageBox._show("Execute protocol",
                                  "Use a job queueing system?",
                                  tkMessageBox.QUESTION, 
                                  tkMessageBox.YESNOCANCEL)
        if (answer=="yes"):
            self.master.update()
            d = MyQueueLaunch(self.master,command)
            self.master.wait_window(d.top)
        elif (answer=="no"):
            import popen2
            command="python "+self.scriptname+' &'
            print "* Executing job with: "+command
            os.system(command)
                
    def GuiLoad(self):
        import tkFileDialog
        import os,shutil
        fname = tkFileDialog.askdirectory()       
        if (len(fname)>0):
            print "* Loading protocol from "+os.path.basename(fname)+" ..."
            fname=fname+'/'+self.scriptname+'_backup'
            shutil.copy(fname,self.scriptname)
            self.master.destroy()
            self.master=Tk()
            self.ScriptRead()
            self.ScriptParseVariables()
            self.SetupGuiParameters()
            self.GuiFill()

    def GuiLanchSetup(self):
        self.GuiSave()
        print "* Launching protocol ..."
        command='python '+str(self.scriptname)+' '+str(self.which_setup.get())+' &'
        os.system(command)
       
    def GuiShowMoreHelp(self):
        import tkMessageBox
        message=str(self.morehelp.get())
        tkMessageBox.showinfo('More Help',message)

# A dialog window to ask for the queueing command
class MyQueueLaunch:
    def __init__(self, parent,command):
        self.command=command
        top = self.top = Toplevel(parent)
        Label(top, text="Job submission command \n (e.g. \"bsub -q 1week\")").grid(row=0,column=0,columnspan=2)
        self.e = Entry(top)
        self.e.grid(row=1,column=0,columnspan=2)
        Button(top, text="Submit", command=self.ok).grid(row=2,column=0)
        Button(top, text="Cancel", command=self.cancel).grid(row=2,column=1)
        

    def ok(self):
        import os
        command=self.e.get()+" "+self.command
        print "* Executing job with: "+command
        os.system(command)
        self.top.destroy()

    def cancel(self):
        self.top.destroy()
      

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
    automated_gui.MakeGui()
