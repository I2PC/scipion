#!/usr/bin/env python
import sys
import os
import string
from Tkinter import *
 
# Create a GUI automatically from a script
class automated_gui_class:

    def __init__(self,
                 scriptname):

        self.scriptname=scriptname
        self.MakeGui()
        
    def MakeGui(self):
        self.master=Tk()
        self.expert_mode=False
        self.is_mastergui=False
        self.ScriptRead()
        self.ScriptParseVariables()
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
            if "{is-master}" in line:
                self.is_mastergui=True
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
        found_comment=False
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
                else:
                    isexpert="normal"
        if (i-j>1):
            line="-----------------------------------------------------------\n"
            self.morehelp+=line
            self.morehelp+=comment
            self.morehelp+=line
            while (j<i-1):
                j=j+1
                self.morehelp+=self.script_header_lines[j]

        return comment,isexpert
        

    def ScriptParseVariables(self):
        self.variables={}
        self.vfields=[]
        self.morehelp=""

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
            
                comment,isexpert=self.ScriptParseComments(i)
                self.variables[args[0]].append(comment)
                self.variables[args[0]].append(isexpert)

    def GuiFill(self):
                       
        # Stuff to make the scrollbars work
        vscrollbar = AutoScrollbar(self.master)
        vscrollbar.grid(row=0, column=1, sticky=N+S)
        hscrollbar = AutoScrollbar(self.master, orient=HORIZONTAL)
        hscrollbar.grid(row=1, column=0, sticky=E+W)
        canvas = Canvas(self.master,
                        width=550,
                        height=650,
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
        
        # Script title
        headertext="GUI for Xmipp "
        programname=sys.argv[1]
        headertext+=programname.replace('.py','')
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
                                     self.variables[var][4])
            elif (self.variables[var][1]=="Boolean"):
                if (self.is_mastergui):
                    self.GuiAddMasterButton(self.variables[var][3])
                else:
                    newvar=BooleanVar()
                    self.variables[var].append(newvar)
                    self.GuiAddBooleanEntry(self.variables[var][3],
                                            self.variables[var][0],
                                            self.variables[var][2],
                                            self.variables[var][4])
            elif (self.variables[var][1]=="FileName"):
                newvar=StringVar()
                self.variables[var].append(newvar)
                self.GuiAddFileNameEntry(self.variables[var][3],
                                     self.variables[var][0],
                                     self.variables[var][2],
                                     self.variables[var][4])
            else:
                print "ERROR",self.variables[var][1]," variable type not recognized"
                exit

        # Add bottom row buttons
        self.buttonrow=(self.frame.grid_size()[1]+1)
        self.GuiAddRestButtons()

        # Launch the window
        canvas.create_window(0, 0, anchor=NW, window=self.frame)
        self.frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))

        # Remove expert options in normal mode
        if (self.expert_mode==False):
            for w in self.widgetexpertlist:
                w.grid_remove()  

        self.master.mainloop()


    def GuiAddSection(self,label,size):
        row=(self.frame.grid_size()[1]+1)
        line="-----------------------------------------------------------"
        if size==1:
            self.l1=Label(self.frame, text=label, font=("Helvetica", 16), fg="blue")
            self.l2=Label(self.frame, text=line, font=("Helvetica", 16), fg="blue")
            self.l1.grid(row=row, column=0,columnspan=5,sticky=EW)
            self.l2.grid(row=row+1, column=0,columnspan=5,sticky=EW)
        if size==2:
            self.l1=Label(self.frame, text=label, fg="blue")
            self.l2=Label(self.frame, text=line, fg="blue")
            self.l1.grid(row=row, column=0,columnspan=3,sticky=E)
            self.l2.grid(row=row+1, column=0,columnspan=3,sticky=E)

    def GuiAddMasterButton(self,label):
        row=(self.frame.grid_size()[1]+1)
        self.bGet = Button(self.frame, text=label, command=self.GuiSaveExecute)
        self.bGet.grid(row=row,columnspan=3,sticky=E)

    def GuiPositionLabel(self,label,default,variable,expert):
        row=(self.frame.grid_size()[1]+1)
        if (expert=="expert"):
            self.l=Label(self.frame, text=label,bg="yellow")
        else:
            self.l=Label(self.frame, text=label)
        self.l.configure(wraplength=350)
        self.l.grid(row=row, column=0,columnspan=3, sticky=E)
        return row,self.l

    def GuiAddBooleanEntry(self,label,default,variable,expert):
        row,self.l=self.GuiPositionLabel(label,default,variable,expert)
        self.r1 = Radiobutton(self.frame, text="Yes", variable=variable, value=True)
        self.r1.grid(row=row, column=3)
        self.r2 = Radiobutton(self.frame, text="No", variable=variable, value=False)
        self.r2.grid(row=row, column=4)
        if (default=="True"):
            self.r1.select()
        else:
            self.r2.select()
        if (expert=="expert"):
            self.widgetexpertlist.append(self.l)
            self.widgetexpertlist.append(self.r1)
            self.widgetexpertlist.append(self.r2)

    def GuiAddTextEntry(self,label,default,variable,expert):
        row,self.l=self.GuiPositionLabel(label,default,variable,expert)
        self.e = Entry(self.frame, text=label, textvariable=variable)
        self.e.delete(0, END) 
        self.e.insert(0,default)
        self.e.grid(row=row, column=3,columnspan=2,sticky=W+E)
        if (expert=="expert"):
            self.widgetexpertlist.append(self.l)
            self.widgetexpertlist.append(self.e)

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

    def GuiAddRestButtons(self):
        if (self.is_mastergui==False):
            self.bGet = Button(self.frame, text="Save & Execute", command=self.GuiSaveExecute)
            self.bGet.grid(row=self.buttonrow,column=4)
            self.bGet = Button(self.frame, text="Save", command=self.GuiSave)
            self.bGet.grid(row=self.buttonrow,column=3)
        if (self.expert_mode==True):
            text2=" Hide expert options "
        else:
            text2="Show expert options"
        self.bGet = Button(self.frame, text=text2, command=self.GuiTockleExpertMode)
        self.bGet.grid(row=self.buttonrow,column=0)
        self.button = Button(self.frame, text="QUIT", command=self.master.quit)
        self.button.grid(row=self.buttonrow,column=2)
        self.hi_there = Button(self.frame, text="More Help", command=self.PrintHelp)
        self.hi_there.grid(row=self.buttonrow,column=1)

    def PrintHelp(self):
        print self.morehelp

    def GuiTockleExpertMode(self):
        if (self.expert_mode==True):
            for w in self.widgetexpertlist:
                w.grid_remove()
            self.expert_mode=False
        else:
            self.expert_mode=True
            for w in self.widgetexpertlist:
                w.grid()  
        self.GuiAddRestButtons()

    def GuiSaveExecute(self):
        self.GuiSave()
        self.GuiExecute()
        
    def GuiExecute(self):
        print "Executing..."
        os.system('python '+self.scriptname+' &')
        
    def GuiSave(self):
        print "Saving..."
        self.ScriptWrite()
        
# A scrollbar that hides itself if it's not needed.
class AutoScrollbar(Scrollbar):
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)


if __name__ == '__main__':

    args=sys.argv[1]
    automated_gui=automated_gui_class(args)
