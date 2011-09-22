#!bin/xmipp_python
'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

from Tkinter import *
import tkFont as font
import ttk
import os
from subprocess import Popen

from protlib_gui_ext import FilePollTextArea, centerWindows
    
class OptionsTab(Frame):
    def __init__(self, master, text, **opts):
        Frame.__init__(self, master)
        self.config(**opts)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=3)
        self.lastRow = 0
        self.options = []
        self.optionsDict = {}
        style = ttk.Style()
        style.configure("Name.TLabel", foreground="red")
        self.normal = font.Font(family='Helvetica', size=10)
        self.bold = font.Font(family='Helvetica', size=10, weight='bold')
        
    def addOption(self, name, comment, default='', cond=None):
        r = self.lastRow
        var = StringVar()
        var.set(default)
        ttk.Label(self, text=comment, font=self.normal).\
        grid(column=0, row=r, padx=5, pady=5, sticky=E)
        ttk.Label(self, text=name, font=self.bold).\
        grid(column=1, row=r, padx=5, pady=5, sticky=E)
        if default in ['yes', 'no']: # Boolean option
            w = ttk.Checkbutton(self, textvariable=var, variable=var, 
                            onvalue='yes', offvalue='no',
                            command=lambda:self.checked(name, var))
            w.grid(column=2, row=r, padx=5, pady=5, sticky=W)
        else:
            w = ttk.Entry(self, width=20, textvariable=var)
            if cond:
                if self.optionsDict[cond].get() == "no":
                    w['state'] = 'disabled'
            w.grid(column=2, row=r, sticky=(W, E), padx=5, pady=5)
        self.options.append((name, default, var, w, cond))
        self.optionsDict[name] = var
        self.lastRow += 1
        
    def setValue(self, name, value):
        self.optionsDict[name].set(value)
        
    def getValue(self, name):
        return self.optionsDict[name].get()
        
    def checked(self, name, var):
        value = var.get()
        if value == 'yes':
            state = 'normal'
        else:
            state = 'disabled'        
        for n, d, v, w, cond in self.options:
            if name == cond:
                w['state'] = state
                
    def addSeparator(self):
        ttk.Separator(self, orient=HORIZONTAL).\
        grid(column=0, row=self.lastRow, columnspan=3, 
             sticky=(W,E), padx=10, pady=5)
        self.lastRow += 1
        
    def getConfigOptions(self):
        optStr = ""
        for n, d, v, w, cond in self.options:
            if v.get() != d:
                optStr += ' %s="%s"' % (n, v.get())
        return optStr
        
class ConfigNotebook(ttk.Notebook):
    def __init__(self, master):
        ttk.Notebook.__init__(self, master)
        self.tabs = []
        
    def addTab(self, text):
        tab = OptionsTab(self, text)
        self.add(tab, text=text)
        self.tabs.append(tab)
        return tab
    
    def getConfigOptions(self):
        '''return string with configs settings '''
        return ' '.join([t.getConfigOptions() for t in self.tabs])


def find(file, pathlist):
  import os
  from os.path import join, getsize
  for path in pathlist:
      for root, dirs, files in os.walk(path):
        if file in files:
          return root
  return None
  
def detectJava():
    from subprocess import Popen, PIPE
    cmd = "readlink -f `which javac`"
    ps = Popen(cmd, shell=True, stdout=PIPE)
    out, err = ps.communicate()
    if err:
        return ''
    else:
        import os
        dirname = os.path.dirname        
        java_home = dirname(dirname(out.strip()))
        return java_home
    
def detectMpi():
    inc_dirs = ['/usr/include', 'usr/local/include']
    lib_dirs = ['/usr/lib64', 'usr/lib']
    inc_mpi = find('mpi.h', inc_dirs + lib_dirs)
    lib_mpi = find('libmpi.so', lib_dirs)
    
    return (inc_mpi, lib_mpi)

def detectQt():
    pass

root = Tk()
root.withdraw()
root.title("Xmipp Install")
root.minsize(width=600, height=350)

nb = ConfigNotebook(root)

tab = nb.addTab("Compiler")
tab.addOption('CC', 'The C compiler', 'gcc')
tab.addOption('CXX', 'The C++ compiler', 'g++')
tab.addOption('LINKERFORPROGRAMS', 'Linker for programs', 'g++')
tab.addOption('CCFLAGS', 'The C compiler flags', '')
tab.addOption('CXXFLAGS', 'The C++ compiler flags', '')

tab = nb.addTab("  MPI  ")
tab.addOption('mpi', 'Build the MPI programs?', 'no')
tab.addOption('MPI_CC', 'MPI C compiler', 'mpicc', cond='mpi')
tab.addOption('MPI_CXX', 'MPI C++ compiler', 'mpiCC', cond='mpi')
tab.addOption('MPI_LINKERFORPROGRAMS', 'MPI Linker for programs', 'mpiCC', cond='mpi')
tab.addOption('MPI_INCLUDE', 'MPI headers dir ', '/usr/include', cond='mpi')
tab.addOption('MPI_LIBDIR', 'MPI libraries dir ', '/usr/lib', cond='mpi')
tab.addOption('MPI_LIB', 'MPI library', 'mpi', cond='mpi')
inc_mpi, lib_mpi = detectMpi()
if inc_mpi:
    tab.setValue('MPI_INCLUDE', inc_mpi)
if lib_mpi:
    tab.setValue('MPI_LIBDIR', lib_mpi)

tab = nb.addTab("Java & QT")
tab.addOption('java', 'Build the java programs?', 'no')
tab.addOption('JAVAC', 'Java compiler', 'javac', cond='java')
tab.addOption('JAVA_HOME', 'Java installation directory', '', cond='java')
java_home = detectJava()
tab.setValue('JAVA_HOME', java_home)
tab.addSeparator()
tab.addOption('qt', 'Build the GUI (qt) programs?', 'yes')
tab.addOption('QTDIR', 'Where is QT installed', '/usr/share/qt3', cond='qt')
tab.addOption('QT_LIB', 'QT library to use', 'qt-mt', cond='qt')
tab.addOption('QT4', 'Use Qt4 instead of Qt3?', 'no', cond='qt')
tabJava = tab

tab = nb.addTab("Advanced")
tab.addOption('debug', 'Build debug version?', 'no')
tab.addOption('profile', 'Build profile version?', 'no')
tab.addOption('warn', 'Show warnings?', 'no')
tab.addOption('fast', 'Fast?', 'no')
tab.addOption('static', 'Prevent dynamic linking?', 'no')
tab.addOption('prepend', 'What to prepend to executable names', 'xmipp')

tab = nb.addTab("Configure")

output = 'build/scons_output.log'
if os.path.exists(output):
    os.remove(output)
text = FilePollTextArea(tab, output, 20, 80, colorOn=False)
text.goEnd()
text.grid(column=0, row=0, sticky=(N, S, E, W), padx=10, pady=10)

root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=5)
root.rowconfigure(0, weight=1)
leftFrame = ttk.Frame(root)
img = PhotoImage(file='resources/xmipp_logo.gif')
label = ttk.Label(leftFrame, image=img)
label['image'] = img
label.grid(column=0, row=0)
leftFrame.grid(column=0, row=0, sticky=(N, S), padx=5, pady=5, rowspan=2)
nb.grid(column=1, row=0, sticky=(N,W,E,S), padx=5, pady=5)

panel = ttk.Frame(root)
panel.grid(column=1, row=1, padx=5, pady=5, sticky=[W,E])
progressVar = IntVar()
progressVar.set(0)
progress = ttk.Progressbar(panel, orient=HORIZONTAL, length=300, mode='determinate', variable=progressVar, maximum="500")
progress.pack(side=LEFT, padx=(0, 30))
from protlib_gui_ext import MyButton, registerCommonFonts
registerCommonFonts()
btn = MyButton(panel, text='Compile')
btn.pack(side=RIGHT,padx=(15, 0))
procVar = IntVar()
procVar.set(2)
procEntry = ttk.Entry(panel, width="5", textvariable=procVar)
procEntry.pack(side=RIGHT,padx=5)
Label(panel, text="Processors").pack(side=RIGHT,padx=5)

proc = None

def setState(compiling):
    if compiling:
        btn.config(command=stopCompile, text="Stop")
    else:
        btn.config(command=runCompile, text="Compile") 
        
def checkProcess():
    print proc.pid
    if proc.poll() is not None:
        text.stopRefresh()
        text.fillTextArea(goEnd=True)
        setState(compiling=False)
        progressVar.set(0)
    else:
        root.after(3000, checkProcess)
        lines = 0
        for line in open(output):
            lines += 1
        progressVar.set(lines)
    
def runCompile(event=None):
    out = output
    procs = procVar.get()
    opts = nb.getConfigOptions()
    cmd = 'echo "*** RUNNING SCONS.CONFIGURE..." > %(out)s \n'
    cmd1 = "./scons.configure %(opts)s >> %(out)s 2>&1 \n" % locals()
    cmd += 'echo "%(cmd1)s" >> %(out)s \n'
    cmd += cmd1
    cmd += 'echo "*** RUNNING SCONS.COMPILE..." >> %(out)s \n'
    cmd2 = "./scons.compile -j %(procs)d >> %(out)s 2>&1 "% locals()
    cmd += 'echo "%(cmd2)s" >> %(out)s \n'
    cmd += cmd2
    #os.system(cmd % locals())
    global proc
    os.environ['JAVA_HOME'] = tabJava.getValue('JAVA_HOME')
    proc = Popen(cmd % locals(), shell=True)    
    text.fillTextArea()
    text.doRefresh(3)
    checkProcess()    
    setState(compiling=True)
    
def stopCompile(event=None):
    proc.terminate()
    setState(compiling=False)
    
btn.config(command=runCompile)
centerWindows(root)
root.deiconify()
root.mainloop()