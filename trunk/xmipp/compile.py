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

import os
from subprocess import Popen

################ Classes for a Console based configuration ######################
class ConsoleOptionsTab():
    def __init__(self, master, text, **opts):
        self.name = text
        self.optionsDict = {}
        self.optionsValue = {}
        
    def addOption(self, name, comment, default='', cond=None):
        self.optionsDict[name] = (comment, default, cond)
        self.optionsValue[name] = default
        
    def setValue(self, name, value):
        self.optionsValue[name] = value
        
    def getValue(self, name):
        return self.optionsValue[name]
        
    def getConfigOptions(self):
        optStr = ""
        for key, comment, default, cond in self.optionsDict.iteritems():
            value = self.optionsValue[key]
            if (cond is None or self.getValue(cond) == 'yes') and default != value: 
                optStr += ' %s="%s"' % (k, value)
        return optStr
    

class ConsoleConfigNotebook():
    def __init__(self):
        self.tabs = {}
        
    def addTab(self, text):
        tab = ConsoleOptionsTab(self, text)
        self.tabs[text] =  tab
        return tab
    
    def getConfigOptions(self):
        '''return string with configs settings '''
        return ' '.join([t.getConfigOptions() for t in self.tabs.values()])
    
    def setValue(self, tab, option, value):
        self.tabs[tab].setValue(option, value)
        
    def getValue(self, tab, option):
        return self.tabs[tab].getValue(option)
    
    def notifyCompile(self, process):
        pass


################### Helper functions ####################
def detectJava():
    from protlib_utils import which
    java_home = None
    javac = which('javac')
    if javac:
        from os.path import dirname
        from protlib_filesystem import findRealFile
        java_home = dirname(dirname(findRealFile(javac)))
    return java_home
    
def detectMpi():
    from protlib_filesystem import findFilePath
    inc_dirs = ['/usr/include', 'usr/local/include']
    lib_dirs = ['/usr/lib64', 'usr/lib']
    inc_mpi = findFilePath('mpi.h', *(inc_dirs + lib_dirs))
    lib_mpi = findFilePath('libmpi.so', *lib_dirs)
    
    return (inc_mpi, lib_mpi)

def detectQt():
    pass

def addTabs(nb):
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
        
def runCompile(notebook, numberOfCpu):
    out = OUTPUT
    procs = numberOfCpu
    opts = notebook.getConfigOptions()
    cmd = 'echo "*** RUNNING SCONS.CONFIGURE..." > %(out)s \n'
    cmd1 = "./scons.configure %(opts)s >> %(out)s 2>&1 \n" % locals()
    cmd += 'echo "%(cmd1)s" >> %(out)s \n'
    cmd += cmd1
    cmd += 'echo "*** RUNNING SCONS.COMPILE..." >> %(out)s \n'
    cmd2 = "./scons.compile -j %(numberOfCpu)d >> %(out)s 2>&1 "% locals()
    cmd += 'echo "%(cmd2)s" >> %(out)s \n'
    cmd += cmd2
    os.environ['JAVA_HOME'] = notebook.getValue('Java & QT', 'JAVA_HOME')
    proc = Popen(cmd % locals(), shell=True)    
    notebook.notifyCompile(proc)    
    
def stopCompile(notebook, proc):
    proc.terminate()
    notebook.notifyStopCompile("Compilation aborted by user")


# Output file
OUTPUT = 'build/scons_output.log'
if os.path.exists(OUTPUT):
    os.remove(OUTPUT)
# Number of CPU
if os.environ.has_key('NUMBER_OF_CPU'):
    numberOfCpu = int(os.environ['NUMBER_OF_CPU'])
else: 
    numberOfCpu = 2
    
# Display Tkinter GUI
GUI = True

if GUI:
    from compile_gui import createGUINotebook
    nb = createGUINotebook(OUTPUT, numberOfCpu, addTabs, runCompile, stopCompile)
else:
    nb = ConsoleConfigNotebook()
    
