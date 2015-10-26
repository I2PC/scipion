# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

import os
import unittest
import threading
import subprocess

import pyworkflow.utils as pwutils



class Command(object):
    def __init__(self, cmd, env=None):
        self.cmd = cmd
        self.process = None
        self.env = env

    def run(self, timeout):
        def target():
            self.process = subprocess.Popen(self.cmd, shell=True, env=self.env)
            self.process.communicate()

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print pwutils.red('ERROR: timeout reached for this process')
            self.process.terminate()
            thread.join()
        self.process = None
    
    def terminate(self):
        if self.process != None:
            self.process.terminate()
            print pwutils.red('Ctrl-c pressed, aborting this test')
            
            
class ProgramTest(unittest.TestCase):
    _testDir = None
    _environ = None
    _timeout = 300
               
    @classmethod
    def setTestDir(cls, newTestDir):
        cls._testDir = newTestDir
    
    @classmethod
    def setEnviron(cls, newEnviron):
        cls._environ = newEnviron
        
    @classmethod
    def setTimeOut(cls, newTimeOut):
        cls._timeout = newTimeOut
    
    def _parseArgs(self, args):
        ''' Expand our tags %o, %p and %d with corresponding values '''
        args = args.replace("%o", self.outputDir)
        args = args.replace("%p", self.program)
        #args = args.replace("%d", self.fnDir)
        return args     
    
    def _runCommands(self, cmdList, cmdType):
        """ Run several commands.
        Params:
            cmdList: the list of commands to execute. 
            cmdType: either 'preruns' or 'postruns'
        """
        pipe = '>'
        outDir = self.outputDir
        for cmd in cmdList:
            if cmd:
                cmd = self._parseArgs(cmd)
                cmd = " %(cmd)s %(pipe)s %(outDir)s/%(cmdType)s_stdout.txt 2%(pipe)s %(outDir)s/%(cmdType)s_stderr.txt" % locals()
                print "    Running %s: " % cmdType, pwutils.cyan(cmd)
                command = Command(cmd, env=self.env)
                command.run(timeout=self._timeout)
                pipe = ">>"
                
    def runCase(self, args, mpi=0, changeDir=False, preruns=None, postruns=None, outputs=None):
        self._testDir = os.path.join(os.environ['SCIPION_TESTS'], 'testXmipp')
        self._counter += 1
        self.outputDir = os.path.join(self._testDir, 'test', '%s_%02d' % (self.program, self._counter))
        pwutils.cleanPath(self.outputDir)
        pwutils.makePath(self.outputDir)
        
        cwd = os.getcwd()
        os.chdir(self._testDir)
        
        if preruns:
            self._runCommands(preruns, 'preruns')
            
        if mpi:
            cmd = "mpirun -np %d `which %s`" % (mpi, self.program)
        else:
            cmd = self.program
        
        args = self._parseArgs(args)
        
        if changeDir:
            cmd = "cd %s ; %s %s > stdout.txt 2> stderr.txt" % (self.outputDir, cmd, args)
        else:
            cmd = "%s %s > %s/stdout.txt 2> %s/stderr.txt" % (cmd, args, self.outputDir, self.outputDir)
        print "    Command: "
        print "       ", pwutils.green(cmd)
            
        #run the test itself
        command = Command(cmd, env=self.env)
        self._command = command
        try:
            command.run(timeout=self._timeout)
        except KeyboardInterrupt:
            command.terminate()
        
        if postruns:
            self._runCommands(postruns, 'postruns')
            
        if outputs:
            for outFile in outputs:
                self.assertTrue(os.path.exists(os.path.join(self.outputDir, outFile)))
            
        os.chdir(cwd)
        
        