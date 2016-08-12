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

from pyworkflow.tests import *
import pyworkflow.utils as pwutils

import xmipp



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
    
    _labels = [WEEKLY]
               
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
                
    def runCase(self, args, mpi=0, changeDir=False, 
                preruns=None, postruns=None, validate=None,
                outputs=None, random=False):
        # Retrieve the correct case number from the test name id
        # We asumme here that 'test_caseXXX' should be in the name
        caseId = unittest.TestCase.id(self)
        if not 'test_case' in caseId:
            raise Exception("'test_case' string should be in the test function name followed by a number")
        _counter = int(caseId.split('test_case')[1])

        self._testDir = self.dataset.getPath()
        self.outputDir = os.path.join('tmpLink', '%s_%02d' % (self.program, _counter))
        self.outputDirAbs = os.path.join(self._testDir, self.outputDir)
        self.goldDir = os.path.join(self._testDir, 'gold', '%s_%02d' % (self.program, _counter))
        
        # Clean and create the program output folder if not exists
        pwutils.cleanPath(self.outputDirAbs)
        pwutils.makePath(self.outputDirAbs)
        
        # Change to tests root folder (self._testDir)
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
            self._checkOutputs(outputs,random)
            
        if validate:
            validate()
            
        os.chdir(cwd)
        
    def _checkOutputs(self, outputs, random=False):
        """ Check that all output files are produced
        and are equivalent to the ones in goldStandard folder.
        """
        for out in outputs:
            outFile = os.path.join(self._testDir, self.outputDir, out)
            fileGoldStd = os.path.join(self.goldDir, out)
            
            # Check the expect output file was produced
            msg = "Missing expected output file:\n  output: %s" % outFile
            self.assertTrue(os.path.exists(outFile), msg)
            
            if random:
                print "WARNING: %s was created using a random seed, check skipped\n " % outFile
            else:
                fnGoldStd = xmipp.FileName(fileGoldStd)
                if fnGoldStd.isImage():
                    im1 = xmipp.Image(fileGoldStd)
                    im2 = xmipp.Image(outFile)
                    msg = "Images are not equal:\n  output: %s\n  gold: %s" % (outFile, fileGoldStd)
                    self.assertTrue(im1.equal(im2, 0.001), msg)
                elif fnGoldStd.isMetaData():
                    msg = "MetaDatas are not equal:\n  output: %s\n  gold: %s" % (outFile, fileGoldStd)
                    self.assertTrue(xmipp.compareTwoMetadataFiles(outFile, fileGoldStd), msg)
                else:
                    msg = "Files are not equal:\n  output: %s\n  gold: %s" % (outFile, fileGoldStd)
                    self.assertTrue(xmipp.compareTwoFiles(outFile, fileGoldStd, 0), msg)
