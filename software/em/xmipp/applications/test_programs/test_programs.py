#!/usr/bin/env xmipp_python
from datetime import datetime, timedelta
from os.path import exists, join, abspath
from pyworkflow.em.packages.xmipp3 import greenStr, warnStr, redStr, XmippScript
from xml.sax import make_parser
from xml.sax.handler import ContentHandler
from xmipp import *
import os
import shutil
import sys
import subprocess
import threading

class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print redStr('ERROR: timeout reached for this process')
            self.process.terminate()
            thread.join()
        self.process = None
    
    def terminate(self):
        if self.process != None:
            self.process.terminate()
            print redStr('Ctrl-c pressed, aborting this test')

class Tester(ContentHandler):    
    #TODO start elemennt  processing must be centralized
    def __init__(self, fnDir, continueFrom, programName, staPath, timeout, programPattern=None):
        self.fnDir = fnDir
        self._continue = continueFrom
        self._programName = programName
        self._statistics = staPath
        self.programPattern = programPattern
        self._timeout = timeout
        self._command = None
        if programName or programPattern:
            self.doit = False
        else:
            self.doit = True
        self.lastProgram = ""
        self.owner = ""
        self.progDict = {}
        self.mpi = False
        self.random = False
        self.prerun = []
        self.postrun = []
        self.testfile = []
        self.changeDirectory = False
        self.error = "init error\n"
        self.errorFlag = False
        self.warning = ""
        self.warningFlag = False
        
    def startElement(self, name, attrs):
        if (name == "XMIPP_TESTS") :
            self.mySuccess = True
            self.myMessage = ""
        elif (name == "PROGRAM") :
            self.programName = attrs.get("name")
            if attrs.has_key("mpi"):
                self.mpi = (attrs.get("mpi") == "TRUE")
            else:
                self.mpi = False
            if attrs.has_key("owner"):
                self.owner = attrs.get("owner")
            else:
                self.owner = "unassigned"

            self.progDict[self.programName] = []

        elif (name == "CASE") :
            self.arguments = attrs.get("arguments")
            self.changeDirectory = (attrs.get("changeDir") == "TRUE")
            self.prerun = []
            self.postrun = []
            self.testfile = []
            if(attrs.has_key("random")):
                self.random = (attrs.get("random") == "TRUE")
            else:
                self.random = False
            
        elif (name == "PRERUN") :
            self.prerun.append(attrs.get("command"))
        elif (name == "POSTRUN") :
            self.postrun.append(attrs.get("command"))    
        elif (name == "FILE") :
            self.testfile.append(attrs.get("filename"))

    def endElement(self, name):
        if (name == "XMIPP_TESTS") :
            #self.funcionToDoTheJob()
            print "xml parsing finished"
        elif (name == "PROGRAM") :
            if(self.programName == ""):
                self.reportError()
        elif (name == "CASE") :
            self.addcase()

    def matchPattern(self):
	if self.programPattern:
	    return self.programPattern in self.programName
	return False

    def addcase(self):
        if self.doit or self._programName==self.programName or self.matchPattern():
	    if self._continue:
	        self.doit=True
	    self.progDict[self.programName].append((self.arguments,
                                               self.mpi,
                                               self.random,
                                               self.prerun,
                                               self.postrun,
                                               self.changeDirectory,
                                               self.testfile,
                                               self.owner))

    def reportError(self):
        self.myMessage += "Program % missing \n" % (self.programName)
        self.mySuccess = False

    def funcionToDoTheJob(self):
        import pprint
        pp = pprint.PrettyPrinter(indent=4, width=20)
        pp.pprint(self.progDict)

    def runAllTests(self):
        for program in sorted(self.progDict):
            self.runProgramTests(program)
            
    def checkResult(self, testfiles, outDir, random):
        error = ""
        try:
            from xmipp import FileName, Image, compareTwoFiles, compareTwoMetadataFiles
            for file in testfiles:
                file = join(outDir, file)
                fileGoldStd = file.replace(self.fnDir, 'goldStandard')
                result = True
                if not exists(file):
                    error += file + " was NOT produced\n"
                else:
                    if (random):
                        self.warning += file + " was created using a random seed, check skipped\n"
                        self.warningFlag = True
                    else:
                        fnGoldStd = FileName(fileGoldStd)
                        if fnGoldStd.isImage():
                            im1 = Image(fileGoldStd)
                            im2 = Image(file)
                            result = im1.equal(im2, 0.001)
                        elif fnGoldStd.isMetaData():
                            result = compareTwoMetadataFiles(file, fileGoldStd)
                        else:
                            result = compareTwoFiles(file, fileGoldStd, 0)
                if not result:
                    error += " file '%s' and '%s' are NOT identical\n" % (file, fileGoldStd)
        except KeyboardInterrupt:
            raise
        except Exception, ex:
            error += " Error checking results: '%s'\n" % str(ex)
        return error
    
    def expandFormat(self, cmd):
        ''' Expand our tags %o, %p and %d with corresponding values '''
        cmd = cmd.replace("%o", self.outDir)
        cmd = cmd.replace("%p", self.program)
        cmd = cmd.replace("%d", self.fnDir)
        return cmd     
    
    def runCommands(self, cmdList, cmdType):
        pipe = '>'
        outDir = self.outDir
        for cmd in cmdList:
            if cmd != "":
                cmd = self.expandFormat(cmd)
                cmd = " %(cmd)s %(pipe)s %(outDir)s/%(cmdType)s_stdout.txt 2%(pipe)s %(outDir)s/%(cmdType)s_stderr.txt"\
                      % locals()
                print "    Running %s: " % cmdType, cmd
                command = Command(cmd)
                command.run(timeout=self._timeout)
                #os.system(cmd)
                #subprocess.call(cmd, shell=True)
                pipe = ">>"
           
    def runProgramTests(self, program):
        tests = self.progDict[program]
        n = len(tests)
        outPath = join(self.fnDir, program)
        outDir = outPath
        testName = ""
        testNo = 1
        
        md = MetaData(self._statistics)
        
        for test, mpi, random, preruns, postruns, changeDirectory, testfiles, owner in tests:
        #for test, mpi, random, preruns, postruns, changeDirectory, testfiles in tests:
            objid = md.addObject()
            if n > 1:
                outDir = outPath + "_%02d" % testNo
                testName = "(%d of %d)" % (testNo, n)
            #test num
            md.setValue(MDL_BLOCK_NUMBER, testNo , objid)
            #program name
            #print type(program)
            md.setValue(MDL_PROGRAM, str(program), objid)
            dtBegin = datetime.now()
            timeStr = str(dtBegin)
            #beginning time
            md.setValue(MDL_DATE, str(timeStr), objid)
            
            print "------------------------------------------------------------------------------------"
            print warnStr(">>> Running test %(testName)s of %(program)s (%(timeStr)s)" % locals())
            print "    Output dir: "
            print "       ", outDir
            print "    Statistics file: "
            print "       ", self._statistics
            print "    Timeout: "
            print "       %d seconds" % self._timeout
            
            if exists(outDir):
                shutil.rmtree(outDir)
            os.makedirs(outDir)
            self.outDir = outDir
            self.program = program
            test = self.expandFormat(test)            
            self.runCommands(preruns, 'prerun')
            
            if mpi:
                cmd = "mpirun -np 3 `which %s`" % program##DO NOT REPLACE SO WE CAN TEST MPI EASILY.replace("xmipp_", "xmipp_mpi_")
            else:
                cmd = program
            if changeDirectory:
                cmd = "cd %s ; " % outDir + cmd + " %s > stdout.txt 2> stderr.txt" % test
            else:
                cmd += " %s > %s/stdout.txt 2> %s/stderr.txt" % (test, outDir, outDir)
            print "    Command: "
            print "       ", greenStr(cmd)
                
            #run the test itself
            command = Command(cmd)
            self._command = command
            try:
                command.run(timeout=self._timeout)
            except KeyboardInterrupt:
                command.terminate()
            #result = os.system(cmd)
            #result = subprocess.call(cmd, shell=True) 

            self.runCommands(postruns, 'postrun')

            tdEnd = (datetime.now() - dtBegin).total_seconds()
            #elapsed time
            md.setValue(MDL_TIME, tdEnd, objid)
            
            print "    Elapsed time: %d seconds" % tdEnd
            error = self.checkResult(testfiles, outDir, random)
            if len(error):
                self.error += error
                self.errorFlag = True
                print redStr("ERRORS:\n" + error)
                md.setValue(MDL_ENABLED, -1, objid)
            else:
                md.setValue(MDL_ENABLED, 1, objid)
            testNo += 1
            md.setValue(MDL_USER, str(owner), objid)
        md.write(self._statistics)

    def listProgramTests(self, pattern=None):
        total = 0
        for program, tests in sorted(self.progDict.iteritems()):
            if pattern is None or pattern in program:
                n = len(tests)
                total = total + n
                print "(%(n)d) %(program)s" % locals()
        print "Total: %d" % total        

    def terminate(self):
        self._command.terminate()

class ScriptProgramsTester(XmippScript):    
    def __init__(self):
        XmippScript.__init__(self,True)
        self._identities = { 'all' : 'xmipp@cnb.csic.es' ,'joton' : 'joton@cnb.csic.es', 'delarosatrevin' : 'jmdelarosa@cnb.csic.es', 'nachofoche' : 'ifoche@cnb.csic.es', 'rmarabini' : 'roberto@cnb.csic.es', 'jvargas' : 'jvargas@cnb.csic.es', 'coss' : 'coss@cnb.csic.es', 'vahid' : 'vabrishami@cnb.csic.es' }
        self._statistics = ''
        self._tester = None
        self._summaryMessage = "Xmipp_goldstandard_report"
        self._outputDirectory  = ""
        self._outputTarFile = join(self._outputDirectory, "execution_output.tar.gz")
        
    def defineParams(self):
        self.addUsageLine('Run program tests')
        ## params
        self.addParamsLine(' [-d <directory=tmp>]   : Output Directory')
        self.addParamsLine('   alias --outputdir;')
        self.addParamsLine(' [-p <program_name>]     : Execute this program')
        self.addParamsLine('   alias --program_name;')
        self.addParamsLine(' [-c]                   : process all programs after programname')
        self.addParamsLine('   alias --continue;')
        self.addParamsLine("   requires --program_name;");
        self.addParamsLine(' [-t <pattern>]     : Execute programs containing this pattern')
        self.addParamsLine('   alias --pattern;')
        self.addParamsLine(' [-l ]     : List all programs tests')
        self.addParamsLine('   alias --list;')
        self.addParamsLine(" [-s <statistics=\"statistics.xmd\">]     : Metadata with statistics")
        self.addParamsLine('   alias --statistics;')
        self.addParamsLine(" [-m]     : send mail report by email")
        self.addParamsLine('   alias --mail;')
        self.addParamsLine(" [-e <edit_subject>]     : subject for the sent mail")
        self.addParamsLine('   alias --edit_subject;')
        self.addParamsLine(' [-o]                   : send ')
        self.addParamsLine('   alias --mail_only_errors;')
        self.addParamsLine("   requires --mail;");
        self.addParamsLine(' [-k]                   : send ')
        self.addParamsLine('   alias --send_output;')
        self.addParamsLine("   requires --mail;");
        self.addParamsLine(' [-b]                   : send ')
        self.addParamsLine('   alias --broadcast;')
        self.addParamsLine("   requires --mail;");
        self.addParamsLine(" [-g]     : List of buggers in a Metadata")
        self.addParamsLine('   alias --list_buggers;')
        self.addParamsLine(" [-f]     : List of failed tests in a Metadata")
        self.addParamsLine('   alias --list_failed_tests;')
        self.addParamsLine(" [-z <timeout=300>]     : Timeout in seconds for each test")
        self.addParamsLine('   alias --timeout;')
        ## examples
        self.addExampleLine('   ./test_programs.py -d directory -p xmipp_metadata_utilities -c')

    def terminate(self):
        self._tester.terminate()

    def run(self):
        import os
        #read command line
        from numpy  import array, dot
        self._outputDirectory = _outputDirectory  = self.getParam('--outputdir')
        _programName = ""
        _summaryMessage = "Xmipp_goldstandard_report"
        _continue = False
        self._statistics = self.getParam('--statistics')
        _mail = self.checkParam("--mail")

        if self.checkParam('--edit_subject'):
            _summaryMessage = self._summaryMessage = self.getParam('--edit_subject')
        
        #if the user wants to list the failed tests in a metadata, we execute the proper function
        if self.checkParam('--list_buggers'):
            print "*" * 50
            print "buggers present on %s file:" % self._statistics
            print "_" * 50
            buggers = set(self.getBuggers(self._statistics))
            if len(buggers) != 0:
                for bugger in buggers:
                    bugmail = self._identities.get(bugger)
                    if bugmail == None:
                        bugmail = "Mail address not available"
                    print "%s - %s" % (bugger, bugmail)
                print "*" * 50
                if self.checkParam("--mail"):
                    self.mailResults(sendOnlyErrors=self.checkParam("--mail_only_errors"), broadcast=self.checkParam("--broadcast"), sendOutput=self.checkParam("--send_output"))
            else:
                print "None! You're lucky"
                print "*" * 50
                if self.checkParam("--mail"):
                    print "No mail sent 'cause no failed test present"
            exit(0)

        if self.checkParam('--list_failed_tests'):
            print "*" * 50
            print "failed tests present on %s file:" % self._statistics
            print "_" * 50
            failed_tests = set(self.getFailedTests(self._statistics))
            if len(failed_tests) != 0:
                for test in failed_tests:
                    print "%s" % test
                print "*" * 50
                if self.checkParam("--mail"):
                    self.mailResults(sendOnlyErrors=self.checkParam("--mail_only_errors"), broadcast=self.checkParam("--broadcast"), sendOutput=self.checkParam("--send_output"))
            else:
                print "None! You're lucky"
                print "*" * 50
                if self.checkParam("--mail"):
                    print "No mail sent 'cause no failed test present"
            exit(0)

        #otherwise, the statistics file will be looked in the output directory and will be the file to place the tests results
        self._statistics = join(_outputDirectory, self._statistics)
    
        if self.checkParam('--program_name'):
            _programName = self.getParam('--program_name')
            _continue = self.checkParam('--continue')
        

        pattern = None
        if self.checkParam('--pattern'):
            pattern = self.getParam('--pattern')
            
        tmpLink = 'tmpLink'
        if not exists(_outputDirectory):
            os.makedirs(_outputDirectory)
        if exists(tmpLink):
            os.remove(tmpLink)
        os.symlink(_outputDirectory, tmpLink)

        # Create tester
        self._tester = Tester(tmpLink,_continue, _programName, self._statistics, int(self.getParam('--timeout')), pattern)
        saxparser = make_parser()
        saxparser.setContentHandler(self._tester)
        testXml = 'test.xml'
        if not exists(testXml):
            globalMessage += "\n cannot read file %s" % (testXml)
        else:
            saxparser.parse(testXml)
    
        if self.checkParam("--list"):
            self._tester.listProgramTests(pattern)
            return
        
    #	programs = ""
            # Remove the output directory if it is not goldStandard
        if _outputDirectory != 'goldStandard' and _outputDirectory != 'goldStandard/':
            if exists(_outputDirectory):
                shutil.rmtree(_outputDirectory)
                os.makedirs(_outputDirectory)

        md = MetaData()
        md.write(self._statistics)

        self._tester.runAllTests()
        programs = '\n'.join(self._tester.progDict.keys())

#        if self._tester.errorFlag:
#            print redStr("ALL ERRORS:\n" + self._tester.error)
#        if self._tester.warningFlag:
#            print warnStr("ALL WARNINGS:\n" + self._tester.warning)

        
        #send the results by email
        if self._tester.errorFlag:
            if _mail:
                self.mailResults(sendOnlyErrors=self.checkParam("--mail_only_errors"), broadcast=self.checkParam("--broadcast"), sendOutput=self.checkParam("--send_output"))
        elif self.checkParam("--mail"):
            print "No mail sent 'cause no failed test present"

    def mailResults(self, sendOnlyErrors=False, broadcast=False, sendOutput=False):
        ''' Send mail to the owner of each test with the report of the tests. 
        With sendOnlyErrors you can select only send errors exists, and with broadcast
        you can select send to everyone in the identities dictionary. Otherwise, each
        owner only receive his/her tests report '''
        lib_path = abspath('./script')
        sys.path.append(lib_path)
        import mail
        import config
        
        summaryMessage = self._summaryMessage 
        globalMessage = "Attached you will find the metadata with the report"
        toaddrs = ""
        attach_path = ""
        attachments = []

        path = self._statistics.split('/')
        extension = path[-1].split('.')[-1]
        last = path[-1].split('.')[0:len(path[-1].split('.'))-1]
        laststr = ""

        laststr = ".".join(last)

        if extension == laststr:
            extension = ""
        previous = path[0:len(path)-1]
        
        attach_path = "/".join(previous)

        print attach_path

        if sendOnlyErrors:
            md = self.extractMDWithError(self._statistics)
            if len(previous):
                attach_path += "/"
            attach_path += laststr + "_onlyerrors." + extension
            md.write(attach_path)
        else:
            if len(previous):
                attach_path += "/"
            attach_path += laststr + "." + extension
        
        if broadcast:
            toaddrs = self._identities.get('all')
        else:
            buggers = set(self.getBuggers(self._statistics))
            if len(buggers) != 0:
                for bugger in buggers:
                    toaddrs += "%s, " % self._identities.get(bugger)
        print "sending %s attached in an email..." % attach_path
        attachments.append(attach_path)
        if sendOutput:
            self.packOutput(self._outputTarFile)
            attachments.append(self._outputTarFile)
            print "sending also %s attached in the email..." % self._outputTarFile
            
        mail.mail(toaddrs,config.fromaddr,summaryMessage,globalMessage,attachments)
        print "...done"

    def packOutput(self, outputTarFile):
        '''Does the packaging of the tests output to be sent by email'''
        import tarfile
        import glob
        if exists(self._outputTarFile):
            os.remove(self._outputTarFile)
        with tarfile.open(self._outputTarFile, "w:gz") as tar:
            for name in glob.glob(self._outputDirectory  + '/xmipp_*/*.txt'):
                tar.add(name)
            tar.close()
        print "tar file done"
        

    def extractMDWithError(self, statistics):
        ''' Return a MD with only the failing tests in the given MD'''
        md = MetaData(statistics)
        md2 = MetaData()
        md2.importObjects(md, MDValueEQ(MDL_ENABLED, -1))
        return md2

    def getBuggers(self, statistics):
        ''' Return the list of owners for the failing tests in the given MD'''
        md2 = self.extractMDWithError(statistics)
        buggers = md2.getColumnValues(MDL_USER)
        return buggers
    
    def getFailedTests(self,statistics):
        ''' Return the list of owners for the failing tests in the given MD'''
        md2 = self.extractMDWithError(statistics)
        failed_tests = md2.getColumnValues(MDL_PROGRAM)
        return failed_tests

if __name__ == '__main__':
    script = ScriptProgramsTester()
    try:
        script.tryRun()
    except KeyboardInterrupt:
        script.terminate()

