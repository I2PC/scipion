#!/usr/bin/env xmipp_python
import datetime, os, shutil, sys
#create log file
from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from protlib_xmipp import greenStr, warnStr

class Tester(ContentHandler):    
    #TODO start elemennt  processing must be centralized
    def __init__(self, _fnDir):
        self.fnDir = _fnDir
        self.lastProgram = ""
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
            if(attrs.has_key("mpi")):
                self.mpi = (attrs.get("mpi") == "TRUE")
            else:
                self.mpi = False
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

    def addcase(self):

                           
        self.progDict[self.programName].append((self.arguments,
	                                       self.mpi,
                                           self.random,
					       self.prerun,
                           self.postrun,
					       self.changeDirectory,
					       self.testfile))

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
        from xmipp import FileName, Image, compareTwoFiles, compareTwoMetadataFiles
        for file in testfiles:
            file = os.path.join(outDir, file)
            fileGoldStd = file.replace(self.fnDir, 'goldStandard')
            result = True
            if not os.path.exists(file):
                self.error += file + " was NOT produced\n"
                self.errorFlag = True
            else:
                if (random):
                    self.warning += file + " was created using a random seed, check skipped\n"
                    self.warningFlag = True
                else:
                    print "comparing '%s' and '%s'" % (file, fileGoldStd)
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
                self.error += " file '%s' and '%s' are NOT identical\n" % (file, fileGoldStd)
                self.errorFlag = True
    
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
               cmd = " %(cmd)s %(pipe)s %(outDir)s/%(cmdType)s_stdout.txt 2%(pipe)s %(outDir)s/%(cmdType)s_stderr.txt" % locals()
               print "    Running %s: " % cmdType, cmd
               os.system(cmd)
               pipe = ">>"
           
    def runProgramTests(self, program):
        tests = self.progDict[program]
        n = len(tests)
        outPath = os.path.join(self.fnDir, program)
        outDir = outPath
        testName = ""
        testNo = 1
        for test, mpi, random, preruns, postruns, changeDirectory, testfiles in tests:
            if n > 1:
                outDir = outPath + "_%02d" % testNo
                testName = "(%d of %d)" % (testNo, n)
            print "------------------------------------------------------------------------------------"
            print warnStr(">>> Running test %(testName)s of %(program)s" % locals())
            print "    Output dir: "
            print "       ", outDir
            if not os.path.exists(outDir):
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
            result = os.system(cmd)
            self.runCommands(postruns, 'postrun')
            self.checkResult(testfiles, outDir, random)
            #print "Result:", result
            testNo += 1
        
if __name__ == '__main__':

    argc = len(sys.argv)

    if argc < 2:
        print "Usage: ./batch.py <directory> [programs]"
        sys.exit()

    fnDir = sys.argv[1]
    tmpLink = 'tmpLink'
    if not os.path.exists(fnDir):
        os.makedirs(fnDir)
    if os.path.lexists(tmpLink):
        os.remove(tmpLink)
    os.symlink(fnDir, tmpLink)
    # Create tester
    tester = Tester(tmpLink)
    saxparser = make_parser()
    saxparser.setContentHandler(tester)
    testXml = 'test.xml'
    if not os.path.exists(testXml):
        globalMessage += "\n cannot read file %s" % (testXml)
    else:
        saxparser.parse(testXml)

    programs = ""
    if argc > 2:
        programs = sys.argv[2:]
        for program in programs:
            tester.runProgramTests(program)
            #programs = program
        #if not os.path.exists(fnDir):
        #    os.makedirs(fnDir)
    else:
        # Remove the output directory if it is not goldStandard
        if fnDir != 'goldStandard' and fnDir != 'goldStandard/':
            if os.path.exists(fnDir):
                shutil.rmtree(fnDir)
            os.makedirs(fnDir)
        tester.runAllTests()
        for program in tester.progDict.keys():
            programs += program + "\n"

    if (tester.errorFlag):
        print "ERROR:"
        print "", tester.error
    if (tester.warningFlag):
        print "WARNING:"
        print "\n\n", tester.warning
    import os, sys
    lib_path = os.path.abspath('./Script')
    sys.path.append(lib_path)
    import config
    import mail

    globalMessage = ""
    if tester.errorFlag:
       summaryMessage = 'XMIPP goldstandard FAILED'
    else:
       summaryMessage = 'XMIPP goldstandard is OK'
       
    if  tester.errorFlag:
       globalMessage += "ERROR:\n" + tester.error
    if  tester.warningFlag:
       globalMessage += "WARNINGS:\n" + tester.warning
    globalMessage += "\nPrograms tested:\n" + str(programs)
    #mail.mail(config.toaddrs,config.fromaddr,summaryMessage,globalMessage)

 
