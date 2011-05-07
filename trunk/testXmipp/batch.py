#!/usr/bin/env python
import datetime, os, shutil,sys
#create log file
from xml.sax.handler import ContentHandler
from xml.sax import make_parser
class Tester(ContentHandler):    
    #TODO start elemennt  processing must be centralized
    def __init__(self, _fnDir):
        self.fnDir = _fnDir
        self.lastProgram = ""
        self.progDict = {}
        self.mpi = False
        self.prerun = []
        self.testfile = []
        self.changeDirectory = False
        self.error=""
        self.errorFlag=False
        
    def startElement(self, name, attrs):
        if (name == "XMIPP_TESTS") :
            self.mySuccess = True
            self.myMessage = ""
        elif (name == "PROGRAM") :
            self.programName = attrs.get("name")
            self.mpi = (attrs.get("mpi") == "TRUE ")
            self.progDict[self.programName] = []

        elif (name == "CASE") :
            self.arguments = attrs.get("arguments")
            self.changeDirectory = (attrs.get("changeDir") == "TRUE ")
            self.prerun = []
            self.testfile = []
            
        elif (name == "PRERUN") :
            self.prerun.append(attrs.get("command"))
            
        elif (name == "FILE") :
            self.testfile.append(attrs.get("filename"))

    def endElement(self, name):
        if (name == "XMIPP_TESTS") :
            #self.funcionToDoTheJob()
            print "xml parsing finished"
        elif (name == "PROGRAM") :
            if(self.programName==""):
                self.reportError()
        elif (name == "CASE") :
            self.addcase()

    def addcase(self):
        self.progDict[self.programName].append((self.arguments,
	                                       self.mpi,
					       self.prerun,
					       self.changeDirectory,
					       self.testfile))

    def reportError(self):
	    self.myMessage += "Program name missing \n"
	    self.mySuccess = False

    def funcionToDoTheJob(self):
        import pprint
        pp = pprint.PrettyPrinter(indent=4,width=20)
        pp.pprint(self.progDict)

    def runAllTests(self):
        for program in self.progDict.keys():
            self.runProgramTests(program)
        
    def checkResult(self,testfiles,outDir):
        import xmipp
        for file in testfiles:
            file = outDir + '/' + file
            if not os.path.exists(file):
                self.error += file, "was NOT produced"
                self.errorFlag=True
            result = xmipp.compareTwoFiles(file,file.replace(self.fnDir,'goldStandard'),0)
            if not result:
                self.error += file + " and " + file.replace(self.fnDir,'goldStandard') + " are NOT identical"
                self.errorFlag=True
           
    def runProgramTests(self, program):
        tests = self.progDict[program]
        n = len(tests)
        outPath = os.path.join(self.fnDir, program)
        outDir = outPath
        testName = ""

        testNo = 1
        for test, mpi, preruns, changeDirectory,testfiles in tests:
            if n > 1:
                outDir = outPath + "_%02d" % testNo
                testName = "(%d of %d)" % (testNo, n)
            print "------------------------------------------------------------------------------------"
            print ">>> Running test", testName, "of", program
            print "    Output dir: "
            print "       ", outDir
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            test = test.replace("%o", outDir)
            test = test.replace("%p", program)
            test = test.replace("%d", self.fnDir)
            pipe=">"
            for prerun in preruns:
                if prerun != "":
                    prerun = prerun.replace("%o", outDir)
                    prerun = prerun.replace("%p", program)
                    prerun = prerun.replace("%d", self.fnDir)
                    cmd = " %s %s %s/prerun_stdout.txt 2%s %s/prerun_stderr.txt" %\
                        (prerun, pipe, outDir, pipe, outDir)
                    print "    Running prerun: ", cmd
                    os.system(cmd)
                    pipe=">>"
            if mpi:
                cmd = "mpirun -np 3 `which %s`" % program.replace("xmipp_", "xmipp_mpi_")
            else:
                cmd = program
            if changeDirectory:
                cmd = "cd %s ; " % outDir + cmd + " %s > stdout.txt 2> stderr.txt" % test
            else:
                cmd += " %s > %s/stdout.txt 2> %s/stderr.txt" % (test, outDir, outDir)
            print "    Command: "
            print "       ", cmd
            result = os.system(cmd)
            self.checkResult(testfiles,outDir)
            #print "Result:", result
            testNo += 1
        
if __name__ == '__main__':

    argc = len(sys.argv)

    if argc < 2 or argc > 3:
        print "Usage: ./batch.py <directory> [program]"
        sys.exit()

    fnDir = sys.argv[1]
    tmpLink = 'tmpLink'
    if not os.path.exists(fnDir):
        os.makedirs(fnDir)
    if os.path.exists(tmpLink):
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

    if argc > 2:
        program = sys.argv[2]
        tester.runProgramTests(program)
        #if not os.path.exists(fnDir):
        #    os.makedirs(fnDir)
    else:
        # Remove the output directory if it is not goldStandard
        if fnDir != 'goldStandard' and fnDir != 'goldStandard/':
            if os.path.exists(fnDir):
                shutil.rmtree(fnDir)
            os.makedirs(fnDir)
        tester.runAllTests()
    if (tester.errorFlag):
        print "ERROR:"
        print "\n\n",tester.error
 
