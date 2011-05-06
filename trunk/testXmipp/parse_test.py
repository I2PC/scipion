#!/usr/bin/env python
import datetime, os, shutil
#create log file
from xml.sax.handler import ContentHandler
from xml.sax import make_parser
class TaskXMLHandler(ContentHandler):    
    #TODO start elemennt  processing must be centralized
    def __init__(self, _fnDir):
        fnDir = _fnDir
        lastProgram = ""
        progDict = {}
        MPI = False
        prerun = []
        testfile = []
        changeDirectory = False
        filesToCheck = []
        
    def startElement(self, name, attrs):
        print "name=",name
        if (name == "XMIPP_TESTS") :
            self.mySuccess = True
            self.myMessage = ""
        elif (name == "PROGRAM") :
            self.programName = attrs.get("name")
            self.mpi = attrs.get("mpi")
            self.progDict[self.programName] = []

        elif (name == "CASE") :
            self.argument = attrs.get("argument")
            self.changeDirectory = attrs.get("changeDirectory")
            self.prerun = []
            self.testfile = []
            
        elif (name == "PRERUN") :
            self.prerun.append(attrs.get("command"))
            
        elif (name == "TESTFILE") :
            self.prerun.append(attrs.get("filename"))

    def endElement(self, name):
        if (name == "XMIPP_TESTS") :
            self.funcionToDoTheJob()
        elif (name == "PROGRAM") :
            if(self.programName==""):
                self.reportError()
        elif (name == "CASE") :
            self.addcase()

    def addcase(self):
        self.progDict[self.programName].append(self.argument,
	                                       self.MPI,
					       self.prerun,
					       self.changeDirectory,
					       self.filesToCheck)

    def reportError(self):
	    self.myMessage += "Program name missing \n"
	    self.mySuccess = False

    def funcionToDoTheJob(self):
        import pprint
        pp = pprint.PrettyPrinter(indent=4,width=20)
        pp.pprint(dict)


print "0"

task = TaskXMLHandler("kkk")
saxparser = make_parser()
saxparser.setContentHandler(task)
testXml = 'test.xml'
if not os.path.exists(testXml):
    globalMessage += "\n cannot read file %s" % (testXml)
else:
    print "here"
    saxparser.parse(testXml)
 
