#!/usr/bin/env python
import datetime, os,shutil,config
import mail,glob
#create log file
from xml.sax.handler import ContentHandler
from xml.sax import make_parser
class TaskXMLHandler(ContentHandler):    
    #TODO start elemennt  processing must be centralized
    nFailures= 0
    inFailures = False
    myMessage=""
    mySuccess=True
    def startElement(self, name, attrs): 
        if (name == "testsuite") :
            self.mySuccess=True
            self.myMessage = ""
        elif (name == "testcase") :
 	    self.inFailures = False   
            self.name      = attrs.get("name")
            self.classname = attrs.get("classname")
        elif (name == "failure") :
            self.message = attrs.get("message")
	    self.inFailures = True   

    def endElement(self, name):
        if (name == "testcase" and self.inFailures) :
            self.reportError()
	elif (name == "testcase"):
	    self.reportOk()

    def reportError(self):
	self.myMessage += "Test %s in testgroup %s FAILED with error %s\n"%(self.name,self.classname,self.message)
	self.mySuccess = False
	
    def reportOk(self):
	self.myMessage += "Test %s in testgroup %s was OK\n"%(self.name,self.classname)

def main(filename):
    """ Parse xml test files """
    from  config import XMIPP_TEST

    XMIPP_OUTPUT=XMIPP_TEST + '/OUTPUT'
    task = TaskXMLHandler()
    saxparser = make_parser()
    saxparser.setContentHandler(task)
    os.chdir(XMIPP_HOME+'/bin')
    testNames = glob.glob('xmipp_test_*')
    globalMessage=""
    success=True
    for testName in testNames:
        globalMessage += "\n"
        filename = XMIPP_OUTPUT+'/'+testName+'.xml'
        if not os.path.exists(filename):
	    globalMessage += "\n FAILED generation file %s by test %s"%(filename,testName)
	    success = False
	else:
	    saxparser.parse(filename)
	    globalMessage += task.myMessage
	    if not task.mySuccess:
	        success = False
    if success:
       summaryMessage='XMIPP compilation was OK'
    else:
       summaryMessage='XMIPP compilation FAILED'
    mail.mail(config.toaddrs,config.fromaddr,summaryMessage,globalMessage)
 
