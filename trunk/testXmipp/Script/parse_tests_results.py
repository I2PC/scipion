#!/usr/bin/env python
import datetime, os,shutil,config
import glob
#create log file
from test_result import *
from xml.sax.handler import ContentHandler
from xml.sax import make_parser

class TaskXMLHandler(ContentHandler):    
    #TODO start elemennt  processing must be centralized
    def __init__(self):
        self.nFailures= 0
        self.inFailures = False
        self.myMessage=""
        self.mySuccess=True
        self.testsResults = []
        self.test = None
        self.group = None
    def startElement(self, name, attrs): 
        if (name == "testsuite") :
            self.mySuccess=True
            self.myMessage = ""
            name = attrs.get("name")
            self.group = XmippTestGroup(name)
            self.testsResults.append(self.group)
            
        elif (name == "testcase") :
            name      = attrs.get("name")
            classname = attrs.get("classname")
            self.test = XmippTestResult(name)
        elif (name == "failure") :
            self.test.message = attrs.get("message")
            self.test.succeed = False  
            
    def endElement(self, name):
        if (name == "testcase"):
            self.group.addTest(self.test)


def parseResults():
    """ Parse xml test files """
    from  config import XMIPP_TEST
    from  config import XMIPP_HOME

    XMIPP_OUTPUT = os.path.join(XMIPP_TEST, 'OUTPUT')
    task = TaskXMLHandler()
    saxparser = make_parser()
    saxparser.setContentHandler(task)
    os.chdir(os.path.join(XMIPP_HOME, 'bin'))
    testNames = glob.glob('xmipp_test_*')
    
    
    for testName in testNames:
        filename = os.path.join(XMIPP_OUTPUT, testName+'.xml')
        saxparser.parse(filename)
        
    return task.testsResults


if __name__ == '__main__':
    parseResults()
 
