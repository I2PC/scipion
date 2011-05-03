#!/usr/bin/env python
import datetime, os,shutil,config
import mail
#create log file
from xml.sax.handler import ContentHandler
from xml.sax import make_parser
class TaskXMLHandler(ContentHandler):    
    #TODO start elemennt  processing must be centralized
    nFailures= 0
    inFailures = False
    
    def startElement(self, name, attrs): 
        if (name == "testcase") :
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
        #print 'testName',self.name
	#print 'testGroup', self.classname
	#print 'message', self.message
	message =self.classname + "\n" +self.name + "\n" +  self.message
	mail.mail(config.toaddrs,config.fromaddr,"xmipp tests failed",message)
	
    def reportOk(self):
	message ="xmipp compilation and test OK"
	mail.mail(config.toaddrs,config.fromaddr,"xmipp compilation and test OK",message)

def main(filename):
    """ Parse xml test files """
    from  config import XMIPP_TEST

    XMIPP_OUTPUT=XMIPP_TEST + '/OUTPUT'
    task = TaskXMLHandler()
    saxparser = make_parser()
    saxparser.setContentHandler(task)
    import glob
    myList= glob.glob(XMIPP_OUTPUT+'/*xml')
    if len(myList)==0:
	message ="No test data was generated. See xmipp@xmipp.cnb.csic.es:" + filename + " for detais"
	mail.mail(config.toaddrs,config.fromaddr,"xmipp compilation failed ",message)
    else:
        for filename in myList: 
	    saxparser.parse(filename)

