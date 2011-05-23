#!/usr/bin/env python

class XmippTestResult:
    def __init__(self, name, succeed=True, message=""):
           self.name = name
           self.succeed = succeed
           self.message = message
           
class XmippTestGroup:    
    def __init__(self, name):
        self.failed = 0
        self.tests = []
        self.name = name
    
    def addTest(self, test):
        self.tests.append(test)
        if not test.succeed:
            self.failed += 1
    
  
