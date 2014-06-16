#!/usr/bin/env python

import os
import logging
import unittest

from pyworkflow.utils.utils import getLineInFile, isInFile
from pyworkflow.utils.log import ScipionLogger, LOG_FILE
from pyworkflow.tests import BaseTest, setupTestOutput


#FIXME:Nacho
# Ok, Nacho and Airen, explain what you have to fix! :)

class TestLogs(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)        

    def testSimpleFileLog(self):
        import random
        logTestCode = random.randint(1, 100000)

        genLogFn = LOG_FILE
        log1 = logging.getLogger('pyworkflow.test.log.test_scipon_log')
        genInfoTest = 'General info [%d]' % logTestCode
        genDebugTest = 'General debug [%d]' % logTestCode
        genWarningTest = 'General warning [%d]' % logTestCode
        genErrorTest = 'General error [%d]' % logTestCode
        log1.info(genInfoTest)
        #log.debug(genDebugTest)
        log1.warning(genWarningTest)

        logFn = 'pyworkflow/tests/fileLog.log'
        log2 = ScipionLogger(logFn)
        fileInfoTest = 'File info [%d]' % logTestCode
        fileDebugTest = 'File debug [%d]' % logTestCode
        fileWarningTest = 'File warning [%d]' % logTestCode
        fileErrorTest = 'File error [%d]' % logTestCode
        log2.info(fileInfoTest)
        #log.debug(fileDebugTest)
        log2.warning(fileWarningTest)
        log3 = logging.getLogger('pyworkflow.tests.log')
        log3.error(genErrorTest)
        
        log4 = ScipionLogger(logFn)
        log4.error(fileErrorTest)
        
        # Check general logs
        lineGenInfoTest = getLineInFile(genInfoTest, genLogFn)
        lineGenWarningTest = getLineInFile(genWarningTest, genLogFn)
        lineGenErrorTest = getLineInFile(genErrorTest, genLogFn)
        
        isFileInfoTest = isInFile(fileInfoTest, genLogFn)
        isFileWarningTest = isInFile(fileWarningTest, genLogFn)
        isFileErrorTest = isInFile(fileErrorTest, genLogFn)     
        
        genLoggerChecked = True
        if lineGenInfoTest is None:
            print ('General info log failed!!!')
            genLoggerChecked = False
        if lineGenWarningTest is None:
            print ('General warning log failed!!!')
            genLoggerChecked = False
        if lineGenErrorTest is None:
            print ('General error log failed!!!')
            genLoggerChecked = False
        
        if not((lineGenInfoTest<lineGenWarningTest) & (lineGenWarningTest<lineGenErrorTest)):
            print ('General logs have an incorrect order!!!')
            genLoggerChecked = False
        
        if (isFileInfoTest | isFileWarningTest | isFileErrorTest):
            print ('File logs in general log!!!')
            genLoggerChecked = False
        
        # Check file logs
        lineFileInfoTest = getLineInFile(fileInfoTest, logFn)
        lineFileWarningTest = getLineInFile(fileWarningTest, logFn)
        lineFileErrorTest = getLineInFile(fileErrorTest, logFn)
        
        isGenInfoTest = isInFile(genInfoTest, logFn)
        isGenWarningTest = isInFile(genWarningTest, logFn)
        isGenErrorTest = isInFile(genErrorTest, logFn)    
        
        fileLoggerChecked = True
        if lineFileInfoTest is None:
            print ('File info log failed!!!')
            fileLoggerChecked = False
        if lineFileWarningTest is None:
            print ('File warning log failed!!!')
            fileLoggerChecked = False
        if lineFileErrorTest is None:
            print ('File error log failed!!!')
            fileLoggerChecked = False
        
        if not((lineFileInfoTest<lineFileWarningTest) & (lineFileWarningTest<lineFileErrorTest)):
            print ('File logs have an incorrect order!!!')
            fileLoggerChecked = False
        
        if (isGenInfoTest | isGenWarningTest | isGenErrorTest):
            print ('General logs in file log!!!')
            fileLoggerChecked = False 
        
        self.assertTrue(genLoggerChecked & fileLoggerChecked)

        print 'Removing log file %s' % logFn
        os.remove(logFn)
        
        
if __name__ == '__main__':
    unittest.main()
