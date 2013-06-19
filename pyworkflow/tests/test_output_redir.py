'''
Created on Jun 6, 2013

@author: antonio
'''
import unittest
import sys
import time


class TestRedirect(unittest.TestCase):


    def testRedirect(self):
        fOut = open('/home/antonio/Documents/Scipion/testFileTransfer/stdout.txt', 'a')
#         fErr = open('/home/antonio/Documents/Scipion/stderr.txt', 'a')
        sys.stdout = fOut
#         sys.stderr = fErr
        numberList = list(range(20))
        for num in numberList:
            print "WACKA " + str(num)
            #self.copyFile('/home/antonio/Documents/Scipion/testFileTransfer/stdout.txt', '/home/antonio/Documents/Scipion/testFileTransfer/stdoutCopy.txt')
            fOut.flush()
            time.sleep(3)
        fOut.close()
        
    def copyFile(self, file1, file2):
        f = open(file1)
        g = open(file2, 'a')
        linea = f.readline()
        while linea != "":
            g.write(linea)
            linea = f.readline()
        g.close()
        f.close()
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()