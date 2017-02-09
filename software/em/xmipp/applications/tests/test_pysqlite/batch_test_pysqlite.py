#!/usr/bin/env python
import unittest, os, sys,shutil
from os.path import join
"""
@summary: This pyUnit test module defines the unit tests for PySqlite
"""
from protlib_filesystem import getXmippPath
from xmipp import *
from pyworkflow.em.packages.xmipp3 import greenStr
from sqlite3 import dbapi2 as sqlite
from distutils.dir_util import mkpath

def equalMd(source, target):
    return MetaData(source) == MetaData(target)

def checkTrue(test, func, source):
    target = source.replace(test.WorkingDir, test.goldWorkingDir)
    #print "Comparing \n   source: %s, \n   target: %s" % (greenStr(source), greenStr(target))
    #print "Calling function: %s" % func.__name__
    test.assertTrue(func(source, target))
    
_sqlCommand1 = "CREATE TABLE runs "\
               "(run_id INTEGER PRIMARY KEY, "\
               "run_name TEXT, "\
               "run_state INT  DEFAULT 0, "\
               "script TEXT, "\
               "init DATE, "\
               "last_modified DATE, "\
               "protocol_name TEXT , "\
               "comment TEXT, "\
               "pid TEXT, "\
               "jobid INTEGER DEFAULT -1,  "\
               "CONSTRAINT unique_workingdir UNIQUE(run_name, protocol_name));"
               #PRIMARY KEY AUTOINCREMENT
_sqlCommand2 = "CREATE TABLE steps "\
               "(step_id INTEGER DEFAULT 0, "\
               "command TEXT, "\
               "parameters TEXT, "\
               "init DATE, "\
               "finish DATE, "\
               "verifyFiles TEXT, "\
               "iter INTEGER DEFAULT 1, "\
               "execution_mode INTEGER, "\
               "passDb BOOL, "\
               "run_id INTEGER REFERENCES runs(run_id) ON DELETE CASCADE, "\
               "parent_step_id INTEGER"\
               ",PRIMARY KEY(step_id, run_id)"\
               ");"

               
_sqlCommand4 = "CREATE TABLE verifyfiles "\
               "(iter INTEGER, "\
               "filename Text, "\
               "alias    Text, "\
               "run_id  INTEGER , "\
               "step_id INTEGER , "\
               "PRIMARY KEY(alias, iter, step_id, run_id)," \
               "FOREIGN KEY(step_id, run_id) REFERENCES steps(step_id, run_id) on delete cascade);" 

class TestPySqlite(unittest.TestCase):
    
    _labels = [WEEKLY]
    
    testsPath = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/applications/tests'
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        # Assume 'xmipp' and 'testXmipp' are at the same level
        os.environ['PROTOCOL_SCRIPT'] = sys.argv[0]
        #uncomment next 5 lines to store database in disk
        #curdir = os.path.dirname(getXmippPath())
        #tmpPath = join(curdir, 'testXmipp/input/test/test_pysqlite')
        #mkpath(tmpPath, 0777, True)
        #os.chdir(tmpPath)
        #tmpDataBaseName="test.sqlite"
        tmpDataBaseName=":memory:"
        if os.path.exists(tmpDataBaseName):
            os.remove(tmpDataBaseName)
        #self.cx = sqlite.connect(":memory:")
        self.cx = sqlite.connect(tmpDataBaseName)
        self.cu = self.cx.cursor()
        #self.cu.execute('pragma foreign_keys=ON')

                
    def test_001CreateTable(self):
        expected_sqls = [_sqlCommand1,_sqlCommand2]
        [self.cu.execute(s) for s in expected_sqls]
        i = self.cx.iterdump()
        actual_sqls = [s for s in i]
        expected_sqls = ['BEGIN TRANSACTION;'] + expected_sqls + \
            [ 'COMMIT;']
        [self.assertEqual(expected_sqls[i], actual_sqls[i])
            for i in xrange(len(expected_sqls))]

    def test_010Insert(self):
        expected_sqls = [_sqlCommand1,
                         'INSERT INTO "runs" VALUES(1,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL,-1);',
                         _sqlCommand2,
                         'INSERT INTO "steps" VALUES(1,NULL,NULL,NULL,NULL,NULL,1,NULL,NULL,1,NULL);',
                         _sqlCommand4
                         ]
        [self.cu.execute(s) for s in expected_sqls]
        i = self.cx.iterdump()
        actual_sqls = [s for s in i]
        expected_sqls = ['BEGIN TRANSACTION;'] + expected_sqls + \
            [ 'COMMIT;']
        [self.assertEqual(expected_sqls[i], actual_sqls[i])
            for i in xrange(len(expected_sqls))]
        
    def test_020Delete(self):
        expected_sqls = [_sqlCommand1,
                         'INSERT INTO "runs" VALUES(1,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL,-1);',
                         _sqlCommand2
                         ]
        [self.cu.execute(s) for s in expected_sqls]
        expected_sqls2 = ['INSERT INTO "steps" VALUES(1,NULL,NULL,NULL,NULL,NULL,1,NULL,NULL,1,NULL);',
                         'pragma foreign_keys=ON',
                         'DELETE FROM "steps" WHERE step_id=1'
                         ]
        [self.cu.execute(s) for s in expected_sqls2]
        i = self.cx.iterdump()
        actual_sqls = [s for s in i]
        expected_sqls = ['BEGIN TRANSACTION;'] + expected_sqls + \
            [ 'COMMIT;']
        [self.assertEqual(expected_sqls[i], actual_sqls[i])
            for i in xrange(len(expected_sqls))]

    def test_021Delete(self):
        #PRAGMA MUSt BE ACTIVE BEFORE TABLE CREATION
        self.cu.execute('pragma foreign_keys=ON')
        expected_sqls = [
                         _sqlCommand1,
                         _sqlCommand2,
                        ]
        [self.cu.execute(s) for s in expected_sqls]
        expected_sqls2 = [
                         'INSERT INTO "runs" VALUES(1,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL,-1);',
                         'INSERT INTO "steps" VALUES(1,NULL,NULL,NULL,NULL,NULL,1,NULL,NULL,1,NULL);',
                         'DELETE FROM "runs" WHERE run_id=1'
                        ]
        [self.cu.execute(s) for s in expected_sqls2]
        i = self.cx.iterdump()
        actual_sqls = [s for s in i]
        expected_sqls = ['BEGIN TRANSACTION;'] + expected_sqls + \
            [ 'COMMIT;']
        [self.assertEqual(expected_sqls[i], actual_sqls[i])
            for i in xrange(len(expected_sqls))]

    def test_022Delete(self):
        #PRAGMA MUSt BE ACTIVE BEFORE TABLE CREATION
        self.cu.execute('pragma foreign_keys=ON')
        expected_sqls = [_sqlCommand1,
                         'INSERT INTO "runs" VALUES(1,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL,-1);',
                         _sqlCommand2,
                         'INSERT INTO "steps" VALUES(2,NULL,NULL,NULL,NULL,NULL,1,NULL,NULL,1,NULL);',
                          _sqlCommand4,
                         'INSERT INTO "verifyfiles" VALUES(1,NULL,NULL,1,2);'
                        ]

        [self.cu.execute(s) for s in expected_sqls]
        expected_sqls2 = ['INSERT INTO "steps" VALUES(1,NULL,NULL,NULL,NULL,NULL,1,NULL,NULL,1,NULL);',
                         'INSERT INTO "verifyfiles" VALUES(1,NULL,NULL,1,1);',
                         'DELETE FROM "steps" WHERE step_id=1'
                        ]
        [self.cu.execute(s) for s in expected_sqls2]
        i = self.cx.iterdump()
        actual_sqls = [s for s in i]
        expected_sqls = ['BEGIN TRANSACTION;'] + expected_sqls + \
            [ 'COMMIT;']
            
#        import pprint
#        pp = pprint.PrettyPrinter(indent=4)
#        pp.pprint(actual_sqls)
#        print "======"
#        pp.pprint(expected_sqls)
        [self.assertEqual(expected_sqls[i], actual_sqls[i])
            for i in xrange(len(expected_sqls))]


from  XmippPythonTestResult import XmippPythonTestResult

                                        
if __name__ == '__main__':
    #unittest.main()   
    argc = len(sys.argv)      
    if  argc > 1:  
        xmlFile = sys.argv[1]
    else: 
        xmlFile = '/dev/null'

    suite = unittest.TestLoader().loadTestsFromTestCase(TestPySqlite)
    result = XmippPythonTestResult()
    result.openXmlReport("TestPySql", xmlFile)    
    suite(result)
    result.closeXmlReport()
    
    if result.testFailed != 0:
       result = unittest.TextTestRunner(verbosity=2).run(suite)    

