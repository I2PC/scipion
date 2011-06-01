#!/usr/bin/env python
import unittest, os, sys
"""
@summary: This pyUnit test module defines the unit tests for the Xmipp Python Interface
"""
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
sys.path.append(scriptdir)
from xmipp import *
from ProjMatchActionsToBePerformedBeforeLoop import *
from ProjMatchActionsToBePerformedInLoop import *
import log, logging
from distutils.dir_util import mkpath

class TestProjMatching(unittest.TestCase):
    testsPath = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/applications/tests'
    def setUp(self):
        """This function performs all the setup stuff.      
        """
        #run this at xmipp level
        curdir = os.path.abspath(os.path.dirname(os.popen('which xmipp_protocols', 'r').read())+ '/../../')
        self.ProjectDir = os.path.join(curdir,'testXmipp/input/Protocol_Projection_Matching')
        self.WorkingDir = 'ProjMatch/new20'
        self.goldWorkingDir = 'ProjMatch/goldStandard'
        self.path = os.path.join(self.ProjectDir,self.WorkingDir)
        mkpath(self.path, 0777, True)
        os.chdir(self.ProjectDir)
        print "Changed directory to", self.ProjectDir
        self.log = log.init_log_system(self.ProjectDir,
                                '/tmp',
                                sys.argv[0],
                                self.WorkingDir)
                
    def test_00execute_ctf_groups(self):
        CtfGroupDirectory = os.path.join(self.path,'CtfGroups')
        CtfGroupRootName  = 'ctf'
        dict = {
                'CTFDatName': 'new_ctf.ctfdat'
                ,'CtfGroupDirectory': CtfGroupDirectory
                ,'CtfGroupMaxDiff': 0.10000000000000001
                ,'CtfGroupMaxResol': 5.5999999999999996
                ,'CtfGroupRootName': CtfGroupRootName
                ,'DataArePhaseFlipped': True
                ,'DoAutoCtfGroup': True
                ,'DoCtfCorrection': True
                ,'PaddingFactor': 2
                ,'SelFileName': 'new20.sel'
                ,'SplitDefocusDocFile': ''
                ,'WienerConstant': -1
           }
        execute_ctf_groups(self.log,dict)

        testName = CtfGroupDirectory+"/"+CtfGroupRootName+'Info.xmd'
        goldName = testName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldName,testName))
        #check two file compare works
        self.assertFalse(compareTwoFiles(goldName,'/etc/passwd'))

        testName1 = CtfGroupDirectory+"/"+CtfGroupRootName+'_ctf.stk'
        goldName = testName1.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldName,testName1))

        testName = CtfGroupDirectory+"/"+CtfGroupRootName+'_wien.stk'
        goldName = testName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldName,testName))
        #test Imgcompare works
        self.assertFalse(ImgCompare(testName,testName1))

        testName = CtfGroupDirectory+"/"+CtfGroupRootName+'_split.doc'
        goldName = testName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldName,testName))

        testName = CtfGroupDirectory+"/"+CtfGroupRootName+'_images.sel'
        goldName = testName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldName,testName))

    def test_10execute_mask(self):
        maskedFileNamesIter='ProjMatch/new20/Iter_01/masked_reference_ref_01.vol'
        dict = {
                'DoMask': True,
                'DoSphericalMask': True,
                'maskRadius': 64,
                'maskedFileName': maskedFileNamesIter ,
                'reconstructedFileName': 'ico.vol',
                'userSuppliedMask': 'mask.vol'
        }
        tmpDirName=os.path.dirname(maskedFileNamesIter)
        if not os.path.exists(tmpDirName):
            os.mkdir(tmpDirName)
        execute_mask(self.log,dict)
        testFileName = maskedFileNamesIter
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(ImgCompare(goldFileName,testFileName))

    def test_20angular_project_library(self):
        dict = {'AngSamplingRateDeg': '1',
                'CtfGroupSubsetFileName': 'ProjMatch/new20/CtfGroup/ctf_images.sel',
                'DoCtfCorrection': True,
                'DoParallel': True,
                'DoRestricSearchbyTiltAngle': False,
                'DocFileInputAngles': 'ProjMatch/new20/original_angles.doc',
                'MaxChangeInAngles': '1000',
                'MpiJobSize': '1',
                'NumberOfMpiProcesses': 3,
                'NumberOfThreads': 1,
                'OnlyWinner': False,
                'PerturbProjectionDirections': False,
                'ProjectLibraryRootName': 'ProjMatch/new20/Iter_01/ReferenceLibrary/gallery_ref_01.stk',
                'SymmetryGroup': 'i3',
                'SymmetryGroupNeighbourhood': '',
                'SystemFlavour': 'TORQUE-OPENMPI',
                'Tilt0': 40,
                'TiltF': 90,
                'maskedFileNamesIter': 'ProjMatch/new20/Iter_01/masked_reference_ref_01.vol'
        }
        tmpDirName ='ProjMatch/new20/Iter_01/ReferenceLibrary'
        if not os.path.exists(tmpDirName):
            os.mkdir(tmpDirName)
        angular_project_library(self.log,dict)
        tmpDirName = os.path.join(tmpDirName,'gallery_ref_01')
        
        #do not use os.path.join because adds an extra /
        testFileName = tmpDirName +'.stk'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(ImgCompare(goldFileName,testFileName))
        
        testFileName = tmpDirName +'.doc'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldFileName,testFileName))
        
        testFileName = tmpDirName +'_sampling.txt'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        print "aaa",goldFileName,testFileName,"aaa"
        self.assertTrue(compareTwoFiles(goldFileName,testFileName))
        
        testFileName = tmpDirName +'_group000001_sampling.txt'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldFileName,testFileName))
        
        testFileName = tmpDirName +'_group000002_sampling.txt'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldFileName,testFileName))
        
    def test_30projection_matching(self):
        tmpDirName ='ProjMatch/new20/Iter_01/ProjMatchClasses'
        dict = {    'AvailableMemory': 2,
                    'CtfGroupDirectory': 'ProjMatch/new20/CtfGroups',
                    'CtfGroupRootName': 'ctf',
                    'DoComputeResolution': True,
                    'DoCtfCorrection': True,
                    'DoParallel': True,
                    'DoScale': False,
                    'InnerRadius': '0',
                    'MaxChangeOffset': '1000',
                    'MpiJobSize': '1',
                    'NumberOfCtfGroups': 2L,
                    'NumberOfMpiProcesses': 3,
                    'NumberOfThreads': 1,
                    'OuterRadius': '64',
                    'PaddingFactor': 2,
                    'ProjMatchRootName': 'ProjMatch/new20/Iter_01/ProjMatchClasses/proj_match_ref_01.doc',
                    'ProjectLibraryRootName': 'ProjMatch/new20/Iter_01/ReferenceLibrary/gallery_ref_01.stk',
                    'ReferenceIsCtfCorrected': '1',
                    'ScaleNumberOfSteps': '3',
                    'ScaleStep': '1',
                    'Search5DShift': '5',
                    'Search5DStep': '2',
                    'SystemFlavour': 'TORQUE-OPENMPI'
        }

        if not os.path.exists(tmpDirName):
            os.mkdir(tmpDirName)
        projection_matching(self.log,dict)
        testFileName = dict['ProjMatchRootName']
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        print "aaaa",testFileName,goldFileName, "aaaa"
        self.assertTrue(ImgCompare(goldFileName,testFileName))

        
from  XmippPythonTestResult import XmippPythonTestResult

                                        
if __name__ == '__main__':
    #unittest.main()   
    argc = len(sys.argv)      
    if  argc > 1:  
        xmlFile = sys.argv[1]
    else: 
        xmlFile = '/dev/null'

    suite = unittest.TestLoader().loadTestsFromTestCase(TestProjMatching)
    result = XmippPythonTestResult()
    result.openXmlReport("TestProjMatching", xmlFile)    
    suite(result)
    result.closeXmlReport()
    
    if result.testFailed != 0:
       result = unittest.TextTestRunner(verbosity=2).run(suite)    

