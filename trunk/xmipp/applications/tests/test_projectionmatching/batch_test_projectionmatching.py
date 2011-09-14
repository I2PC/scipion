#!/usr/bin/env xmipp_python
import unittest, os, sys,shutil
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
        self.src = 'ProjMatch/goldStandard/'
        self.dst = 'ProjMatch/new20/'

                
    def test_000execute_ctf_groups(self):
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
        self.assertTrue(ImgCompare(goldName,testName1))

        testName = CtfGroupDirectory+"/"+CtfGroupRootName+'_wien.stk'
        goldName = testName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(ImgCompare(goldName,testName))

        testName = CtfGroupDirectory+"/"+CtfGroupRootName+'_split.doc'
        goldName = testName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldName,testName))

        testName = CtfGroupDirectory+"/"+CtfGroupRootName+'_images.sel'
        goldName = testName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldName,testName))

    def test_010execute_mask(self):
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

    def test_020angular_project_library(self):
        dict = {'AngSamplingRateDeg': '1',
                'CtfGroupSubsetFileName': 'ProjMatch/new20/CtfGroups/ctf_images.sel',
                'DoCtfCorrection': True,
                'DoParallel': True,
                'DoRestricSearchbyTiltAngle': False,
                'DocFileInputAngles': 'ProjMatch/new20/original_angles.doc',
                'MaxChangeInAngles': '1000',
                'MpiJobSize': '1',
                'NumberOfMpi': 3,
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
        src = self.src + 'original_angles.doc'
        dst = self.dst 

        shutil.copy(src, dst)
        tmpDirName = self.dst + 'Iter_01/ReferenceLibrary'
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
        auxGold=MetaData(goldFileName)
        auxTest=MetaData(testFileName)
        self.assertTrue(auxGold==auxTest)
        
        testFileName = tmpDirName +'_sampling.xmd'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        auxGold=MetaData("extra@"+goldFileName)
        auxTest=MetaData("extra@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        auxGold=MetaData("neighbors@"+goldFileName)
        auxTest=MetaData("neighbors@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        auxGold=MetaData("projectionDirections@"+goldFileName)
        auxTest=MetaData("projectionDirections@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        
        testFileName = tmpDirName +'_group000001_sampling.xmd'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        auxGold=MetaData("extra@"+goldFileName)
        auxTest=MetaData("extra@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        auxGold=MetaData("neighbors@"+goldFileName)
        auxTest=MetaData("neighbors@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        auxGold=MetaData("projectionDirections@"+goldFileName)
        auxTest=MetaData("projectionDirections@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        
        testFileName = tmpDirName +'_group000002_sampling.xmd'
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        auxGold=MetaData("extra@"+goldFileName)
        auxTest=MetaData("extra@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        auxGold=MetaData("neighbors@"+goldFileName)
        auxTest=MetaData("neighbors@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        auxGold=MetaData("projectionDirections@"+goldFileName)
        auxTest=MetaData("projectionDirections@"+testFileName)
        self.assertTrue(auxGold==auxTest)
        
    def test_030projection_matching(self):
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
                    'NumberOfMpi': 3,
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
        testFileName = dict['ProjMatchRootName']
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        md1= MetaData(testFileName)
        md2=MetaData(goldFileName)
        projection_matching(self.log,dict)
        self.assertTrue(md1==md2)

    def test_040assign_images_to_references(self):
        dict = {   
                'DocFileInputAngles': 'ProjMatch/new20/Iter_01/current_angles.doc',
                'NumberOfCtfGroups': 2L,
                'NumberOfReferences': 3,
                'ProjMatchRootName': [   None,
                    'ProjMatch/new20/Iter_01/ProjMatchClasses/proj_match_ref_01.doc',
                    'ProjMatch/new20/Iter_01/ProjMatchClasses/proj_match_ref_02.doc',
                    'ProjMatch/new20/Iter_01/ProjMatchClasses/proj_match_ref_03.doc']}
        #cp from goldstandard
        src = self.src + 'Iter_01/ProjMatchClasses/'
        dst = self.dst + 'Iter_01/ProjMatchClasses'
        shutil.copy(src+'proj_match_ref_01.doc', dst)
        shutil.copy(src+'proj_match_ref_02.doc', dst)
        shutil.copy(src+'proj_match_ref_03.doc', dst)
        assign_images_to_references(self.log,dict)
        testFileName = dict['DocFileInputAngles']
        goldFileName = testFileName.replace(self.WorkingDir,self.goldWorkingDir)
        self.assertTrue(compareTwoFiles(goldFileName,testFileName))
        
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

