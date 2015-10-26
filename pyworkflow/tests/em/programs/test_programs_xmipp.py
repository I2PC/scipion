# ***************************************************************************
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************/


import os
from base import ProgramTest

import pyworkflow.em.packages.xmipp3 as xmipp3

RM = 'rmarabini'
COSS = 'coss'
JMRT = 'delarosatrevin'


class XmippProgramTest(ProgramTest):
    _counter = 0
    
    @classmethod
    def setUpClass(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        cls.program = cls.getProgram()
        cls.env = xmipp3.getEnviron()
        #cls._counter = 0 # count number of cases per test
        
    def runCase(self, *args, **kwargs):
        if 'mpi' not in kwargs:
            kwargs['mpi'] = 2 if self.program.startswith('xmipp_mpi') else 0
        ProgramTest.runCase(self, *args, **kwargs)
        
        
class AngularDiscreteAssign(XmippProgramTest):
    _owner = RM
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_angular_discrete_assign'
        
    def test_case1(self):
        self.runCase("-i input/aFewProjections.sel -o %o/assigned_angles.xmd --ref %o/reference.doc",
                     preruns=["xmipp_angular_project_library -i input/phantomBacteriorhodopsin.vol -o %o/reference.stk --sampling_rate 10"],
                     outputs=['assigned_angles.xmd'])
        
        
class AngularDistance(XmippProgramTest):
    _owner = COSS
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_angular_distance'
        
    def test_case1(self):
        self.runCase("--ang1 input/discreteAssignment.xmd --ang2 input/aFewProjections.sel --oroot %o/angleComparison --sym c3 -v 2",
                     outputs=['angleComparison.xmd'])        


class AngularDistributionShow(XmippProgramTest):
    _owner = RM
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_angular_distribution_show'
        
    def test_case1(self):
        self.runCase("-i input/randomAngularDistribution.sel -o %o/distribution.xmd histogram",
                     outputs=['distribution.xmd'])
        
    def test_case2(self):
        self.runCase("-i input/randomAngularDistribution.sel -o %o/distribution.bild chimera",
                     outputs=['distribution.bild'])
        
        
class AngularNeighbourhood(XmippProgramTest):
    _owner = RM
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_angular_neighbourhood'
        
    def test_case1(self):
        self.runCase("--i1 input/randomAngularDistribution.sel --i2 input/aFewProjections.sel -o %o/neighborhood.sel",
                     outputs=['neighborhood.sel'])  
        
  
class AngularProjectLibrary(XmippProgramTest):
    _owner = RM
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_angular_project_library'
        
    def test_case1(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5",
                     outputs=["output_projections.doc", "output_projections.stk"])         
        
    def test_case2(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5 "
                     "--compute_neighbors --experimental_images input/aFewProjections.sel  --angular_distance -1",
                     outputs=["output_projections.doc", "output_projections.stk", "output_projections_sampling.xmd"]) 
         
    def test_case3(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5 "
                     "--compute_neighbors --experimental_images input/aFewProjections.sel  --angular_distance 10 ",
                     outputs=["output_projections.doc", "output_projections.stk", "output_projections_sampling.xmd"])                
         
    def test_case4(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5 "
                     "--compute_neighbors --experimental_images input/aFewProjections.sel  --angular_distance 10  --near_exp_data",
                     outputs=["output_projections.doc", "output_projections.stk", "output_projections_sampling.xmd"])                
        
    def test_case5(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5 --method real_space",
                     outputs=["output_projections.doc", "output_projections.stk"]) 
         

class AngularProjectLibraryMpi(AngularProjectLibrary):
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_mpi_angular_project_library'
    
    def test_case6(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5",
                     outputs=["output_projections.doc", "output_projections.stk"])
        
    def test_case7(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5 --compute_neighbors --experimental_images input/aFewProjections.sel  --angular_distance -1",
                     outputs=["output_projections.doc", "output_projections.stk", "output_projections_sampling.xmd"])
        
    def test_case8(self):
        self.runCase("-i input/phantomBacteriorhodopsin.vol -o %o/output_projections.stk --sym c6 --sampling_rate 5 --compute_neighbors --experimental_images input/aFewProjections.sel  --angular_distance 10",
                     outputs=["output_projections.doc", "output_projections.stk", "output_projections_sampling.xmd"])
    

class AngularProjectionMatching(XmippProgramTest):
    _owner = RM
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_angular_projection_matching'
    
    def test_case1(self):
        self.runCase("-i input/aFewProjections.sel -o %o/assigned_angles.xmd --ref %o/reference.stk --thr 3",
                     preruns=["xmipp_angular_project_library -i input/phantomBacteriorhodopsin.vol --experimental_images input/aFewProjections.sel -o %o/reference.stk --sampling_rate 10 --compute_neighbors --angular_distance -1"],
                     outputs=["reference.doc", "reference_sampling.xmd", "reference.stk", "assigned_angles.xmd"])       
     
     
class AngularProjectionMatchingMpi(AngularProjectionMatching):
    @classmethod
    def getProgram(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        return 'xmipp_mpi_angular_projection_matching'     
    
    def test_case2(self):
        self.runCase("-i input/aFewProjections.sel -o %o/assigned_angles.xmd --ref %o/reference.stk --thr 3",
                     preruns=["xmipp_angular_project_library -i input/phantomBacteriorhodopsin.vol --experimental_images input/aFewProjections.sel -o %o/reference.stk --sampling_rate 10 --compute_neighbors --angular_distance -1"],
                     outputs=["reference.doc", "reference_sampling.xmd", "reference.stk", "assigned_angles.xmd"])    