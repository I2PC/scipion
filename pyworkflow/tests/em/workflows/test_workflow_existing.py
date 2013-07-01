import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.utils.graph import Graph, Node
    
    
class TestXmippWorkflow(unittest.TestCase):
    
    def testXmippWorkflow(self):
        projName = "TestXmippWorkflow"
        project = Manager().loadProject(projName) # Now it will be loaded if exists
        
        g = project.getRunsGraph()
        
        g.printNodes()

       
       
if __name__ == "__main__":
    unittest.main()
