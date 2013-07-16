import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.utils.graph import Graph, Node
    
    
class TestXmippWorkflow(unittest.TestCase):

    def testXmippWorkflow(self):
        projName = "TestXmippWorkflow"
        project = Manager().loadProject(projName) # Now it will be loaded if exists
        
#        g = project.getRunsGraph()
#        g.printNodes()

        particleList = project.mapper.selectByClass('SetOfParticles')
        
        for p in particleList:
            p.printAll()
            print "hasCTF: ", p.evalCondition('hasCTF and samplingRate > 4.0')
            
        print "=" * 100
        f = ConditionFilter('hasCTF and samplingRate > 4.0')
        for p in project.mapper.selectByClass('SetOfParticles', iterate=True, objectFilter=f):
            p.printAll()
            
             
class ConditionFilter():
    def __init__(self, condition):
        self.condition = condition
        
    def __call__(self, obj):
        return obj.evalCondition(self.condition)

       
       
if __name__ == "__main__":
    unittest.main()
