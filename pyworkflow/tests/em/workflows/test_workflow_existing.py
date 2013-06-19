import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.xmipp3 import *
from pyworkflow.utils.graph import Graph, Node
    
    
class TestXmippWorkflow(unittest.TestCase):
    
    def testXmippWorkflow(self):
        projName = "TestXmippWorkflow"
        project = Manager().loadProject(projName) # Now it will be loaded if exists
        
        outputDict = {} # Store the output dict
        
        runs = project.getRuns()
        g = Graph()
        
        for r in runs:
            key = r.getName()
            n = g.addNode(key)
            n.run = r
            for key, attr in r.iterOutputAttributes(EMObject):
                outputDict[attr.getName()] = n # mark this output as produced by r
                #print "   %s: %s" % (key, attr.getName())
            
        for r in runs:
            node = g.getNode(r.getName())
            #print '\n=========================\n', r.getName()
            #print "> Inputs:"
            for key, attr in r.iterInputAttributes():
                attrName = attr.get().getName()
                if attrName in outputDict:
                    parentNode = outputDict[attrName]
                    parentNode.addChilds(node)
                    
        for r in runs:
            node = g.getNode(r.getName())
            print "Node: ", node
            print " Childs: ", ','.join([str(c) for c in node.getChilds()])
                    
#            print "   %s: %s" % (key, attr.get().getName())
#            print "> Outputs:"

       
       
if __name__ == "__main__":
    unittest.main()
