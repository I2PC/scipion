import unittest, sys
import subprocess

import pyworkflow as pw
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.utils.graph import Graph, Node
import pyworkflow.em.packages.eman2 as eman2
    
    
class TestXmippWorkflow(unittest.TestCase):

    def testXmippWorkflow(self):
        projName = "TestXmippWorkflow"
        project = Manager().loadProject(projName) # Now it will be loaded if exists
        
#        g = project.getRunsGraph()
#        g.printNodes()

        particleList = project.mapper.selectByClass('SetOfParticles')
        
        eman2.loadEnvironment()
        
        program = pw.join('em', 'packages', 'eman2', 'e2converter.py')
        
        cmd = eman2.getEmanCommand(program, 'myimage.hdf')
        
        gcmd = greenStr(cmd)
        print "** Running: '%s'" % gcmd
        proc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
        
        p = particleList[0]
        for i, img in enumerate(p):
            loc = img.getLocation()
            index = loc.getIndex()
            fn = loc.getFileName()
            print >> proc.stdin, index, fn
        p.wait()
        
            #print " reading with EMAN2"
            #imageData = EMData('myimage.hdf', index - 1, False)
            #print " writing with EMAN2"
            #imageData.write_image(eman_stack, i,EMUtil.ImageType.IMAGE_HDF,True)
                
            
#        print "=" * 100

#        f = ConditionFilter('hasAlignment and samplingRate > 4.0')
#        
#        for p in project.mapper.selectByClass('SetOfParticles', iterate=True, objectFilter=f):
#            p.printAll()
            
        
        
        
            
             
class ConditionFilter():
    def __init__(self, condition):
        self.condition = condition
        
    def __call__(self, obj):
        return obj.evalCondition(self.condition)

       
       
if __name__ == "__main__":
    unittest.main()
