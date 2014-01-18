import unittest, sys
import subprocess

import pyworkflow as pw
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.utils.graph import Graph, Node
import pyworkflow.em.packages.eman2 as eman2
    
    
class TestXmippWorkflow(unittest.TestCase):

    def a_testEmanConvert(self):
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

    def test_findRowById(self):
        projName = "TestXmippWorkflow"
        project = Manager().loadProject(projName) # Now it will be loaded if exists
        
        coords = project.mapper.selectByClass('SetOfCoordinates')[0]

        for c in coords.iterCoordinates():
            print "coord: ", c.getPosition(), " from mic: ", c.getMicrograph().getFileName()        
            
        
    def test_Sets(self):
        projName = "TestXmippWorkflow"
        project = Manager().loadProject(projName) # Now it will be loaded if exists
        
        sets = project.mapper.selectByClass('Set')
        
        for s in sets:
            print s
            
        mics = sets[2]
        
        for m in mics:
            print "mic: ", m
            #m.printAll()
            
        m1 = mics[1]
        m1.printAll()
            
        imgs = sets[5]
        
        for i, img in enumerate(imgs):
            print "img: ", img
            if i == 10:
                break # Limit to 10 images
            
    def testXmippWriteImages(self):
        projName = "relion_ribo"
        project = Manager().loadProject(projName) # Now it will be loaded if exists
        
        sets = project.mapper.selectByClass('SetOfParticles')
        
        print sets[0]
#        print "writing set to .xmd"
#        from pyworkflow.em.packages.xmipp3 import writeSetOfParticles
#        
#        writeSetOfParticles(sets[0], "images.xmd")
        print "iterating set:"
        s = sets[0]
        dbFn = s.getFileName()
        print "FileName: ", dbFn
        from sqlite3 import dbapi2 as sqlite
        conn = sqlite.Connection(dbFn, check_same_thread = False)
        conn.row_factory = sqlite.Row
        cur = conn.cursor()
        cur.execute('select * from Objects;')
        #conn2 = sqlite.Connection(':memory:')
        #cur2 = conn2.cursor()
        i = 0
#        for line in conn.iterdump():
        row = cur.fetchone()
        while row is not None:
            #yield row
            row = cur.fetchone()
            #cur2.executescript(line)
            i += 1
                
        print "Total lines: ", i
            
    def testObjDict(self):
        img = Image()
        img.setLocation(1, 'image.spi')
        img.setSamplingRate(3.5)
        ctf = CTFModel()
        ctf.defocusU.set(1000)
        ctf.defocusAngle.set(90)
        img.setCTF(ctf)
        img.setAttributeValue('_ctfModel.defocusV', 1000)
        
        objDict = img.getObjDict()
        
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(objDict)
        
        objDict['_ctfModel.defocusAngle'] = '10'
        objDict['_samplingRate'] = '1.0'
        for k, v in objDict.iteritems():
            img.setAttributeValue(k, v)
        img.printAll()
        
   
             
class ConditionFilter():
    def __init__(self, condition):
        self.condition = condition
        
    def __call__(self, obj):
        return obj.evalCondition(self.condition)

       
       
if __name__ == "__main__":
    unittest.main()
