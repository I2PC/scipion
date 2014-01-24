import unittest, sys
import subprocess

import pyworkflow as pw
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.utils import cleanPath
from pyworkflow.utils.graph import Graph, Node
import pyworkflow.em.packages.eman2 as eman2
from pyworkflow.mapper.sqlite import SqliteFlatMapper
    
    
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
            
        print s.getObjDict()
        
    def testFlatDb(self):
        projName = "relion_ribo"
        self.proj = Manager().loadProject(projName) # Now it will be loaded if exists
        
#        protImport = ProtImportParticles(pattern=getInputPath('Images_Vol_ML3D/phantom_images', '*.xmp'), checkStack=True, samplingRate=1.237)
#        self.proj.launchProtocol(protImport, wait=True)
#        
        setOfPart = self.proj.mapper.selectByClass('SetOfParticles')[0]
        
        setOfPart2 = SetOfParticles()
        setOfPart2._mapper = SqliteFlatMapper("partFlat.sqlite", globals())
        for p in setOfPart:
            print "takarras"
            setOfPart2.append(p.clone())
        setOfPart2.write()
        
        
    def creatingFlatDb(self):
        n = 1000
        imgSet = SetOfParticles(filename='particles_flat.sqlite', mapperClass=SqliteFlatMapper)
        
        for i in range(n):
            img = self._createImage()
            img.setIndex(i+1)
            imgSet.append(img)
            a = 1
            
        imgSet.write()

    def nestedFlatDb(self): 
        fn = 'classes_flat.sqlite'
        cleanPath(fn)
        
        images = SetOfImages()
        images.setSamplingRate(1.2)
        classes2DSet = SetOfClasses2D(filename=fn)
        classes2DSet.setImages(images)
        averages = classes2DSet.createAverages()
    
        for ref in range(1, 11):
            print "class: ", ref
            class2D = Class2D()
            class2D.setObjId(ref)
            print "append class to set, ref=", ref
            classes2DSet.append(class2D)
            avg = Particle()
            avg.setLocation(ref, 'averages.stk')
            print "   populating class "
            for i in range(1, 101):         
                img = Particle()
                img.setSamplingRate(5.3)
                class2D.append(img)
                    
            print "   writing class "
            class2D.write()
            print "   append avg"
            averages.append(avg)
        
        classes2DSet.write()
        
    def selectingFlatDb(self):
        imgSet = SetOfParticles(mapperClass=SqliteFlatMapper)
        imgSet._filename.set('particles_flat.sqlite')
        
        for img in imgSet:
            #img.printAll()
            #img = self._createImage()
            a = img.getSamplingRate()

    def _createImage(self):
        img = Image()
        img.setLocation(1, 'image.spi')
        img.setSamplingRate(3.5)
        ctf = CTFModel()
        ctf.defocusU.set(1000)
        ctf.defocusAngle.set(90)
        img.setCTF(ctf)
        img.setAttributeValue('_ctfModel.defocusV', 1000)
        return img
              
    def testObjDict(self):
        img = self._createImage()
        img.printObjDict()
        
        objDict = img.getObjDict()
        objDict['_ctfModel.defocusAngle'] = '10'
        
        def printTrace():
            print "sampling changed to: ", img.getSamplingRate()
            
        img._samplingRate.trace(printTrace)
        
        objDict['_samplingRate'] = '1.0'
        for k, v in objDict.iteritems():
            img.setAttributeValue(k, v)
        img.printAll()
        
        img.setSamplingRate(5.0)
        
    
   
             
class ConditionFilter():
    def __init__(self, condition):
        self.condition = condition
        
    def __call__(self, obj):
        return obj.evalCondition(self.condition)

       
       
if __name__ == "__main__":
    unittest.main()
