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
    
    _labels = [WEEKLY]
    
    def loadProject(self, projName):
        """ Try to load an existing project. """
        manager = Manager()
        return manager.loadProject(projName)

    def a_testEmanConvert(self):
        project = self.loadProject("TestXmippWorkflow")

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
        project = self.loadProject("TestXmippWorkflow") # Now it will be loaded if exists
        
        coords = project.mapper.selectByClass('SetOfCoordinates')[0]

        for c in coords.iterCoordinates():
            print "coord: ", c.getPosition(), " from mic: ", c.getMicrograph().getFileName()        
            
        
    def test_sets(self):
        project = self.loadProject("TestXmippWorkflow") # Now it will be loaded if exists
        
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
            
    def test_xmippWriteImages(self):
        projName = "relion_ribo"
        project = self.loadProject(projName) # Now it will be loaded if exists
        
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
        
    def test_flatDb(self):
        projName = "relion_ribo"
        self.proj = Manager().loadProject(projName) # Now it will be loaded if exists
        
#        protImport = ProtImportParticles(filesPath=getInputPath('Images_Vol_ML3D/phantom_images', '*.xmp'), checkStack=True, samplingRate=1.237)
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
        ctf.setDefocusU(1000)
        ctf.setDefocusAngle(90)
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
        
        
    def test_relations(self):
        projName = "TestXmippWorkflow"
        project = self.loadProject(projName)
        
        project.getTransformGraph()
    
    def test_exportEMX(self):
        projName = "TestXmippWorkflow"
        project = self.loadProject(projName)
        
        protExtract = project.mapper.selectByClass('XmippProtExtractParticles')[0]
        
        #protExtract2 = project.copyProtocol(protExtract)
        #project.launchProtocol(protExtract2, wait=True)
        
        protEmx = ProtEmxExport()
        protEmx.inputSet.set(protExtract.outputParticles)
        
        project.launchProtocol(protEmx, wait=True)
        
    def test_classesWard(self):
        projName = "high throughput"
        project = self.loadProject(projName)
        
        protName = 'XmippProtCL2D'
        protName = 'SpiderProtClassifyWard'
        prot = project.mapper.selectByClass(protName)[-1]
        print "prot.getObjId(): ", prot.getObjId()
        
        classes = prot.outputClasses
        averages = prot.outputClasses.getRepresentatives()
        import pyworkflow.em.packages.xmipp3 as xmipp3
        #xmipp3.writeSetOfClasses2D(classes, "classes.xmd")
        
        print "=" * 100
        
        i = 0
        for cls, avg in zip(classes, averages):
            #for img in cls:
            #    print img.getLocation()
            #cls.printAll()
            print "class id: ", cls.getObjId(), "avg: ", avg.getLocation()
            if i == 9:
                break
            i += 1
        
        i = 0           
        for cls in classes:
            print "class id: ", cls.getObjId()
            if i == 9:
                break
            i += 1    
            
        i = 0
        for avg in averages:
            print "avg: ", avg.getLocation()
            if i == 9:
                break
            i += 1
            
            
    def test_mpiStepsExecution(self):
        project = self.loadProject("TestXmippWorkflow")
        protCTF1 = project.mapper.selectByClass('XmippProtCTFMicrographs')[0]   
        protCTF2 = project.copyProtocol(protCTF1)
        protCTF2.setObjLabel('ctf - Day2')
        protCTF2.numberOfMpi.set(3)
        protCTF2.inputMicrographs.set(protCTF1.inputMicrographs.get())
        project.launchProtocol(protCTF2, wait=True)                 
    
    def test_cleanDay2(self):
        """ Delete all runs from Day2. """
        project = self.loadProject("HighThroughputTest")
        protocols = project.getRuns()
        for prot in reversed(protocols):
            if 'Day2' in prot.getObjLabel():
                print "Deleting protocol: ", prot.getObjLabel()
                project.deleteProtocol(prot)
        
        
    def test_autoDay2(self):
        project = self.loadProject("HighThroughputTest")
        
        def getProtocol(className):
            """ Return the first protocol found from a give className. """
            return project.mapper.selectByClass(className)[0]

        protImport1 = getProtocol('ProtImportMovies')
        pattern = protImport1.pattern.get()
        protImport2 = project.copyProtocol(protImport1)
        protImport2.setObjLabel('import movies - Day2')
        protImport2.pattern.set(pattern.replace('1', '2'))
        project.launchProtocol(protImport2, wait=True)
        
        # Copy of align movies
        protAlignMov1 = getProtocol('ProtOpticalAlignment')
        protAlignMov2 = project.copyProtocol(protAlignMov1)
        protAlignMov2.setObjLabel('Align Movies - Day2')
        protAlignMov2.inputMovies.set(protImport2.outputMovies)
        project.launchProtocol(protAlignMov2, wait=True)
        
        # Copy of preprocess
        protPrep1 = getProtocol('XmippProtPreprocessMicrographs')
        protPrep2 = project.copyProtocol(protPrep1)
        protPrep2.setObjLabel('crop mics 50 px - Day2')
        protPrep2.inputMicrographs.set(protAlignMov2.outputMicrographs)
        project.launchProtocol(protPrep2, wait=True)
        
        # Copy of CTF estimation.
        protCTF1 = getProtocol('XmippProtCTFMicrographs')
        protCTF2 = project.copyProtocol(protCTF1)
        protCTF2.setObjLabel('ctf - Day2')
        protCTF2.numberOfThreads.set(3)
        protCTF2.inputMicrographs.set(protPrep2.outputMicrographs)
        project.launchProtocol(protCTF2, wait=True)
        
        # Launch an automatic picking using the previous train
        protPick1 = getProtocol('XmippProtParticlePicking')
        protPick2 = XmippParticlePickingAutomatic()
        protPick2.setObjLabel('autopicking - Day2')
        protPick2.xmippParticlePicking.set(protPick1)
        protPick2.micsToPick.set(1) # OTHER
        protPick2.inputMicrographs.set(protPrep2.outputMicrographs)
        project.launchProtocol(protPick2, wait=True)
        
        # Extract particles of day 2
        protExtract1 = getProtocol('XmippProtExtractParticles')
        protExtract2 = project.copyProtocol(protExtract1)
        protExtract2.setObjLabel('extract particles - Day2')
        protExtract2.numberOfThreads.set(3)
        protExtract2.inputCoordinates.set(protPick2.outputCoordinates)
        protExtract2.ctfRelations.set(protCTF2.outputCTF)
        protExtract2.inputMicrographs.set(protPrep2.outputMicrographs)
        project.launchProtocol(protExtract2, wait=True)
        
        # Run Spider-filter for day 2
        protFilter1 = getProtocol('SpiderProtFilter')
        protFilter2 = project.copyProtocol(protFilter1)
        protFilter2.setObjLabel('spi filter - Day2')
        protFilter2.inputParticles.set(protExtract2.outputParticles)
        project.launchProtocol(protFilter2, wait=True)
        
        # align for day 2
        protAlign1 = getProtocol('XmippProtCL2DAlign')
        protAlign2 = project.copyProtocol(protAlign1)
        protAlign2.setObjLabel('cl2d align - Day2')
        protAlign2.numberOfMpi.set(8)
        protAlign2.inputParticles.set(protFilter2.outputParticles)
        project.launchProtocol(protAlign2, wait=True)
        
        # capca for day 2
        protCAPCA1 = getProtocol('SpiderProtCAPCA')
        protCAPCA2 = project.copyProtocol(protCAPCA1)
        protCAPCA2.setObjLabel('spi capca - Day2')
        protCAPCA2.inputParticles.set(protAlign2.outputParticles)
        project.launchProtocol(protCAPCA2, wait=True)
        
        # ward for day 2
        protWard1 = getProtocol('SpiderProtClassifyWard')
        protWard2 = project.copyProtocol(protWard1)
        protWard2.setObjLabel('spi ward - Day2')
        protWard2.inputParticles.set(protAlign2.outputParticles)
        protWard2.pcaFilePointer.set(protCAPCA2.imcFile)
        project.launchProtocol(protWard2, wait=True)        
             
    def test_browseRelations(self):
        """ Check that relation are correctly retrieved. """
        project = self.loadProject("TestXmippExtractParticles")
        graph = project.getTransformGraph()
        obj = project.mapper.selectByClass('SetOfMicrographs')[1]
        
        #print "obj: ", obj
        #print "graph:"
        #graph.printDot()
        
        objs = project._getConnectedObjects(obj, graph)
        for o in objs.values():
            print "o.label = ", o
            
        related = project.getRelatedObjects(RELATION_CTF, obj)
        for r in related:
            print "r: ", r
             
    def test_importEMX(self):
        manager = Manager()
        project = manager.createProject("emx import")
        
        prot = ProtEmxImport(objLabel='emx import mics')
        prot.inputEMX.set('/home/josem/emxDataMicrographs/data.emx')
        project.launchProtocol(prot, wait=True)

        prot = ProtEmxImport(objLabel='emx import parts')
        prot.inputEMX.set('/home/josem/emxDataParticles/data.emx')
        project.launchProtocol(prot, wait=True)

        prot = ProtEmxImport(objLabel='emx import coords')
        prot.inputEMX.set('/home/josem/emxDataCoordinates/data.emx')
        project.launchProtocol(prot, wait=True)        
        
                    
class ConditionFilter():
    def __init__(self, condition):
        self.condition = condition
        
    def __call__(self, obj):
        return obj.evalCondition(self.condition)

       
       
if __name__ == "__main__":
    unittest.main()
