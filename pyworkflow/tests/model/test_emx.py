#!/usr/bin/env python
# To run only the tests in this file, use:
# python -m unittest test_mappers -v
# To run a single test,
# python -m unittest -v test_mappers.TestMappers.test_connectUsing

import os
import unittest
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
    
    
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for _elem in elem:
            indent(_elem, level+1)
        if not _elem.tail or not _elem.tail.strip():
            _elem.tail = i
               
                
class TestEMX(unittest.TestCase):

    def test_writeXML(self):
        fn = 'test.xml'
        if os.path.exists(fn):
            os.remove(fn)
        f = open(fn, 'w+')  
        html = ET.Element('html')
        body = ET.SubElement(html, 'body')
        body.text = 'abc'
        indent(html)
        
        ET.ElementTree(html).write(f)
        ET.ElementTree(body).write(f)
        f.close()
        
        
    def test_writeMicrographs(self):
        from pyworkflow.em.data import Acquisition, Micrograph, SetOfMicrographs, CTFModel
        import pyworkflow.em.emx as emx
        from pyworkflow.utils import cleanPath
        
        fnXml = 'test_micrographs.xml'
        fnSql = 'test_micrographs.sqlite'
        cleanPath(fnXml, fnSql)
        
        micSet = SetOfMicrographs(filename=fnSql)
        acquisition = Acquisition(magnification=60000, voltage=300,
                          sphericalAberration=2., amplitudeContrast=0.07)
        micSet.setAcquisition(acquisition)
        micSet.setSamplingRate(1.2)
        ctf = CTFModel()
        ctf.setDefocusAngle(15)
        ctf.setDefocusU(1000)
        ctf.setDefocusV(2000)
        
        for i in range(10):
            mic = Micrograph()
            #mic.setLocation(1, "mic%03d.mrc" % (i+1))
            mic.setLocation(i+1, "mics.stk")
            mic.setCTF(ctf)
            micSet.append(mic)
            
        emx.writeSetOfMicrographs(micSet, fnXml)
        
        
    def test_readEMX(self):
        import pyworkflow.em.emx as emx
        fn = 'emxData/data.emx'
        
        def handleMicrograph(elem):
            print 'elem.tag: ', elem.tag
            print 'elem.attrib', elem.attrib
         
        def handleParticle(elem):
            print 'elem.tag: ', elem.tag     
            print 'elem.attrib', elem.attrib
              
        emx._iterXml(fn, {'micrograph': handleMicrograph, 
                          'particle': handleParticle})
