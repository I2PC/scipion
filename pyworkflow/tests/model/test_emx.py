#!/usr/bin/env python
# To run only the tests in this file, use:
# python -m unittest test_mappers -v
# To run a single test,
# python -m unittest -v test_mappers.TestMappers.test_connectUsing

import os
import unittest

from pyworkflow.tests import *
from pyworkflow.em.data import Acquisition, Micrograph, SetOfMicrographs, CTFModel
import pyworkflow.em.packages.emxlib as emx

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
               
#FIXME:Jose Miguel                
class TestEMX(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')  

    def test_writeXML(self):
        fn = self.getOutputPath('test.xml')
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
        
        
        emxDir = self.getOutputPath('emx')
        
        fnXml = self.getOutputPath('test_micrographs.xml')
        fnSql = self.getOutputPath('test_micrographs.sqlite')
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
        
        for i in range(3):
            mic = Micrograph()
            micFn = self.dataset.getFile("mic%d" % (i+1))
            print micFn
            mic.setLocation(micFn)
            #mic.setLocation(i+1, "mics.stk")
            mic.setCTF(ctf)
            micSet.append(mic)
            
        emx.exportData(emxDir, micSet)
        
        
    def test_Transform(self):
        from pyworkflow.em.packages.xmipp3.convert import xmippGeoFromMatrix
        
        m = Manager()
        p = m.loadProject('emx')
        alignment = p.mapper.selectByClass('SetOfAlignment')[-1]
        
        for t in alignment:
            print "="*100
            print "matrix: \n", t.getMatrix()
            shifts, angles = xmippGeoFromMatrix(t.getMatrix()._matrix)
            print "shifts: ", shifts
            print "angles: ", angles
        
        
        
        
    def test_readEMX(self):
        pass
