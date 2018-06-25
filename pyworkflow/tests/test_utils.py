#!/usr/bin/env python
# coding: latin-1
'''
Created on Mar 25, 2014

@author: airen
'''
import pyworkflow as pw

from subprocess import Popen
import pyworkflow.utils as pwutils
from pyworkflow.utils.process import killWithChilds
from pyworkflow.tests import *
from pyworkflow.utils import utils, prettyDict


    
class TestBibtex(BaseTest):
    """ Some minor tests to the bibtexparser library. """

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        
        
    def test_Parsing(self):
        bibtex = """

@article{delaRosaTrevin2013,
title = "Xmipp 3.0: An improved software suite for image processing in electron microscopy ",
journal = "Journal of Structural Biology ",
volume = "184",
number = "2",
pages = "321 - 328",
year = "2013",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/j.jsb.2013.09.015",
url = "http://www.sciencedirect.com/science/article/pii/S1047847713002566",
author = "J.M. de la Rosa-Trevín and J. Otón and R. Marabini and A. Zaldívar and J. Vargas and J.M. Carazo and C.O.S. Sorzano",
keywords = "Electron microscopy, Single particles analysis, Image processing, Software package "
}

@incollection{Sorzano2013,
title = "Semiautomatic, High-Throughput, High-Resolution Protocol for Three-Dimensional Reconstruction of Single Particles in Electron Microscopy",
booktitle = "Nanoimaging",
year = "2013",
isbn = "978-1-62703-136-3",
volume = "950",
series = "Methods in Molecular Biology",
editor = "Sousa, Alioscka A. and Kruhlak, Michael J.",
doi = "10.1007/978-1-62703-137-0_11",
url = "http://dx.doi.org/10.1007/978-1-62703-137-0_11",
publisher = "Humana Press",
keywords = "Single particle analysis; Electron microscopy; Image processing; 3D reconstruction; Workflows",
author = "Sorzano, CarlosOscar and Rosa Trevín, J.M. and Otón, J. and Vega, J.J. and Cuenca, J. and Zaldívar-Peraza, A. and Gómez-Blanco, J. and Vargas, J. and Quintana, A. and Marabini, Roberto and Carazo, JoséMaría",
pages = "171-193",
}
"""

        prettyDict(utils.parseBibTex(bibtex))
        

def wait (condition, timeout=30):
    """ Wait until "condition" returns False or return after timeout (seconds) param"""
    t0 = time.time()

    while condition():
        time.sleep(1)

        # Check timeout
        tDelta = time.time()

        if tDelta - t0 >= timeout:
            print("Wait timed out after ", timeout, " seconds")
            return


class TestProccess(BaseTest):
    """ Some tests for utils.process module. """

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        
    def test_Process(self):
        prog = pw.join('apps', 'pw_sleep.py')
        p = Popen('xmipp_python %s 500' % prog, shell=True)
        print "pid: ", p.pid
        time.sleep(5)
        killWithChilds(p.pid)


class TestGetListFromRangeString(BaseTest):

    def test_getListFromRangeString(self):
        inputStrings = ["1,5-8,10",         "2,6,9-11",        "2 5, 6-8"]
        outputLists = [[1, 5, 6, 7, 8, 10], [2, 6, 9, 10, 11], [2, 5, 6, 7, 8]]

        for s, o in zip(inputStrings, outputLists):
            self.assertEqual(o, pwutils.getListFromRangeString(s))
            # Check that also works properly with spaces as delimiters
            s2 = s.replace(',', ' ')
            self.assertEqual(o, pwutils.getListFromRangeString(s2))


if __name__ == '__main__':
    unittest.main()        
