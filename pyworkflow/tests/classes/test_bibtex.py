#!/usr/bin/env python
# coding: latin-1

import unittest
from os.path import join, dirname, exists
from StringIO import StringIO
from bibtexparser.bparser import BibTexParser   
    
class TestPyworkflow(unittest.TestCase):
    """ Some minor tests to the bibtexparser library. """

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
        f = StringIO()
        f.write(bibtex)
        f.seek(0, 0)
        parser = BibTexParser(f)
        from pyworkflow.utils import prettyDict
        prettyDict(parser.get_entry_dict())
        
        
        
if __name__ == '__main__':
    unittest.main()
    
