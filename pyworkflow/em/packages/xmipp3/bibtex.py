# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for Xmipp package.
"""

_bibtexStr = """

@article{Otsu1979,
title = "A Threshold Selection Method from Gray-Level Histograms",
journal = "Systems, Man and Cybernetics, IEEE Transactions",
volume = "9",
number = "2",
pages = "62 - 66",
year = "1979",
issn = "0018-9472",
doi = "http://dx.doi.org/10.1109/TSMC.1979.4310076",
author = "de la Rosa-Trevín, J.M.  and Otón, J. and R. Marabini and A. Zaldívar and J. Vargas and J.M. Carazo and Sorzano, C.O.S.",
}

@article{Abrishami2013,
author = {Abrishami, V. and Zald�var-Peraza, A. and de la Rosa-Trev�n, J. M. and Vargas, J. and Ot�n, J. and Marabini, R. and Shkolnisky, Y. and Carazo, J. M. and Sorzano, C. O. S.}, 
title = {A pattern matching approach to the automatic selection of particles from low-contrast electron micrographs},
volume = {29}, 
number = {19}, 
pages = {2460-2468}, 
year = {2013}, 
doi = {http://dx.doi.org/10.1093/bioinformatics/btt429}, 
url = {http://bioinformatics.oxfordjournals.org/content/29/19/2460.abstract}, 
journal = {Bioinformatics} 
}
}

@incollection{Sorzano2013,
title = "Semiautomatic, High-Throughput, High-Resolution Protocol for Three-Dimensional Reconstruction of Single Particles in Electron Microscopy",
booktitle = "Nanoimaging",
year = "2013",
isbn = "978-1-62703-136-3",
volume = "950",
journal = "Methods in Molecular Biology",
editor = "Sousa, Alioscka A. and Kruhlak, Michael J.",
doi = "http://dx.doi.org/10.1007/978-1-62703-137-0_11",
publisher = "Humana Press",
keywords = "Single particle analysis; Electron microscopy; Image processing; 3D reconstruction; Workflows",
author = "Sorzano, C.O.S. and de la Rosa-Trevín, J.M. and Otón, J. and Vega, J.J. and Cuenca, J. and Zaldívar-Peraza, A. and Gómez-Blanco, J. and Vargas, J. and Quintana, A. and Marabini, Roberto and Carazo, JoséMaría",
pages = "171-193",
}

@ARTICLE{Vargas2013a,
author = {Vargas, J. and Ot�n, J. and Marabini, R. and Jonic, S. and {de la
  Rosa-Trev�n}, J. M. and et.al.},
title = {{FASTDEF}: Fast defocus and astigmatism estimation for high-throughput
  transmission electron microscopy.},
journal = "JSB",
doi = "http://dx.doi.org/10.1016/j.jsb.2012.12.006",
year = {2013},
volume = {181},
pages = {136--148},
number = {2},
month = {Feb},
}

@article{Vargas2013b,
title = "Particle quality assessment and sorting for automatic and semiautomatic particle-picking techniques ",
journal = "JSB",
volume = "183",
number = "3",
pages = "342 - 353",
year = "2013",
note = "",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/j.jsb.2013.07.015",
url = "http://www.sciencedirect.com/science/article/pii/S1047847713001950",
author = "J. Vargas and V. Abrishami and R. Marabini and J.M. de la Rosa-Trev�n and A. Zaldivar and J.M. Carazo and C.O.S. Sorzano",
keywords = "Electron microscopy, Particle picking, Machine learning, Single particle analysis "
}

@article{delaRosaTrevin2013,
title = "Xmipp 3.0: An improved software suite for image processing in electron microscopy ",
journal = "JSB",
volume = "184",
number = "2",
pages = "321 - 328",
year = "2013",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/j.jsb.2013.09.015",
url = "http://www.sciencedirect.com/science/article/pii/S1047847713002566",
author = "de la Rosa-Trevín, J.M.  and Otón, J. and R. Marabini and A. Zaldívar and J. Vargas and J.M. Carazo and Sorzano, C.O.S.",
keywords = "Electron microscopy, Single particles analysis, Image processing, Software package "
}

@article{chen2013fast,
  title={Fast and accurate reference-free alignment of subtomograms},
  author={Chen, Yuxiang and Pfeffer, Stefan and Hrabe, Thomas and Schuller, Jan Michael and F{\"o}rster, Friedrich},
  journal={Journal of structural biology},
  year={2013},
  publisher={Elsevier}
}
"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
