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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for Xmipp package.
"""

_bibtexStr = """
@article{Heymann2007,
title = "Bsoft: Image processing and molecular modeling for electron microscopy ",
journal = "Journal of Structural Biology ",
volume = "157",
number = "1",
pages = "3 - 18",
year = "2007",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/j.jsb.2006.06.006",
url = "http://www.sciencedirect.com/science/article/pii/S1047847706001997",
author = "J. Bernard Heymann and David M. Belnap",
keywords = "Single particle analysis Tomography",
}
"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
