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

@article{Mindell2003,
title = "Accurate determination of local defocus and specimen tilt in electron microscopy ",
journal = "Journal of Structural Biology ",
volume = "142",
number = "3",
pages = "334 - 347",
year = "2003",
note = "",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/S1047-8477(03)00069-8",
url = "http://www.sciencedirect.com/science/article/pii/S1047847703000698",
author = "Joseph A. Mindell and Nikolaus Grigorieff",
keywords = "Electron microscopy, Contrast transfer function, Algorithm, Tilt determination "
}


"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
