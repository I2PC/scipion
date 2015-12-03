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
Bibtex string file for Gctf program.
"""

_bibtexStr = """

@article{Zhang2015,
title = "Gctf: Real-time CTF determination and correction ",
journal = "JSB ",
volume = "",
number = "",
pages = "",
year = "2015",
note = "in press",
issn = "",
doi = "http://dx.doi.org/10.1016/j.jsb.2015.11.003",
url = "http://www.sciencedirect.com/science/article/pii/S1047847715301003",
author = "Zhang, Kai",
keywords = "Contrast transfer function, Cryo-electron microscopy, GPU program, CTF determination "
}

"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
