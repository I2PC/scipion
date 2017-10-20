# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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

#   -+-+-+- C H A N G E   T H I S   S T R I N G -+-+-+-
_bibtexStr = """
@article{Tang2007,
title = "EMAN2: An extensible image processing suite for electron microscopy ",
journal = "JSB",
volume = "157",
number = "1",
pages = "38 - 46",
year = "2007",
note = "Software tools for macromolecular microscopy ",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/j.jsb.2006.05.009",
url = "http://www.sciencedirect.com/science/article/pii/S1047847706001894",
author = "Guang Tang and Liwei Peng and Philip R. Baldwin and Deepinder S. Mann and Wen Jiang and Ian Rees and Steven J. Ludtke",
keywords = "EMAN, Single particle analysis , cryoEMTEM, Software, Image processing, Electron microscopy}


"""




from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
