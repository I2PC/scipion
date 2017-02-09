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

_bibtexStr = """

@article{Kivioja2000,
title = "Local Average Intensity-Based Method for Identifying Spherical Particles in Electron Micrographs",
journal = "Journal of Structural Biology",
volume = "131",
number = "2",
pages = "126 - 134",
year = "2000",
note = "",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1006/jsbi.2000.4279",
url = "http://www.sciencedirect.com/science/article/pii/S1047847700942795",
author = "Teemu Kivioja and Janne Ravantti and Anatoly Verkhovsky and Esko Ukkonen and Dennis Bamford",
keywords = "spherical viruses",
keywords = "cryoelectron microscopy",
keywords = "automated object detection",
keywords = "image processing",
keywords = "pattern matching software",
abstract = "Amethod is presented that reliably detects spherical viruses from a wide variety of noisy low-contrast electron micrographs. Such detection is one of the first image analysis steps in the computer-aided reconstruction of three-dimensional density distribution models of viruses. Particle detection is based on the comparison of intensity in a circular area and in the surrounding ring followed by a number of tests to validate the potential particles. The only required input from the user in addition to the micrograph is an approximate radius of the particle. The method has been implemented as program ETHAN that has been tested for several different data sets. ETHAN has also successfully been used to detect DNA-less virus particles for an actual reconstruction."
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
