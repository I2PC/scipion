# coding: latin-1
# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
title = {Gctf: Real-time CTF determination and correction},
abstract = {Accurate estimation of the contrast transfer function (CTF) is critical for a near-atomic resolution cryo electron microscopy (cryoEM) reconstruction. Here, a GPU-accelerated computer program, Gctf, for accurate and robust, real-time CTF determination is presented. The main target of Gctf is to maximize the cross-correlation of a simulated CTF with the logarithmic amplitude spectra (LAS) of observed micrographs after background subtraction. Novel approaches in Gctf improve both speed and accuracy. In addition to GPU acceleration (e.g. 10-50x), a fast "1-dimensional search plus 2-dimensional refinement (1S2R)" procedure further speeds up Gctf. Based on the global CTF determination, the local defocus for each particle and for single frames of movies is accurately refined, which improves CTF parameters of all particles for subsequent image processing. Novel diagnosis method using equiphase averaging (EPA) and self-consistency verification procedures have also been implemented in the program for practical use, especially for aims of near-atomic reconstruction. Gctf is an independent program and the outputs can be easily imported into other cryoEM software such as Relion (Scheres, 2012) and Frealign (Grigorieff, 2007). The results from several representative datasets are shown and discussed in this paper.},
journal = {JSB},
volume = {},
number = {},
pages = {},
year = {2015},
note = {in press},
issn = {},
doi = {http://dx.doi.org/10.1016/j.jsb.2015.11.003},
url = {http://www.sciencedirect.com/science/article/pii/S1047847715301003},
author = {Zhang, Kai},
keywords = {Contrast transfer function, Cryo-electron microscopy, GPU program, CTF determination}
}

"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
