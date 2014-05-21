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
Bibtex string file for Spider protocols.
"""

_bibtexStr = """

@article{Shaikh2008,
  author = {Shaikh, T. R. and Gao, H. and Baxter, W. and Asturias, F. J. and
    Boisset, N. and Leith, A. and Frank, J.},
  title = {SPIDER image processing for single-particle reconstruction of biological
    macromolecules from electron micrographs},
  journal = {Nature Protocols},
  year = {2008},
  volume = {3},
  pages = {1941--1974},
  doi = {http://dx.doi.org/10.1038/nprot.2008.156}
}

@article{Baxter2007,
  author = {Baxter, W.T. and Leith, A. and Joachim Frank},
  title = {SPIRE: the SPIDER reconstruction engine.},
  journal = JSB,
  year = {2007},
  volume = {157},
  pages = {56--63},
  doi = {10.1016/j.jsb.2006.07.019},
  keywords = {Computational Biology; Image Processing, Computer-Assisted; Software; Software Design},
  owner = {coss},
  pii = {S1047-8477(06)00235-8},
  pmid = {17055743},
  timestamp = {2009.11.04},
  doi = {http://dx.doi.org/10.1016/j.jsb.2006.07.019}
}

@article{Frank1996b,
  author = {J. Frank and M. Radermacher and P. Penczek and J. Zhu and Y. Li and
    et.al.},
  title = {{SPIDER} and {WEB}: {P}rocessing and visualization of images in 3{D}
    electron microscopy and related fields.},
  journal = JSB,
  year = {1996},
  volume = {116},
  pages = {190-9},
  doi = {http://dx.doi.org/10.1006/jsbi.1996.0030}
}

@article{Marco1996,
  author = {Marco, S. and Chagoyen, M. and {de la Fraga}, {L. G.} and Carazo,
    {J. M.} and Carrascosa, {J. L.}},
  title = {A variant to the "random approximation" of the reference-free alignment
    algorithm},
  journal = Ultramicroscopy,
  year = {1996},
  volume = {66},
  pages = {5-10},
  abstract = {Angular refinement. 2D alignment},
  doi = {http://dx.doi.org/10.1016/S0304-3991(96)00083-6}
}


"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
