# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Tapu Shaikh            (shaikh@ceitec.muni.cz)
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
  author = {Frank, J. and Radermacher, M. and Penczek, P. and Zhu, J. and Li, Y. and
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

@article{Penczek1992,
        abstract = {Single particles embedded in ice pose new challenges for image processing because of the intrinsically low signal-to-noise ratio of such particles in electron micrographs. We have developed new techniques that address some of these problems and have applied these techniques to electron micrographs of the Escherichia coli ribosome. Data collection and reconstruction follow the protocol of the random-conical technique of Radermacher et al. [J. Microscopy 146 (1987) 113]. A reference-free alignment algorithm has been developed to overcome the propensity of reference-based algorithms to reinforce the reference motif in very noisy situations. In addition, an iterative 3D reconstruction method based on a chi-square minimization constraint has been developed and tested. This algorithm tends to reduce the effects of the missing angular range on the reconstruction, thereby facilitating the merging of random-conical data sets obtained from differently oriented particles}, 
        number = {1}, 
        month = {Jan}, 
        year = {1992}, 
        author = {Penczek, P and Radermacher, M and Frank, J}, 
        howpublished = {Journal}, 
        journal = {Ultramicroscopy}, 
        volume = {40}, 
        address = {Wadsworth Center, New York State Department of Health, Albany 12201-0509}, 
        pages = {33-53}, 
        url = {http://view.ncbi.nlm.nih.gov/pubmed/1580010}, 
        title = {Three-dimensional reconstruction of single particles embedded in ice}, 
        doi = {http://dx.doi.org/10.1016/0304-3991(92)90233-A}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
