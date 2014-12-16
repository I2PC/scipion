# coding: latin-1
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for Simple package.
"""

_bibtexStr = """

@Article{Elmlund2010,
  Title                    = {Ab initio structure determination from electron microscopic images of single molecules coexisting in different functional states.},
  Author                   = {Elmlund, D. and Davis, R. and Elmlund, H.},
  Journal                  = {Structure},
  Year                     = {2010},
  Month                    = {Jul},
  Number                   = {7},
  Pages                    = {777--786},
  Volume                   = {18},
  Doi                      = {http://dx.doi.org/10.1016/j.str.2010.06.001},
  Keywords                 = {Algorithms; Cryoelectron Microscopy, methods; Escherichia coli; Fourier Analysis; Image Processing, Computer-Assisted, methods; Models, Molecular; Nanoparticles, chemistry; RNA Polymerase II, chemistry; Ribosomes, ultrastructure; Yeasts},
  Pii                      = {S0969-2126(10)00192-9},
  Pmid                     = {20637414},
  Url                      = {http://dx.doi.org/10.1016/j.str.2010.06.001}
}

@Article{Elmlund2012,
  Title                    = {SIMPLE: Software for ab initio reconstruction of heterogeneous single-particles.},
  Author                   = {Elmlund, D. and Elmlund, H.},
  Journal                  = {J Struct Biol},
  Year                     = {2012},
  Month                    = {Dec},
  Number                   = {3},
  Pages                    = {420--427},
  Volume                   = {180},
  Doi                      = {http://dx.doi.org/10.1016/j.jsb.2012.07.010},
  Keywords                 = {Algorithms; Computer Simulation; Cryoelectron Microscopy; Fourier Analysis; Image Processing, Computer-Assisted; Imaging, Three-Dimensional; Models, Molecular; Software},
  Pii                      = {S1047-8477(12)00218-3},
  Pmid                     = {22902564},
  Url                      = {http://dx.doi.org/10.1016/j.jsb.2012.07.010}
}

@Article{Elmlund2013,
  Title                    = {PRIME: probabilistic initial {3D} model generation for single-particle cryo-electron microscopy.},
  Author                   = {Elmlund, Hans and Elmlund, Dominika and Bengio, Samy},
  Journal                  = {Structure},
  Year                     = {2013},

  Month                    = {Aug},
  Number                   = {8},
  Pages                    = {1299--1306},
  Volume                   = {21},

  Abstract                 = {Low-dose electron microscopy of cryo-preserved individual biomolecules (single-particle cryo-EM) is a powerful tool for obtaining information about the structure and dynamics of large macromolecular assemblies. Acquiring images with low dose reduces radiation damage, preserves atomic structural details, but results in low signal-to-noise ratio of the individual images. The projection directions of the two-dimensional images are random and unknown. The grand challenge is to achieve the precise three-dimensional (3D) alignment of many (tens of thousands to millions) noisy projection images, which may then be combined to obtain a faithful 3D map. An accurate initial 3D model is critical for obtaining the precise 3D alignment required for high-resolution (<10 Å) map reconstruction. We report a method (PRIME) that, in a single step and without prior structural knowledge, can generate an accurate initial 3D map directly from the noisy images.},
  Doi                      = {http://dx.doi.org/10.1016/j.str.2013.07.002},
  Institution              = {Department of Structural Biology, Stanford University Medical School, Fairchild Building, 1st Floor, 299 Campus Drive, Stanford, CA 94305-5126, USA. hael@stanford.edu},
  Keywords                 = {Cryoelectron Microscopy, methods; Imaging, Three-Dimensional, methods; Macromolecular Substances, ultrastructure; Models, Molecular; Models, Statistical; Ribosomes, ultrastructure; Software},
  Pii                      = {S0969-2126(13)00250-5},
  Pmid                     = {23931142},
  Timestamp                = {2014.06.26},
  Url                      = {http://dx.doi.org/10.1016/j.str.2013.07.002}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
