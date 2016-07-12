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
List of related references in Bibtex format for dosefgpu program
developed by Xueming Li at Yifan Cheng lab and
motioncor2 program developed by Shawn Zheng at David Agard lab.
"""

_bibtexStr = """

@article{Li2013,
  title="Electron counting and beam-induced motion correction enable near-atomic-resolution single-particle cryo-EM",
  author="Li, Xueming and Mooney, Paul and Zheng, Shawn and Booth, Christopher R and Braunfeld, Michael B and Gubbens, Sander and Agard, David A and Cheng, Yifan",
  journal="Nature methods",
  volume="10",
  number="6",
  pages="584-590",
  year="2013",
  publisher="Nature Publishing Group",
  doi = "http://dx.doi.org/10.1038/nmeth.2727"
}

@article{Zheng2016,
  title="Anisotropic Correction of Beam-induced Motion for Improved Single-particle Electron Cryo-microscopy",
  author="Zheng, Shawn and Palovcak, Eugene and Armache, Jean-Paule and Cheng, Yifan and Agard, David",
  year="2016",
  journal="submitted",
  doi = "http://dx.doi.org/10.1101/061960"
}

"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
