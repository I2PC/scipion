# coding: latin-1
# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
Bibtex string file for Imagic protocols.
"""

from pyworkflow.utils import parseBibTex

_bibtexStr = """

@article{Borland1990,
  author = {Borland, L. and van Heel, M.},
  title = {Classification of image data in conjugate representation spaces},
  journal = {J.Opt.Soc.Am. A},
  year = {1990},
  volume = {7},
  issue = {4},
  pages = {601--610},
  doi = {http://dx.doi.org/10.1364/JOSAA.7.000601}
}

@article{vanHeel1981,
  author = {Marin van Heel and Wilko Keegstra},
  title = {IMAGIC: A fast, flexible and friendly image analysis software system},
  journal = {Ultramicroscopy},
  year = {1981},
  volume = {7},
  issue = {2},
  pages = {113--129},
  doi = {http://dx.doi.org/10.1016/0304-3991(81)90001-2}
}

@article{vanHeel1984,
  author = {Marin van Heel},
  title = {Multivariate statistical classification of noisy images (randomly oriented biological macromolecules)},
  journal = {Ultramicroscopy},
  year = {1984},
  volume = {13},
  issue = {1-2},
  pages = {165--183},
  doi = {http://dx.doi.org/10.1016/0304-3991(84)90066-4}
}

@article{vanHeel1989,
  author = {Marin van Heel},
  title = {Classification of very large electron microscopical image data sets},
  journal = {Optik},
  year = {1989},
  volume = {82},
  pages = {114--126},
  doi = {}
}

@article{vanHeel1996,
  author = {Marin van Heel, George Harauz, Elena V. Orlova},
  title = {A New Generation of the IMAGIC Image Processing System},
  journal = {J.Str.Biol},
  year = {1996},
  volume = {116},
  issue = {1}
  pages = {17-24},
  doi = {http://dx.doi.org/10.1006/jsbi.1996.0004}
}

@article{vanHeel2012,
  author = {M. van Heel, R. Portugal, A. Rohou, C. Linnemayr, C. Bebeacua, R. Schmidt, T. Grant and M. Schatz},
  title = {Four-dimensional cryo-electron microscopy at quasi-atomic resolution: IMAGIC 4D},
  journal = {International Tables for Crystallography},
  year = {2012},
  volume = {F},
  chapter = {19.9}
  pages = {624-628},
  doi = {http://dx.doi.org/10.1107/97809553602060000875}
}

"""

_bibtex = parseBibTex(_bibtexStr)
