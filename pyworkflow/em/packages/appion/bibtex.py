# coding: latin-1
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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
# *  e-mail address 'jlvilas@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for DoGPiker.
"""

_bibtexStr = """

@Article{Voss2009,
  Title                    = {BDoG Picker and TiltPicker: software tools to facilitate particle selection in single particle electron microscopy.},
  Author                   = {Voss, N.R. Yoshioka, C.K. Radermacher, M. Potter, C.S. and Carragher, B.},
  Journal                  = {Journal Structutal Biology},
  Year                     = {2009},
  Month                    = {May},
  Number                   = {2},
  Pages                    = {205--213},
  Volume                   = {166},
  Doi                      = {http://dx.doi.org/10.1016/j.jsb.2009.01.004},
  Keywords                 = {EM; Cryo-electron microscopy; Particle picking; Random conical tilt; Orthogonal tilt reconstruction},
  Pii                      = {S1047847709000197},
  Pmid                     = {19374019},
  Url                      = {http://dx.doi.org/10.1016/j.jsb.2009.01.004}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
