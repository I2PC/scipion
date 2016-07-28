# coding: latin-1
# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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

@Article{Serban2015,
  Title                    = {Localized reconstruction of subunits from electron cryomicroscopy images.},
  Author                   = {Ilca, Serban L. and Kotecha, A. Sun, Xiaoyu and Poranen, Minna M. and Stuart, David I. and Huiskonen, Juha T.},
  Journal                  = {Nat Commun},
  Year                     = {2015},
  Month                    = {November},
  Volume                   = {6},
  Doi                      = {http://dx.doi.org/10.1038/ncomms9843},
  Url                      = {http://www.nature.com/ncomms/2015/151104/ncomms9843/full/ncomms9843.html}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
