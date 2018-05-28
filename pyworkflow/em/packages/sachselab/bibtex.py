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

_bibtexStr = """

@Article{Jakobi2017,
  Title                    = {Model-based local density sharpening of cryo-EM maps},
  Author                   = {Jakobi, Arjen J and Wilmanns, Matthias and Sachse, Carsten},
  Journal                  = {eLife},
  Year                     = {2017},
  Month                    = {October},
  Volume                   = {6},
  Doi                      = {http://doi.org/10.7554/eLife.27131},
  Url                      = {http://elifesciences.org/articles/27131},
  Citation                 = {eLife 2017;6:e27131},
  Issn                     = {2050-084X},
  Publisher                = {eLife Sciences Publications, Ltd}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  


