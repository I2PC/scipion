# coding: latin-1
# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for Gautomatch package.
"""

_bibtexStr = """
@Article{Emsley_2004,
Author="Emsley, P.  and Cowtan, K. ",
Title="{{C}oot: model-building tools for molecular graphics}",
Journal="Acta Crystallogr. D Biol. Crystallogr.",
Year="2004",
Volume="60",
Number="Pt 12 Pt 1",
Pages="2126--2132",
Month="Dec",
doi = "http://doi.org/10.1107/S0907444904019158",
url = "http://scripts.iucr.org/cgi-bin/paper?S0907444904019158"
}
"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
