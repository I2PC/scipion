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
@Article{Pettersen2004,
   Author="Pettersen, E. F.  and Goddard, T. D.  and Huang, C. C.  and Couch, G. S.  and Greenblatt, D. M.  and Meng, E. C.  and Ferrin, T. E. ",
   Title="{{U}{C}{S}{F} {C}himera--a visualization system for exploratory research and analysis}",
   Journal="J Comput Chem",
   Year="2004",
   Volume="25",
   Number="13",
   Pages="1605--1612",
   Month="Oct",
   doi="https://doi.org/10.1002/jcc.20084"
}
"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
