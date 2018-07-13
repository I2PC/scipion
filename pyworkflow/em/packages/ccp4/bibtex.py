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
@article{Winn_2011,
   Author="Winn, M. D.  and Ballard, C. C.  and Cowtan, K. D.  and Dodson, E. J.  and Emsley, P.  and Evans, P. R.  and Keegan, R. M.  and Krissinel, E. B.  and Leslie, A. G.  and McCoy, A.  and McNicholas, S. J.  and Murshudov, G. N.  and Pannu, N. S.  and Potterton, E. A.  and Powell, H. R.  and Read, R. J.  and Vagin, A.  and Wilson, K. S. ",
   Title="{{O}verview of the {C}{C}{P}4 suite and current developments}",
   Journal="Acta Crystallogr. D Biol. Crystallogr.",
   Year="2011",
   Volume="67",
   Number="Pt 4",
   Pages="235--242",
   Month="Apr",
   doi = "http://doi.org/10.1107/S0907444910045749",
   url = "http://scripts.iucr.org/cgi-bin/paper?S0907444910045749" 
}
@article{Emsley_2004,
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
@article{Vagin_2004,
   Author="Vagin, A. A.  and Steiner, R. A.  and Lebedev, A. A.  and Potterton, L.  and McNicholas, S.  and Long, F.  and Murshudov, G. N. ",
   Title="{{R}{E}{F}{M}{A}{C}5 dictionary: organization of prior chemical knowledge and guidelines for its use}",
   Journal="Acta Crystallogr. D Biol. Crystallogr.",
   Year="2004",
   Volume="60",
   Number="Pt 12 Pt 1",
   Pages="2184--2195",
   Month="Dec",
   doi = "http://doi.org/10.1107/S0907444904023510",
   url = "http://scripts.iucr.org/cgi-bin/paper?ba5073"
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
