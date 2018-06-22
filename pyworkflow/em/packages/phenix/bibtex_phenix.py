# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *
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
@Article{Adams_2010,
Author="Adams, P.D., Afonine, P.V., Bunkoczi, G., Chen, V.B., Davis, I.W., 
Echols, N., Headd, J.J., Hung, L.-W., Kapral, G.J., Grosse-Kunstleve, R.W., 
McCoy, A.J., Moriarty, N.W., Oeffner, R., Read, R.J., Richardson, D.C., 
Richardson, J.S., Terwilliger, T.C., and Zwart, P.H.",
Title="{{PHENIX}: a comprehensive Python-based system for macromolecular
structure solution}",
Journal="Acta Crystallogr. D Struc. Biol.",
Year="2010",
Volume="66",
Number="Pt 2",
Pages="213--221",
Month="Feb",
doi = "http://doi.org/10.1107/S0907444909052925",
url = "https://scripts.iucr.org/cgi-bin/paper?dz5186"
}

@Article{Barad_2015,
Author="Barad, B.A., Echols, N., Wang, R.Y., Cheng, Y., DiMaio, F., 
Adams, P.D., and Fraser, J.S.",
Title="{{EMR}inger: side chain-directed model and map validation for 3D 
cryo-electron microscopy}",
Journal="Nat. Methods",
Year="2015",
Volume="12",
Number="",
Pages="213--221",
Month="Aug",
doi = "https://doi.org/10.1038/nmeth.3541",
url = "https://www.nature.com/articles/nmeth.3541"
}

@Article{Chen_2010,
Author="Chen, V.B., Arendall, W.B., Headd, J.J., Keedy, D.A., Immormino, R.M.,
Kapral, G.J., Murray, L.W., Richardson, J.S., and Richardson, D.C.",
Title="{{M}ol{P}robity: all-atom structure validation for macromolecular 
crystallography}",
Journal="Acta Crystallogr. D Biol. Crystallogr.",
Year="2010",
Volume="66",
Number="Pt 1",
Pages="12--21",
Month="Dec",
doi = "http://doi:10.1107/S0907444909042073",
url = "https://scripts.iucr.org/cgi-bin/paper?S0907444909042073"
}
"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)
