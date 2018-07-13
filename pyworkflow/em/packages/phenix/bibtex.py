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
   Author="Adams, P. D.  and Afonine, P. V.  and Bunkoczi, G.  and Chen, V. B.  and Davis, I. W.  and Echols, N.  and Headd, J. J.  and Hung, L. W.  and Kapral, G. J.  and Grosse-Kunstleve, R. W.  and McCoy, A. J.  and Moriarty, N. W.  and Oeffner, R.  and Read, R. J.  and Richardson, D. C.  and Richardson, J. S.  and Terwilliger, T. C.  and Zwart, P. H. ",
   Title="{{P}{H}{E}{N}{I}{X}: a comprehensive {P}ython-based system for macromolecular structure solution}",
   Journal="Acta Crystallogr. D Biol. Crystallogr.",
   Year="2010",
   Volume="66",
   Number="Pt 2",
   Pages="213--221",
   Month="Feb",
   doi = "http://doi.org/10.1107/S0907444909052925",
   url = "http://scripts.iucr.org/cgi-bin/paper?dz5186"
}


@Article{Barad_2015,
   Author="Barad, B. A.  and Echols, N.  and Wang, R. Y.  and Cheng, Y.  and DiMaio, F.  and Adams, P. D.  and Fraser, J. S. ",
   Title="{{E}{M}{R}inger: side chain-directed model and map validation for 3{D} cryo-electron microscopy}",
   Journal="Nat. Methods",
   Year="2015",
   Volume="12",
   Number="10",
   Pages="943--946",
   Month="Oct",
   doi = "http://doi.org/10.1038/nmeth.3541",
   url = "http://www.nature.com/articles/nmeth.3541"
}

@Article{Chen_2010,
   Author="Chen, V. B.  and Arendall, W. B.  and Headd, J. J.  and Keedy, D. A.  and Immormino, R. M.  and Kapral, G. J.  and Murray, L. W.  and Richardson, J. S.  and Richardson, D. C. ",
   Title="{{M}ol{P}robity: all-atom structure validation for macromolecular crystallography}",
   Journal="Acta Crystallogr. D Biol. Crystallogr.",
   Year="2010",
   Volume="66",
   Number="Pt 1",
   Pages="12--21",
   Month="Jan",
   doi = "http://dx.doi.org/10.1016/j.jsb.2010.03.011",
   url = "http://www.sciencedirect.com/science/article/pii/S1047847710000882"

}
"""
from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)
