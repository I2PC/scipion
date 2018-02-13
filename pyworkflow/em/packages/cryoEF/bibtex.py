# coding: latin-1
# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
Bibtex string file for cryoEF package.
"""

_bibtexStr = """

@article{naydenova2017,
  title={Measuring the effects of particle orientation to improve the efficiency of electron cryomicroscopy},
  author={Naydenova, K. and Russo, C.J.},
  journal={Nature communications},
  volume={8},
  number={1},
  pages={},
  year={2017},
  publisher={Nature Publishing Group},
  doi = {http://dx.doi.org/10.1038/s41467-017-00782-3}
}


"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
