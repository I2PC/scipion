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
Bibtex string file for 3DFSC package.
"""

_bibtexStr = """

@article{tan2017,
  title={Addressing preferred specimen orientation in single-particle cryo-EM through tilting},
  author={Tan, Y.Z. and Baldwin, P.R. and Davis, J.H. and Williamson, J.R. and Potter, C.S. and Carragher, B. and Lyumkis, D.},
  journal={Nature methods},
  volume={14},
  number={8},
  pages={793--796},
  year={2017},
  publisher={Nature Publishing Group},
  doi = {http://dx.doi.org/10.1038/nmeth.4347}
}


"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
