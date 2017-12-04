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

#   -+-+-+- C H A N G E   T H I S   S T R I N G -+-+-+-
_bibtexStr = """
@article {Jakobi2017,
article_type 	= {journal},
title 			= {Model-based local density sharpening of cryo-EM maps},
author 			= {Jakobi, Arjen J and Wilmanns, Matthias and Sachse, Carsten},
editor 			= {Brunger, Axel T},
volume 			= 6,
year 			= 2017,
month 			= {oct},
pub_date 		= {2017-10-23},
pages 			= {e27131},
citation 		= {eLife 2017;6:e27131},
doi 			= {10.7554/eLife.27131},
url 			= {https://doi.org/10.7554/eLife.27131},
keywords 		= {cryo-EM, model building, B-factor sharpening, contrast improvement, amplitude scaling},
journal 		= {eLife},
issn 			= {2050-084X},
publisher 		= {eLife Sciences Publications, Ltd},
}
"""


from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  


