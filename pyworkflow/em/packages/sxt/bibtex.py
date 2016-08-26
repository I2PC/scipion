# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for SXT.
"""

_bibtexStr = """

@article{J.Oton2015,
title = "Measurement of the modulation transfer function of an X-ray microscope based on multiple Fourier orders analysis of a Siemens star",
journal = "Optics Express",
volume = "23",
number = "8",
pages = "",
year = "2015",
note = "",
issn = "",
doi = "10.1364/OE.23.009567",
url = "http://www.ncbi.nlm.nih.gov/pubmed/25968993",
author = "J. Oton, C.O.S. Sorzano, R. Marabini, Eva Pereiro. J.M. Carazo",
keywords = "",
keywords = "",
keywords = "",
keywords = ""
}

@article{J.Oton2016,
title = "XTEND: Extending the depth of field in soft X-ray tomography",
journal = "Journal",
volume = "",
number = "",
pages = "",
year = "2016",
note = "",
issn = "",
doi = "",
url = "http://www.biocomp.cnb.csic.es/publications",
author = "J. Oton, E. Pereiro, J.J. Conesa, F.J. Chichon, D. Luque, J.M. Rodriguez, A.J. Perez-Berna, C.O.S. Sorzano, J. Klukowska, G.T. Herman, J. Vargas, R. Marabini, J.L. Carrascosa, J.M. Carazo",
keywords = "",
keywords = "",
keywords = "",
keywords = ""
}

"""


from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
