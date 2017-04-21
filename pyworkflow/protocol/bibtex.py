# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Yaiza Rancel (cyrancel@cnb.csic.es)
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
"""
Bibtex string file for Scipion package.
"""

_bibtexStr = """

@article{delaRosaTrevin201693,
title = "Scipion: A software framework toward integration, reproducibility and validation in 3D electron microscopy ",
journal = "Journal of Structural Biology",
volume = "195",
number = "1",
pages = "93 - 99",
year = "2016",
note = "",
issn = "1047-8477",
doi = "http://doi.org/10.1016/j.jsb.2016.04.010",
url = "http://www.sciencedirect.com/science/article/pii/S104784771630079X",
author = "J.M. de la Rosa-Trevín and A. Quintana and L. del Cano and A. Zaldívar and I. Foche and J. Gutiérrez and J. Gómez-Blanco and J. Burguet-Castell and J. Cuenca-Alba and V. Abrishami and J. Vargas and J. Otón and G. Sharov and J.L. Vilas and J. Navas and P. Conesa and M. Kazemi and R. Marabini and C.O.S. Sorzano and J.M. Carazo",
keywords = "Electron microscopy",
keywords = "Single particle analysis",
keywords = "Image processing",
keywords = "Software package",
keywords = "Workflows",
keywords = "Reproducibility "
}
"""



from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
