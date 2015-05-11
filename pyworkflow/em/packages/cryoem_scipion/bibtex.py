# coding: latin-1
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for Cryoem package.
"""

_bibtexStr = """

@Article{Joubert2015,
  Title                    = {Bayesian Inference of Initial Models in Cryo-Electron Microscopy Using Pseudo-atoms.},
  Author                   = {Joubert, P. and Habeck, M.},
  Journal                  = {Biophysical Journal},
  Year                     = {2015},
  Month                    = {March},
  Number                   = {5},
  Pages                    = {1165--1175},
  Volume                   = {108},
  Doi                      = {http://dx.doi.org/10.1016/j.bpj.2014.12.054},
  Keywords                 = {Algorithms; Cryoelectron Microscopy, methods; Pseudo atoms, methods; Ribosomes;},
  Pii                      = {S0006349515000648},
  Pmid                     = {25762328},
  Url                      = {http://dx.doi.org/10.1016/j.bpj.2014.12.054}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
