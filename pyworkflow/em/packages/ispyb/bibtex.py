# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es) [1]
# *              Kevin Savage (kevin.savage@diamond.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] Diamond Light Source, Ltd
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

@article{Delageniere2011,
  author    = {Delageniere, Solange and
               Brenchereau, Patrice  and
               Launer, Ludovic  and
               Ashton, Alun W.  and
               Leal, Ricardo  and
               Veyrier, Stephanie  and
               Gabadinho, Jose  and
               Gordon, Elspeth J.  and
               Jones, Samuel D.  and
               Levik, Karl Erik  and
               McSweeney, Sean M.  and
               Monaco, Stephanie  and
               Nanao, Max  and
               Spruce, Darren  and
               Svensson, Olof  and
               Walsh, Martin A.  and
               Leonard, Gordon A.},
  title     = {ISPyB: an information management system for synchrotron macromolecular
               crystallography},
  journal   = {Bioinformatics},
  volume    = {27},
  number    = {22},
  pages     = {3186--3192},
  year      = {2011},
  url       = {http://dx.doi.org/10.1093/bioinformatics/btr535},
  doi       = {http://dx.doi.org/10.1093/bioinformatics/btr535},
  timestamp = {Mon, 14 Nov 2011 15:11:21 +0100}
}
"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  

'''

'''