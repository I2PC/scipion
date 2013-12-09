# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This sub-package will contains Relion protocols
"""

_logo = "relion_logo.png"
_references = [
               'Sjors H.W. Scheres, A Bayesian View on Cryo-EM Structure Determination, Journal of Molecular Biology, Volume 415, Issue 2, 13 January 2012, Pages 406-418, ISSN 0022-2836, http://dx.doi.org/10.1016/j.jmb.2011.11.010.',
               'Sjors H.W. Scheres, RELION: Implementation of a Bayesian approach to cryo-EM structure determination, Journal of Structural Biology, Volume 180, Issue 3, December 2012, Pages 519-530, ISSN 1047-8477, http://dx.doi.org/10.1016/j.jsb.2012.09.006.'
               ]

from protocol_classify3d import Relion3DClassification