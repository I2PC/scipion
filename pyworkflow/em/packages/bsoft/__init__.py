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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains the protocols and data for BSOFT
"""

from bibtex import _bibtex # Load bibtex dict with references
 
_logo = "bsoft_logo.png"
_references = ['Heymann2007']


from protocol_particle_pick import BsoftProtParticlePicking
from protocol_bfilter import BsoftProtBfilter
from protocol_blocres import BsoftProtBlocres

from wizard import BsoftFilterParticlesWizard


from convert import getEnviron, getVersion
_environ = getEnviron()

# Since bsoft is not installed by default,
# remove for now the viewer
from viewer import BsoftViewerBlocres
