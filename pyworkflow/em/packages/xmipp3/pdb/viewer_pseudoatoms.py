# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Nov 2014
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

from pyworkflow.em.viewer import DataView, ChimeraView
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer

from protocol_pseudoatoms import XmippProtConvertToPseudoAtoms



class XmippPseudoAtomsViewer(XmippViewer):
    """ Visualize the output of protocol Convert to PseudoAtoms """
    _label = 'pseudoatoms viewer'
    _targets = [XmippProtConvertToPseudoAtoms]
    
    def _visualize(self, obj, **args):
        self._views.append(ChimeraView(obj.outputPdb._chimeraScript))
        self._views.append(DataView(obj._getExtraPath('pseudoatoms_approximation.mrc')))

