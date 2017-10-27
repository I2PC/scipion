# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
from pyworkflow.gui.project.constants import STATUS_COLORS
from pyworkflow.protocol import STATUS_FAILED


def getStatusColorFromNode(node):

    # If it is a run node (not PROJECT)
    return getStatusColorFromRun(node.run)


def getStatusColorFromRun(prot):
    """
    Returns the color associated with the status
    -------

    """

    if prot:
        return getStatusColor(prot.status.get(STATUS_FAILED))
    else:
        return getStatusColor()


def getStatusColor(status=None, default='#ADD8E6'):
    """
    Parameters
    ----------
    status status of the protocol

    Returns the color associated with he status
    -------

    """
    if status:
        return STATUS_COLORS[status]
    else:
        return default