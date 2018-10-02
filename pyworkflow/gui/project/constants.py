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


from pyworkflow.protocol import (STATUS_SAVED, STATUS_LAUNCHED, STATUS_RUNNING,
                                 STATUS_FINISHED, STATUS_FAILED,
                                 STATUS_INTERACTIVE, STATUS_ABORTED,
                                 STATUS_SCHEDULED)

STATUS_COLORS = {
               STATUS_SAVED: '#D9F1FA',
               STATUS_LAUNCHED: '#D9F1FA',
               STATUS_RUNNING: '#FCCE62',
               STATUS_FINISHED: '#D2F5CB',
               STATUS_FAILED: '#F5CCCB',
               STATUS_INTERACTIVE: '#F3F5CB',
               STATUS_ABORTED: '#F5CCCB',
               STATUS_SCHEDULED: '#F3F5CB'
               }

# For protocols with warnings
WARNING_COLOR = '#848484'