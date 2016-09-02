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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************


from pyworkflow.protocol import (STATUS_SAVED, STATUS_LAUNCHED, STATUS_RUNNING,
                                 STATUS_FINISHED, STATUS_FAILED, STATUS_INTERACTIVE,
                                 STATUS_ABORTED)

STATUS_COLORS = {
               STATUS_SAVED: '#D9F1FA',
               STATUS_LAUNCHED: '#D9F1FA',
               STATUS_RUNNING: '#FCCE62',
               STATUS_FINISHED: '#D2F5CB',
               STATUS_FAILED: '#F5CCCB',
               STATUS_INTERACTIVE: '#F3F5CB',
               STATUS_ABORTED: '#F5CCCB',
               #STATUS_SAVED: '#124EB0',
               }

# Based on: http://paletton.com/palette.php?uid=75C1d0kleqtbzEKgVuIpcmGtdhZ
# STATUS_COLORS = {
#     STATUS_SAVED: '#424093',
#     STATUS_LAUNCHED: '#424093',
#     STATUS_RUNNING: '#d3b147',
#     STATUS_FINISHED: '#3fab3a',
#     STATUS_FAILED: '#d1464a',
#     STATUS_INTERACTIVE: '#f5d573',
#     STATUS_ABORTED: '#d1464a',
#     # STATUS_SAVED: '#124EB0',
# }