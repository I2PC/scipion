# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
Some Imagic protocol base classes.
"""

from pyworkflow.em import EMProtocol
from ..imagic import runTemplate


class ImagicProtocol(EMProtocol):
    """ Sub-class of EMProtocol to group some common Imagic utils. """

    def _getFileName(self, key):
        """ Give a key, append the img extension
        and prefix the protocol working dir.
        """
        template = '%(' + key + ')s' + '.img'

        return self._getPath(template % self._params)

    def runTemplate(self, inputScript, paramsDict):
        """ This function will create a valid Imagic script
        by copying the template and replacing the values in dictionary.
        After the new file is read, the Imagic interpreter is invoked.
        """
        self._enterWorkingDir()

        log = getattr(self, '_log', None)
        runTemplate(inputScript, paramsDict, log)

        self._leaveWorkingDir()

        f = open(self.getLogPaths()[0], 'r')
        for line in f.readlines():
            if '**ERROR' in line:
                raise Exception('IMAGIC script error!')
        f.close()
