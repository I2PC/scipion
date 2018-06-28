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
This modules handles Plugin management
"""

import os
import pkg_resources
from email import message_from_string
from collections import OrderedDict

from pyworkflow.utils import Environ


class Plugin(object):

    def __init__(self, name, version=None, configVars=None,
                 logo=None, bibtex=None):
        # Plugin name - only mandatory attribute
        self.name = name
        # Plugin version
        self.version = version
        # dict with default values for env vars needed
        self.configVars = configVars
        # path to the logo, relative to each plugin's plugin.py file
        self.logo = logo
        # List with the default plugin references e.g. []
        self.bibtex = bibtex
        # Set default env vars
        self.setDefaultEnviron()
        # plugin pip metadata
        self.metadata = self.getMetadata()

    def setDefaultEnviron(self):
        for k in self.configVars:
            os.environ.setdefault(k, self.configVars[k])
        environ = Environ(os.environ)
        return environ

    def registerPluginBinaries(self, env):
        """Overwrite in subclass"""
        pass

    def getMetadata(self):
        pipPackage = pkg_resources.get_distribution(self.name)
        metadataLines = [l for l in pipPackage._get_metadata(pipPackage.PKG_INFO)]
        metadataTuples = message_from_string('\n'.join(metadataLines))
        metadata = OrderedDict()
        for v in metadataTuples.items():
            if v[0] == 'Keywords':
                break
            metadata[v[0]] = v[1]
        return metadata

    def getAuthor(self):
        return self.metadata.get('Author', "")

    def getAuthorEmail(self):
        return self.metadata.get('Author-email', '')

    def getHomePage(self):
        return self.metadata.get('Home-page', '')

    def getKeywords(self):
        return self.metadata.get('Keywords', '')
