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
This modules serve to define some Configuration classes
mainly for project GUI
"""

import json
import os
import sys
from ConfigParser import ConfigParser
from collections import OrderedDict

import pyworkflow as pw
import pyworkflow.object as pwobj


class DownloadRecord(pwobj.OrderedObject):
    """ Store information about Scipion downloads. """
    def __init__(self, **kwargs):
        pwobj.OrderedObject.__init__(self, **kwargs)
        
        self.fullName = pwobj.String(kwargs.get('fullName', None))
        self.organization = pwobj.String(kwargs.get('organization', None))
        self.email = pwobj.String(kwargs.get('email', None))
        self.subscription = pwobj.String(kwargs.get('subscription', None))
        self.country = pwobj.String(kwargs.get('country', None))
        self.version = pwobj.String(kwargs.get('version', None))
        self.platform = pwobj.String(kwargs.get('platform', None))

