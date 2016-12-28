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

import os

HOME = os.path.abspath(os.path.dirname(__file__))

def join(*paths):
    """ join paths from HOME . """
    return os.path.join(HOME, *paths)

RESOURCES = [join('resources')]
WEB_RESOURCES = join('web', 'pages', 'resources')

if "SCIPION_USER_DATA" not in os.environ:
    raise Exception("SCIPION_USER_DATA is not defined as environment variable")

SCIPION_HOME = os.environ['SCIPION_HOME']
SCIPION_USER_DATA = os.environ["SCIPION_USER_DATA"]
PYTHON = os.environ.get("SCIPION_PYTHON", 'python')
SCIPION_PYTHON = PYTHON



SETTINGS = os.path.join(SCIPION_USER_DATA, 'settings.sqlite')


_logo = 'scipion_logo.png'


#----------------- Some functions to centralize path access -------------------

def findResource(filename):
    """ This function will search for a give
    resource filename in the paths specified
    in pyworkflow.RESOURCES path list.
    """
    from utils.path import findFile

    return findFile(filename, *RESOURCES)

# Following are a set of functions to centralize the way to get
# files from several scipion folder such as: config or apps

def getScipionPath(*paths):
     return os.path.join(SCIPION_HOME, *paths)


def getScipionScript():
    return getScipionPath('scipion')


def getConfigPath(*paths):
    return getScipionPath('config', *paths)


def getTemplatePath(*paths):
    return getConfigPath('templates', *paths)
