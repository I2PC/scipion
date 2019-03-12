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
import ast
import os
from .utils.path import findFile


# Versions
VERSION_1 = 'v1.0'
VERSION_1_1 = 'v1.1'
VERSION_1_2 = 'v1.2'
VERSION_2_0 = 'v2.0'
LAST_VERSION = VERSION_2_0
OLD_VERSIONS = (VERSION_1, VERSION_1_1)


HOME = os.path.abspath(os.path.dirname(__file__))

PYTHON = os.environ.get("SCIPION_PYTHON", 'python')


class Config:
    __get = os.environ.get  # shortcut
    SCIPION_HOME = __get('SCIPION_HOME', '')
    SCIPION_USER_DATA = __get('SCIPION_USER_DATA', '')
    SCIPION_SUPPORT_EMAIL = __get('SCIPION_SUPPORT_EMAIL',
                                  'scipion@cnb.csic.es')
    SCIPION_LOGO = __get('SCIPION_LOGO',
                         'scipion_logo.png')
    # Where is the input data for tests...also where it will be downloaded
    SCIPION_TESTS = __get('SCIPION_TESTS',
                          os.path.join(SCIPION_HOME, 'data', 'tests'))

    # Where the output of the tests will be stored
    SCIPION_TESTS_OUTPUT = __get('SCIPION_TESTS_OUTPUT',
                                 os.path.join(SCIPION_USER_DATA, 'Tests'))

    SCIPION_CONFIG_MAIN = __get('SCIPION_CONFIG_MAIN', 'scipion.conf')
    SCIPION_CONFIG_HOSTS = __get('SCIPION_CONFIG_HOSTS', 'hosts.conf')
    SCIPION_CONFIG_PROTOCOLS = __get('SCIPION_CONFIG_PROTOCOLS', 'protocols.conf')

    SCIPION_PLUGIN_JSON = __get('SCIPION_PLUGIN_JSON', None)
    SCIPION_PLUGIN_REPO_URL = __get('SCIPION_PLUGIN_REPO_URL',
                                    'http://scipion.i2pc.es/getplugins/')

    # Get general log file path
    LOG_FILE = os.path.join(__get('SCIPION_LOGS'), 'scipion.log')

    try:
        VIEWERS = ast.literal_eval(__get('VIEWERS', "{}"))
    except Exception as e:
        VIEWERS = {}
        print("ERROR loading preferred viewers, VIEWERS variable will be ignored")
        print(e)

def join(*paths):
    """ join paths from HOME . """
    return os.path.join(HOME, *paths)


__resourcesPath = [join('resources')]


def findResource(filename):
    from utils.path import findFile

    return findFile(filename, *__resourcesPath)

# Following are a set of functions to centralize the way to get
# files from several scipion folder such as: config or apps

def getScipionPath(*paths):
     return os.path.join(Config.SCIPION_HOME, *paths)


def getScipionScript():
    return getScipionPath('scipion')


def getConfigPath(*paths):
    return getScipionPath('config', *paths)


def getTemplatePath(*paths):
    return join('templates', *paths)
