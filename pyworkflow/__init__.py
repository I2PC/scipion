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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
import os.path

HOME = os.path.abspath(os.path.dirname(__file__))
RESOURCES = [os.path.join(HOME, 'resources')]

if "SCIPION_HOME" not in os.environ:
    raise Exception("SCIPION_HOME is not defined as environment variable")

SCIPION_HOME = os.environ["SCIPION_HOME"]

PROJECTS = os.path.join(SCIPION_HOME, 'projects')
SETTINGS = os.path.join(SCIPION_HOME, 'settings.sqlite')
TESTS = os.path.join(SCIPION_HOME, 'tests')

from utils.path import findResource

