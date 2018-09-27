# **************************************************************************
# *
# * Authors:     David Maluenda Niubo (dmaluenda@cnb.csic.es)
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

import traceback

# FIXME: This is a bypass until xmipp is replaced to xmippLib
# FIXME:   in all imports related to the C++ binding

print("\n >>> WARNING: 'import xmipp' is deprecated for the xmipp binding.")
stackList = traceback.extract_stack()
for stackLine in stackList:
    if 'import xmipp' in stackLine[3]:
        print("  > Please change 'import xmipp' to 'import xmippLib' in\n"
              "    File: '%s', line %d.\n" % (stackLine[0], stackLine[1]))
    if 'from xmipp import' in stackLine[3]:
        print("  > Please change 'from xmipp import ...' to 'from xmippLib import ...' in\n"
              "    File: '%s', line %d.\n" % (stackLine[0], stackLine[1]))

from xmippLib import *
