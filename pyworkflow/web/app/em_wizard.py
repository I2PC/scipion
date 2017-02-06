# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

# Imports for web wizards
from pyworkflow.web.app.wizards.xmipp_wizard import *
from pyworkflow.web.app.wizards.spider_wizard import *
from pyworkflow.web.app.wizards.relion_wizard import *
from pyworkflow.web.app.wizards.brandeis_wizard import *
from pyworkflow.web.app.wizards.resmap_wizard import *

#===============================================================================
# Wizard base function (to call the others)
#===============================================================================

def wizard(request):
    from views_util import loadProtocolProject
    from views_protocol import updateProtocolParams
    
    # Get the post-dictionary
    requestDict = getattr(request, "POST")

    # Get and instance the wizard class
    className = requestDict.get("wizClassName")
    wizClass = globals().get(className, None)()

    # Get the protocol object
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)

    print "======================= in wizard: " + str(wizClass)
    
    # Obtain the parameters for the wizard
    return wizClass._run(protocol, request)

