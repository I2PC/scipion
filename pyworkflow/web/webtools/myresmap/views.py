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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from os.path import exists, join, basename
from pyworkflow.web.app.views_util import loadProject, getResourceCss, getResourceJs
from pyworkflow.web.app.views_base import base_grid, base_flex
from pyworkflow.web.app.views_project import contentContext
from pyworkflow.web.app.views_protocol import contextForm
from django.shortcuts import render_to_response
from pyworkflow.web.pages import settings as django_settings
from pyworkflow.manager import Manager
from django.http import HttpResponse
from pyworkflow.tests.tests import DataSet
from pyworkflow.utils import copyFile

def resmap_projects(request):
   
    if 'projectName' in request.session: request.session['projectName'] = ""
    if 'projectPath' in request.session: request.session['projectPath'] = ""

    myresmap_utils = join(django_settings.STATIC_URL, "js/", "myresmap_utils.js")

    context = {'projects_css': getResourceCss('projects'),
               'project_utils_js': getResourceJs('project_utils'),
               'myfirstmap_utils': myresmap_utils,
               'hiddenTreeProt': True,
               }
    
    context = base_grid(request, context)
    return render_to_response('resmap_projects.html', context)

