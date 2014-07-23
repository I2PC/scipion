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

import json
from views_base import getResourceCss, getResourceIcon, getResourceJs 
from views_base import loadProject, base_flex, render_to_response 
from views_tree import loadObjTree
from pyworkflow.manager import Manager
from django.http import HttpResponse


def data_content(request):        
    projectName = request.GET.get('projectName', None)
    
    manager = Manager()
    request.session['projectPath'] = manager.getProjectPath(projectName)
    project = loadProject(projectName)
    
    # load the object tree 
    root = loadObjTree(project)
    
    context = {'projectName': projectName,
               'editTool': getResourceIcon('edit_toolbar'),
               'graph_utils': getResourceJs('graph_utils'),
               'project_content_utils': getResourceJs('project_content_utils'),
               'jquery_cookie': getResourceJs('jquery_cookie'),
               'jquery_treeview': getResourceJs('jquery_treeview'),
               'project_content_css':getResourceCss('project_content'),
               'sections': root.childs,
               'view':'data'
               }
    
    context = base_flex(request, context)
    
    return render_to_response('data_content/data_content.html', context)


def object_info(request):
    if request.is_ajax():
        projectName = request.session['projectName']
        objId = request.GET.get('objectId', None)
        project = loadProject(projectName)
        obj = project.mapper.selectById(objId)
        
        ioDict = {'info': str(obj),
                  'created': obj.getObjCreation(),
                  'comment': obj.getObjComment()
                  }
        
        jsonStr = json.dumps(ioDict, ensure_ascii=False)
        
    return HttpResponse(jsonStr, mimetype='application/javascript')

