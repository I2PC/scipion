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

import os
import json
import shutil
from django.shortcuts import render_to_response, HttpResponse
from django.template import RequestContext
from pyworkflow.web.pages import settings as django_settings
from models import Document
from forms import DocumentForm
from views_base import base_form
from views_util import getResourceIcon, getResourceJs

def upload(request, form=None):
    # Load documents for the list page
    mode = request.GET.get('mode', None)
    path = os.path.join(request.session['projectPath'],'Uploads')
    split_path = path.split("/ScipionUserData/")
    relative_path = "ScipionUserData/" + split_path[1]

    context = {'relative_path': relative_path,
               'form': DocumentForm(),
               'logo_scipion_small': getResourceIcon('logo_scipion_small'),
               "upload_utils": getResourceJs('upload_utils'),
               "mode": mode,
               }

    context = base_form(request, context)
    
    # Render list page with the documents and the form
    return render_to_response('upload/upload.html', context, 
                              context_instance=RequestContext(request))

def executeUpload(request):
    # Save the files
    form = DocumentForm(request.POST, request.FILES)
    
    file_new = request.FILES['docfile']
    
    if form.is_valid():
        #Save temporary file
        newdoc = Document(docfile = file_new)
        newdoc.save()
    
        fn = file_new.name
        fn = fn.replace (" ", "_")
    
        #Move the file to the new folder
        src = os.path.join(django_settings.FILE_UPLOAD_TEMP_DIR, 'uploads', fn)
        file_upload = src
        path = os.path.join(request.session['projectPath'],'Uploads')
        target = os.path.join(path, fn)
        if os.path.exists(target):
            os.remove(target)
        shutil.move(src, path)
        
        #Delete the temporary file
        newdoc.delete()
        
def doUpload(request):
    form = DocumentForm(request.POST, request.FILES)
    
    try:
        executeUpload(request)
    except Exception, ex:
        print "Error: %s" % ex
        return HttpResponse("error", mimetype='application/javascript')

    return upload(request, form)

# """ File Browser Utils """
def getPath(request):
    # action = request.GET.get('action')
    # time = request.GET.get('time')
    path = request.GET.get('path')
    ioDict = []
    
    for f in os.listdir(path): 
        file_path = os.path.join(path, f)
        folder = os.path.isdir(file_path)
        ext = f.split('.').pop();
        
        ob = {'name': f,
              'isFolder': folder,
              'isError': False,
              'icon': getExtIconPath(ext)}
        
        ioDict.append(ob)
            
    jsonStr = json.dumps(ioDict, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')


def getExtIcon(request):
    ext = request.GET.get('ext')
    res = getExtIconPath(ext)
    return HttpResponse(res, mimetype='application/javascript')

def getExtIconPath(ext):

    txt = {'txt', 'log', 'out', 'err', 'stdout', 'stderr', 'emx'}
    img = {'png', 'gif', 'jpg', 'jpeg'}
    py = {'py', 'pyc'} 
    java = {'java'}
    md = {'xmd', 'star', 'pos'}
    sqlite = {'sqlite', 'db'}
    particle = {'xmp', 'tif', 'tiff', 'spi', 'mrc', 'map', 'raw', 
                'inf', 'dm3', '.em', 'pif', 'psd', 'spe', 
                'ser', 'img', 'hed'}
    vol = {'vol'}
    stk = {'stk', 'mrcs', 'st', 'pif'}
    
    if ext == 'folder':
        res = getResourceIcon('folder')
#     elif ext == 'unknown':
#         res = getResourceIcon('file_normal')
    else:
        if ext in txt:
            res = getResourceIcon('file_text')
        elif ext in img or ext in particle:
            res =  getResourceIcon('file_image')
        elif ext in py:
            res = getResourceIcon('file_python')
        elif ext in java:
            res = getResourceIcon('file_java')
        elif ext in md:
            res = getResourceIcon('file_md')
        elif ext in sqlite:
            res = getResourceIcon('file_sqlite')
        elif ext in vol:
            res = getResourceIcon('file_vol')
        elif ext in stk:
            res = getResourceIcon('file_stack')
        else:
            res = getResourceIcon('file_normal')
    
    return res

