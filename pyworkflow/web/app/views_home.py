#!/usr/bin/python
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
from django.shortcuts import render_to_response
from pyworkflow.web.app.views_util import getResourceCss, getResourceIcon, getResourceJs, parseText
from pyworkflow.web.app.views_base import base_grid
from django.http import HttpResponse, HttpResponseNotFound , HttpResponseBadRequest
from django.core.context_processors import csrf
from pyworkflow.web.pages import settings as django_settings
from pyworkflow.mapper.sqlite import SqliteFlatMapper
from django.core.servers.basehttp import FileWrapper
import mimetypes

import pyworkflow as pw
from pyworkflow.config import DownloadRecord
DB_PATH_DOWNLOAD = os.path.join(pw.HOME, 'web', 'home', "download_statistics.db")

def home(request):
    context = {}
    context = base_grid(request, context)
    return render_to_response('home/index.html', context)

def download_form(request):
    
    desktop_utils = django_settings.STATIC_URL + "js/download_utils.js"
    context = {"download_utils": desktop_utils }
    context = base_grid(request, context)
    context.update(csrf(request))
    return render_to_response('home/download_form.html', context)

def doDownload(request):
    
    fullName = request.POST.get('fullName')
    organization = request.POST.get('organization')
    email = request.POST.get('email')
    mailoption = request.POST.get('mailoption')
    country = request.POST.get('country')
    version = request.POST.get('version')
    platform = request.POST.get('platform')
    
    errors = ""
    
    if not len(fullName) > 0:
       errors += "Please fill in the fullName field.\n"
    if not len(organization) > 0:
        errors += "Please fill in the Organization field.\n"
    if not len(email) > 0:
        errors += "Please fill in the Email field.\n"
#     if not len(mailoption) > 0:
#         errors += "Please choose one into the Country field.\n"
    if not len(version) > 0:
        errors += "Please fill in the Scipion Version field.\n"
    if not len(platform) > 0:
        errors += "Please fill in the Platform field.\n"

    if len(errors) == 0:
        dbName = os.path.join(os.environ['SCIPION_HOME'], 'downloads.sqlite')
        #dbName = '/tmp/downloads.sqlite'
        
        mapper = SqliteFlatMapper(dbName, globals())
        mapper.enableAppend()
        download = DownloadRecord(fullName = fullName,
                organization = organization,
                email = email,
                subscription = mailoption,
                country = country,
                version = version,
                platform = platform)
        

        mapper.store(download)
        mapper.commit()
        mapper.close()
        "Return a response with the scipion download file"
        if platform == 'linuxbin':
            path = os.path.join(os.environ['SCIPION_HOME'], 'pyworkflow', 'web', 'pages', 'resources', 'install', 'scipion_all_packages_2015-06-29.tgz')
        else:
            path = os.path.join(os.environ['SCIPION_HOME'], 'pyworkflow', 'web', 'pages', 'resources', 'install', 'scipion_source_2015-06-29.tgz')
        if not os.path.exists(path):
            return HttpResponseNotFound('Path not found: %s' % path)
    
        response = HttpResponse(FileWrapper(open(path)),
                                content_type=mimetypes.guess_type(path)[0])
        response['Content-Length'] = os.path.getsize(path)
        response['Content-Disposition'] = 'attachment; filename=%s'%os.path.basename(path) 

        return response
    else:
        jsonStr = json.dumps({'errors' : parseText(errors)}, ensure_ascii=False)
    
        return HttpResponse(jsonStr, mimetype='application/javascript')   


