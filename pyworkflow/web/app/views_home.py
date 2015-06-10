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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
import json
from django.shortcuts import render_to_response
from pyworkflow.web.app.views_util import getResourceCss, getResourceIcon, getResourceJs, parseText
from pyworkflow.web.app.views_base import base_grid
from django.http import HttpResponse
from django.core.context_processors import csrf
from pyworkflow.web.pages import settings as django_settings

import pyworkflow as pw
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
        errors += "Please choose one into the Scipion Version field.\n"
    if not len(platform) > 0:
        errors += "Please choose one into the Platform field.\n"

    if len(errors) == 0:
        # Save statistics into DB
        data = {"name": fullName,
                "org" : organization,
                "mail": email,
                "subscription" : mailoption,
                "country": country,
                "version": version,
                "platform": platform }
        
        # Update database with the new data
#         createDB()
        updateDB(data)
        
    jsonStr = json.dumps({'errors' : parseText(errors)}, ensure_ascii=False)
    
    return HttpResponse(jsonStr, mimetype='application/javascript')   

def createDB():
    import sqlite3

    conn = sqlite3.connect(DB_PATH_DOWNLOAD)
    
    conn.execute('''CREATE TABLE DOWNLOAD (
        NAME           TEXT    NOT NULL,
        ORG            TEXT     NOT NULL,
        MAIL           TEXT    NOT NULL,
        SUBSCRIPTION   INTEGER    NOT NULL,
        COUNTRY    TEXT    NOT NULL,
        VERSION     TEXT    NOT NULL,
        PLATFORM    TEXT    NOT NULL
        );''')
    print "Table created successfully";
    conn.close()

def updateDB(data):
    import sqlite3

    conn = sqlite3.connect(DB_PATH_DOWNLOAD)
    print "Opened database successfully";
    
    print "data: ", data
    
    name = data['name']
    org = data['org']
    mail = data['mail']
    subscription = data['subscription']
    country = data['country']
    version = data['version']
    platform = data['platform']
    
#     table = "downloads (NAME, ORG, MAIL, SUBSCRIBE, COUNTRY, VERSION, PLATFORM)"
    table = "DOWNLOAD"
    values = "VALUES ('"+ name +"', '"+ org +"', '"+ mail +"', '"+ subscription +"', '"+ country +"', '"+ version +"', '"+ platform +"')"
    
    conn.execute("INSERT INTO " + table + " " + values);
    conn.commit()
    
    print "Records created successfully";
    conn.close()
    
