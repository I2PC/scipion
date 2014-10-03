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
from django.shortcuts import render_to_response, HttpResponse
from django.template import RequestContext

from models import Document
from forms import DocumentForm

def upload(request):
    return renderUpload(request, DocumentForm())
        

def renderUpload(request, form):
    # Load documents for the list page

    documents = Document.objects.all()

    context = {'documents': documents, 
               'form': form}

    # Render list page with the documents and the form
    return render_to_response('upload.html', context, 
                              context_instance=RequestContext(request))

def doUpload(request):
    # Save the files
    form = DocumentForm(request.POST, request.FILES)
    
    if form.is_valid():
        newdoc = Document(docfile = request.FILES['docfile'])
        newdoc.save()
        
    return renderUpload(request, form)


def getPath(request):
    # action = request.GET.get('action')
    # time = request.GET.get('time')
    path = request.GET.get('path')
    ioDict = []
    
    for f in os.listdir(path): 
        file_path = os.path.join(path, f)
        folder = os.path.isdir(file_path)
        
        ob = {'name': f,
              'isFolder': folder,
              'isError': False}
        
        ioDict.append(ob)
            
    jsonStr = json.dumps(ioDict, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')


