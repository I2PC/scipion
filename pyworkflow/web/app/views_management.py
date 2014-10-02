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
    # Not Used
    # action = request.GET.get('action')
    # time = request.GET.get('time')
    
    path_browse = request.GET.get('path')
    
    ioDict = []
    for fi in os.listdir(path_browse): 
        file_path = os.path.join(path_browse , fi)
        
        folder = False
        if os.path.isdir(file_path):
            folder = True
            
        ob = {'name': fi,'isFolder': folder,'isError': False}
        ioDict.append(ob)
        
            
    jsonStr = json.dumps(ioDict, ensure_ascii=False)
    return HttpResponse(jsonStr, mimetype='application/javascript')


