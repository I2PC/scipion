from django.shortcuts import render_to_response
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