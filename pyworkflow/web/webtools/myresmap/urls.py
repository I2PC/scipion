import os
from django.conf.urls import url
import pyworkflow as pw

MEDIA_MYRESMAP = os.path.join(pw.HOME, 'web', 'webtools', 'myresmap', 'resources')

urls = [
    (r'^resources_myresmap/(?P<path>.*)$', 
        'django.views.static.serve', 
        {'document_root': MEDIA_MYRESMAP}
    ),
    
    url(r'^myresmap/', 'app.views_webtools.resmap_projects'),
    url(r'^create_resmap_project/$', 'app.views_webtools.create_resmap_project'),
#     url(r'^get_testdata/$', 'app.views_webtools.get_testdata'),
    url(r'^r_form/$', 'app.views_webtools.resmap_form'),
    url(r'^r_content/$', 'app.views_webtools.resmap_content')
    
]