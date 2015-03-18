import os
from django.conf.urls import url
import pyworkflow as pw

MEDIA_MYRESMAP = os.path.join(pw.HOME, 'web', 'webtools', 'myresmap', 'resources')

urls = [
    (r'^resources_resmap/(?P<path>.*)$', 
        'django.views.static.serve', 
        {'document_root': MEDIA_MYRESMAP}
    ),
    
    url(r'^myresmap/', 'app.views_webtools.resmap_projects'),
#     url(r'^intro/', 'app.views_webtools.service_projects'),
#     url(r'^check_project_id/$', 'app.views_webtools.check_project_id'),
#     url(r'^create_service_project/$', 'app.views_webtools.create_service_project'),
#     url(r'^get_testdata/$', 'app.views_webtools.get_testdata'),
#     url(r'^my_form/$', 'app.views_webtools.myfirstmap_form'),
#     url(r'^content/$', 'app.views_webtools.service_content')
    
]