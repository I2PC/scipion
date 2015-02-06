import os
from django.conf.urls import url

MEDIA_MYFIRSTMAP = os.path.join(pw.HOME, 'web', 'webservices', 'myfirstmap', 'resources')

urls = [
    (r'^resources_myfirstmap/(?P<path>.*)$', 
        'django.views.static.serve', 
        {'document_root': MEDIA_MYFIRSTMAP}
    ),
    
    url(r'^service_projects/', 'app.views_webservice.service_projects'),
    url(r'^check_project_id/$', 'app.views_project.check_project_id'),
    url(r'^create_service_project/$', 'app.views_project.create_service_project'),
    url(r'^get_testdata/$', 'app.views_project.get_testdata'),
    url(r'^service_content/$', 'app.views_project.service_content')
    
]