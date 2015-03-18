import os
from django.conf.urls import url
import pyworkflow as pw

MEDIA_MOVIES = os.path.join(pw.HOME, 'web', 'webtools', 'movies', 'resources')

urls = [
    (r'^resources_movies/(?P<path>.*)$', 
        'django.views.static.serve', 
        {'document_root': MEDIA_MOVIES}
    ),
    
    url(r'^mymovies/', 'app.views_webtools.service_movies'),
    url(r'^create_movies_project/$', 'app.views_webtools.create_movies_project'),
    url(r'^mov_form/$', 'app.views_webtools.movies_form'),
    url(r'^m_content/$', 'app.views_webtools.movies_content'),
    url(r'^upload_movies/', 'app.views_webtools.upload_movies'),
    
]