import os
from django.conf.urls import url
import pyworkflow as pw
from pyworkflow.web.webtools.myparticlevalidation.views import MYPVAL_SERVICE, MYPVAL_FORM_URL

MEDIA_MYPARTICLE_VALIDATION = os.path.join(pw.HOME, 'web', 'webtools', 'myparticlevalidation', 'resources')

urls = [
    (r'^resources_mypval/(?P<path>.*)$',
     'django.views.static.serve',
     {'document_root': MEDIA_MYPARTICLE_VALIDATION}
     ),
    url('^' + MYPVAL_SERVICE + '/', 'app.views_webtools.particlevalidation_projects'),
    url(r'^create_pval_project/$', 'app.views_webtools.create_particlevalidation_project'),
    url('^' + MYPVAL_FORM_URL + '/$', 'app.views_webtools.particlevalidation_form'),
    url(r'^p_content/$', 'app.views_webtools.particlevalidation_content')

]
