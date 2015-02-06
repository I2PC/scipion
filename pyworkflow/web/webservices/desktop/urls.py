import os
from django.conf.urls import url

urls = [
#     url(r'^desktop/', 'app.views_desktop.desktop'),
    url(r'^download_form/', 'app.views_webservice.download_form'),
    url(r'^doDownload/', 'app.views_webservice.doDownload'),
    
]