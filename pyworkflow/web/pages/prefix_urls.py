from django.conf.urls import include, url, patterns
from pyworkflow.web.pages import settings as django_settings

urlpatterns = [
    # ... snip ...
    url('^' + django_settings.ABSOLUTE_URL, include('pages.urls')),
]

patterns(urlpatterns)


