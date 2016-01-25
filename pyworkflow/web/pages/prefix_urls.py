from django.conf.urls import include, url, patterns
from pyworkflow.web.pages import settings as django_settings
from django.views.generic.base import RedirectView


urlpatterns = [
    # ... snip ...
    url('^' + django_settings.ABSOLUTE_URL, include('pages.urls')),
    # url('^.*', RedirectView.as_view(url=django_settings.ABSOLUTE_URL + 'home/')), This creates an infinite redirections
]

patterns(urlpatterns)


