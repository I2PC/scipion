from django.conf.urls import include, url, patterns
from pyworkflow.web.pages import settings as django_settings
from django.views.generic.base import RedirectView


urlpatterns = [
    # To publish in a subpath
    url('^' + django_settings.ABSOLUTE_URL, include('pages.urls')),
    url(r'^favicon\.ico$', RedirectView.as_view(url=django_settings.STATIC_URL + 'favicon.ico')),

    # Anything not starting with the absolute url goes to home:
    # url('^(?!' + django_settings.ABSOLUTE_URL + ')', RedirectView.as_view(url='/' + django_settings.ABSOLUTE_URL + 'home/'), name='outside')
]

patterns(urlpatterns)


