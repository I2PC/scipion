from django.conf.urls import patterns, include, url
# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()
from django.conf import settings

urlpatterns = patterns('',
    (r'^resources/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT}),
    (r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_ROOT}),
                
    # Examples:
    # url(r'^pages/', include('pages.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^pages/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
    
    url(r'^projects/', 'app.views.projects'),
    url(r'^create_project/$', 'app.views.create_project'),
    url(r'^delete_project/$', 'app.views.delete_project'),
    url(r'^project_content/$', 'app.views.project_content'),
    url(r'^delete_protocol/$', 'app.views.delete_protocol'),
    url(r'^form/$', 'app.views.form'),
    url(r'^copy_run/$', 'app.views.copy_run'),
    url(r'^protocol/$', 'app.views.protocol'),
    url(r'^browse_objects/$', 'app.views.browse_objects'),
    url(r'^openHostsConfig/$', 'app.views.openHostsConfig'),
    url(r'^getHost/$', 'app.views.getHost'),
    url(r'^updateHostsConfig', 'app.views.updateHostsConfig'),
    url(r'^showj/', 'app.views.showj'),
    url(r'^table/', 'app.views.table'),
    url(r'^get_image/', 'app.views.get_image'), # Load images dynamically
    url(r'^get_table/', 'app.views.get_table'), # Load table dynamically
)
