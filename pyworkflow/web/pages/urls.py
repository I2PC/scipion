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
    url(r'^protocol_io/$', 'app.views.protocol_io'),
    url(r'^protocol_summary/$', 'app.views.protocol_summary'),
    url(r'^project_graph/$', 'app.views.project_graph'),
    url(r'^form/$', 'app.views.form'),
    #url(r'^copy_run/$', 'app.views.copy_run'),
    url(r'^save_protocol/$', 'app.views.save_protocol'),
    url(r'^protocol/$', 'app.views.protocol'),
    url(r'^browse_objects/$', 'app.views.browse_objects'),
    url(r'^view_hosts', 'app.views.viewHosts', name='view_hosts'),    
    
#     url(r'^getHost/$', 'app.views.getHost'),
    url(r'^update_host', 'app.views.updateHost'),
    url(r'^delete host/$', 'app.views.deleteHost'),
    url(r'^host_form/$', 'app.views.hostForm'),
#    url(r'^table/', 'app.views.table'),
    url(r'^visualize_object/$', 'app.views.visualizeObject'),
    url(r'^visualize_volume/$', 'app.views.visualizeVolume'),

    #SHOWJ
    url(r'^showj/', 'app.views_showj.showj'), #Load web
    
    #url(r'^showj/$', 'app.views_showj.showj'), #Load web
    
    url(r'^get_image/', 'app.views_showj.get_image'), # Load images dynamically
    url(r'^save_showj_table/', 'app.views_showj.save_showj_table'), # Save table to session variable dynamically
    url(r'^showVolVisualization/', 'app.views.showVolVisualization'),
)
