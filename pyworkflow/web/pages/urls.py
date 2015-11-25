import os, sys
import pyworkflow as pw
from django.conf.urls import patterns, include, url
# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()
from django.conf import settings
from pyworkflow.web.pages.settings import WS_ROOT, serviceFolders
from django.views.generic import TemplateView

#===============================================================================
# URL ASSOCIATION
#===============================================================================

mainUrls = ['',
    # To serve different static files
    (r'^resources/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT}),
    (r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.STATIC_ROOT}),
    
    url(r'^admin/', include(admin.site.urls)),
    # url(r'^pages/doc/', include('django.contrib.admindocs.urls')),
    
    # If no path given, load the projects view
#     url(r'^$', 'app.views_project.projects'),
    url(r'^$', 'app.views_home.home'),
    url(r'^services/', 'app.views_project.webservice_projects'),
    
    #PROJECT (CONTENT, RUNTABLE AND GRAPH)
    url(r'^projects/', 'app.views_project.projects'),
    url(r'^create_project/$', 'app.views_project.create_project'),
    
    url(r'^check_project_id/$', 'app.views_project.check_project_id'),
    url(r'^delete_project/$', 'app.views_project.delete_project'),
    url(r'^project_content/$', 'app.views_project.project_content'),
    url(r'^get_protocols/$', 'app.views_project.get_protocols'),
    url(r'^search_protocol/$', 'app.views_project.search_protocol'),
    url(r'^tree_prot_view/$', 'app.views_project.tree_prot_view'),
    url(r'^run_table_graph/$', 'app.views_project.run_table_graph'),
    url(r'^protocol_info/$', 'app.views_project.protocol_info'),
    url(r'^update_graph_view/$', 'app.views_project.update_graph_view'),
    url(r'^update_prot_tree/$', 'app.views_project.update_prot_tree'),
    url(r'^save_selection/$', 'app.views_project.save_selection'),
    
    #GRAPH (RUN & DATA)
    url(r'^project_graph/$', 'app.views_graph.project_graph'),
    url(r'^object_graph/$', 'app.views_graph.object_graph'),
    url(r'^elements_graph/', 'app.views_graph.elements_graph'),
    
    #DATA (CONTENT)
    url(r'^data_content/$', 'app.views_data.data_content'),
    url(r'^object_info/$', 'app.views_data.object_info'),
    url(r'^object_tree/$', 'app.views_data.object_tree'),
    
    #UTILS
    url(r'^error/', 'app.views_util.error'), # Launch error page
    url(r'^render_column/', 'app.views_util.render_column'), # Load images dynamically
    url(r'^get_image_plot/', 'app.views_util.get_image_plot'), # Load plots images dynamically
    url(r'^get_image_path/', 'app.views_util.get_image_path'), # Load images dynamically
    url(r'^get_image/', 'app.views_util.get_image'), # Load images dynamically
    url(r'^get_slice/', 'app.views_util.get_slice'), # Load slices dynamically
    url(r'^browse_objects/$', 'app.views_util.browse_objects'), # Browse objects from the database
    url(r'^browse_relations/$', 'app.views_util.browse_relations'), # Browse relation objects from the database
    url(r'^browse_protocol_class/$', 'app.views_util.browse_protocol_class'), # Browse objects from the database
    url(r'^get_attributes/$', 'app.views_util.get_attributes'), # Get Label and Comment for an Object
    url(r'^set_attributes/$', 'app.views_util.set_attributes'), # Set Label and Comment for an Object
    url(r'^download_output/$', 'app.views_util.download_output'), # Download output files

    #PROTOCOL (INCLUDE FORM)
    url(r'^save_protocol/$', 'app.views_protocol.save_protocol'),
    url(r'^protocol/$', 'app.views_protocol.protocol'),
    url(r'^stop_protocol/$', 'app.views_protocol.stop_protocol'),
    url(r'^delete_protocol/$', 'app.views_protocol.delete_protocol'),
    url(r'^copy_protocol/$', 'app.views_protocol.copy_protocol'),
    url(r'^form/$', 'app.views_protocol.form'),

    #WIZARDS
    url(r'^wizard/$', 'app.em_wizard.wizard'),
    url(r'^get_image_mask/$', 'app.wizards.tools.get_image_mask'),
    url(r'^get_image_psd/$', 'app.wizards.tools.get_image_psd'),
    url(r'^get_image_bandpass/$', 'app.wizards.tools.get_image_bandpass'),
    url(r'^get_image_gaussian/$', 'app.wizards.tools.get_image_gaussian'),
    url(r'^get_image_filter_spider/$', 'app.wizards.spider_wizard.get_image_filter_spider'),
    url(r'^run_custom_mask_spider/$', 'app.wizards.spider_wizard.run_custom_mask_spider'),
    url(r'^get_resmap_plot/$', 'app.wizards.resmap_wizard.get_resmap_plot'),
    
    #VIEWERS
    url(r'^launch_viewer/$', 'app.em_viewer.launch_viewer'),
    url(r'^viewer_element/$', 'app.em_viewer.viewer_element'),
    url(r'^file_viewer/$', 'app.views_util.file_viewer'),
    url(r'^get_log/$', 'app.views_util.get_log'),
    
    #SHOWJ 
    url(r'^showj/$', 'app.views_showj.showj'), #Load showj web
    url(r'^update_session_table/$', 'app.views_showj.updateSessionTable'),
    url(r'^jsmol/$', 'app.views_showj.jsmol'),
    url(r'^get_chimera_html/$', 'app.views_showj.get_chimera_html'),
    
    #BROWSER & UPLOAD FILES
    url(r'^upload/', 'app.views_management.upload', name='upload'),
    url(r'^doUpload/', 'app.views_management.doUpload'),
    url(r'^upload_service/', 'app.views_management.upload_service'),
    url(r'^getPath/', 'app.views_management.getPath'),
    url(r'^getExtIcon/$', 'app.views_management.getExtIcon'),
    url(r'^get_file/$', 'app.views_util.get_file'),
    url(r'^get_image_dim/$', 'app.views_util.get_image_dim'),
    
    #TESTING
#    url(r'^testingSSH/', 'app.views_showj.testingSSH'), #Load web

    
    url(r'^home/', 'app.views_home.home'),
    url(r'^download_form/', 'app.views_home.download_form'),
    url(r'^doDownload/', 'app.views_home.doDownload'),

]

# Load URLS for webtools
from pyworkflow.utils.reflection import getModules
toolModules = getModules(WS_ROOT)

for tm in toolModules.values():
    mainUrls += tm.urls
    
 

# handler404 = "app.views_util.error"
# handler500 = "app.views_util.error"

urlpatterns = patterns(*mainUrls)

