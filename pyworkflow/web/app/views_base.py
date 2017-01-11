# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from views_util import getResourceCss, getResourceIcon, getResourceJs
from pyworkflow.web.pages import settings as django_settings

def VARS_base(request, context):
    from pyworkflow.web.app.properties import MessageWeb, ColorWeb, IconWeb

    # Properties class
    messages = MessageWeb()
    colors = ColorWeb()
    icons = IconWeb()
    
    context_base = {
                    #ABSOLUTE PATH URL CONFIG
                    'abs_url': django_settings.ABSOLUTE_URL, 
                    'config': getResourceJs('config'),
                    'js_root': getResourceJs(),
                    'css_root': getResourceCss(),

                    #OTHER
                    'msg': messages,
                    'color': colors,
                    'icon': icons,
                    }
    
    context.update(context_base)
    return context
    

def base(request, context):
    
    context_base = {
                    #ICON
                    'favicon': getResourceIcon('favicon'),
                    'logo_scipion': getResourceIcon('logo_scipion'),
                    'logo_scipion_transparent' : getResourceIcon('logo_scipion_transparent'),
                    #CSS
                    'general_css': getResourceCss('general'),
                    'messi_css': getResourceCss('messi'),
                    'font_awesome': getResourceCss('font_awesome'),
                    #JS
                    'jquery': getResourceJs('jquery'),
                    'messi_js': getResourceJs('messi'),
                    #JS
                    'utils': getResourceJs('utils'),
                    }
    
    context = VARS_base(request, context)
    context.update(context_base)
    return context


def base_form(request, context):
    
    context_base = {
                    #Folder
                    'project_folder': request.session['projectPath'],
                    #CSS
                    'form_css': getResourceCss('form'),
                    'jquery_ui_css': getResourceCss('jquery_ui'),
                    #JS
                    'jquery_ui': getResourceJs('jquery_ui'),
                    'jquery_ui_touch': getResourceJs('jquery_ui_touch'),
                    'jquery_browser': getResourceJs('jquery_browser'),
                    "upload_utils": getResourceJs('upload_utils'),
                    }

    context = base(request, context)
    context.update(context_base)
    return context


def base_flex(request, context):
    
    context_base = {
                    #CSS
                    'general_flex': getResourceCss('general_flex'),
                    'jquery_ui_css': getResourceCss('jquery_ui'),
                    #JS
                    'jquery_ui': getResourceJs('jquery_ui'),
                    'jquery_ui_touch': getResourceJs('jquery_ui_touch'),
                    'jsplumb': getResourceJs('jsplumb'),
                    'jquery_sizes': getResourceJs('jquery_sizes'),
                    'jlayout_border': getResourceJs('jlayout_border'),
                    'jquery_layout': getResourceJs('jquery_layout'),
                    }
    
    context = base(request, context)
    context.update(context_base)
    return context


def base_grid(request, context):
    
    context_base = {
                    #CSS
                    'general_grid': getResourceCss('general_grid'),
                    #JS
                    'jquery_sizes': getResourceJs('jquery_sizes'),
                    'jlayout_border': getResourceJs('jlayout_border'),
                    'jquery_layout': getResourceJs('jquery_layout')
                    }
    
    context = base(request, context)
    context.update(context_base)
    return context


def base_showj(request, context):
    
    context_base = {
                    #CSS
                    'showj_css': getResourceCss('showj'),
                    'smoothness': getResourceCss('ui_smoothness'),
                    'demo_table_jui': getResourceCss('showj_demo_table_jui'),
                    #JS
                    'jquery_datatable': getResourceJs('jquery_datatables'),
                    'jquerydataTables_colreorder': getResourceJs('jquery_colreorder'),
                    'jeditable': getResourceJs('jquery_editable'),
                    'jquery_waypoints':getResourceJs('jquery_waypoints'),
                    'jquery_hover_intent':getResourceJs('jquery_hover_intent'),
                    'showj_js':getResourceJs('showj_utils'),
                    'jquery_ui':getResourceJs('jquery_ui'),
                    'transpose_lib':getResourceJs('transpose'),
                    'showj_menu_utils': getResourceJs('showj_menu_utils'),
                    }
    
    context = base_grid(request, context)
    context.update(context_base)
    return context



def base_wiz(request, context):
    
    context_base = {'general_style': getResourceCss('general'),
                    'wizard_style': getResourceCss('wizard'),
                    'jquery_ui_style': getResourceCss('jquery_ui'),
                    'font_awesome': getResourceCss('font_awesome'),
                    'jquery': getResourceJs('jquery'),
                    'jquery_ui': getResourceJs('jquery_ui'),
                    'jquery_ui_touch': getResourceJs('jquery_ui_touch'),
                    'wizard_utils': getResourceJs('wizard_utils'),
                    'raphael':getResourceJs('raphael'),
                    'projectName': request.session['projectName'],
                    'loading' : getResourceIcon('loading'),
                    }
    
    context = VARS_base(request, context)
    context.update(context_base)
    return context