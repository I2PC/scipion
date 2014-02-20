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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.em import *
from views_util import * 

from pyworkflow.web.app.properties import MessageWeb


def base(request, context):
    
    # Messages properties class
    messages = MessageWeb()
    
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
                    'utils': getResourceJs('utils'),
                    #OTHER
                    'msg': messages
                    }
    
    context.update(context_base)
    return context


def base_form(request, context):
    
    context_base = {
                    #CSS
                    'form_css': getResourceCss('form'),
                    'jquery_ui_css': getResourceCss('jquery_ui'),
                    #JS
                    'jquery_ui': getResourceJs('jquery_ui'),
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
                    'jsplumb': getResourceJs('jsplumb'),
                    'jquery_sizes': getResourceJs('jquery_sizes'),
                    'jlayout_border': getResourceJs('jlayout_border'),
                    'jquery_layout': getResourceJs('jquery_layout'),
                    #OTHER
                    'contentConfig': 'divided'
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
                    'jquerydataTables_colreorder_resize': getResourceJs('jquery_colreorder_resize'),
                    'jeditable': getResourceJs('jquery_editable'),
                    'jquery_waypoints':getResourceJs('jquery_waypoints'),
                    'jquery_hover_intent':getResourceJs('jquery_hover_intent'),
                    'showj_js':getResourceJs('showj_utils'),
                    'jquery_ui':getResourceJs('jquery_ui'),
                    'transpose_lib':getResourceJs('transpose'),
                    }
    
    context = base_grid(request, context)
    context.update(context_base)
    return context

