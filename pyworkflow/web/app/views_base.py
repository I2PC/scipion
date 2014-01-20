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

from pyworkflow.web.app.messages_properties import MessageWeb

def base(request, context):
    
    # Messages properties class
    messages = MessageWeb()
    
    context_base = {
                'favicon': getResourceIcon('favicon'),
                'general_css': getResourceCss('general'),
                'messi_css': getResourceCss('messi'),
                'font_awesome': getResourceCss('font_awesome'),
                'jquery': getResourceJs('jquery'),
                'messi_js': getResourceJs('messi'),
                'utils': getResourceJs('utils'),
                'messages': messages
               }
    
    context.update(context_base)
    return context

def base_form(request, context):
    context_base = {
                    'form_css': getResourceCss('form'),
                    'jquery_ui': getResourceJs('jquery_ui'),
                    'jquery_ui_css': getResourceCss('jquery_ui')
                    }

    context = base(request, context)
    context.update(context_base)
    return context

def base_grid(request, context):
    
    context_base = {
                    'general_grid': getResourceCss('general_grid'),
                    'jquery_sizes': getResourceJs('jquery_sizes'),
                    'jlayout_border': getResourceJs('jlayout_border'),
                    'jquery_layout': getResourceJs('jquery_layout')
                    }
    
    context = base(request, context)
    context.update(context_base)
    return context

def base_flex(request, context):
    
    context_base = {
                    'general_flex': getResourceCss('general_flex'),
                    'jquery_ui': getResourceJs('jquery_ui'),
                    'jquery_ui_css': getResourceCss('jquery_ui'),
                    'jsplumb': getResourceJs('jsplumb'),
                    'jquery_sizes': getResourceJs('jquery_sizes'),
                    'jlayout_border': getResourceJs('jlayout_border'),
                    'jquery_layout': getResourceJs('jquery_layout'),
                    'contentConfig': 'divided'
                    }
    
    context = base(request, context)
    context.update(context_base)
    return context

