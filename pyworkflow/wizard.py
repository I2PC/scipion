# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
"""
This module is mainly for the Viewer class, which 
serve as base for implementing visualization tools(Viewer sub-classes).
"""

DESKTOP_TKINTER = 'tkinter'
WEB_DJANGO = 'django'


class Wizard(object):
    """ This is a special case of GUI to help the user
    selecting important parameters.
    The _targets will serve to define to which Definition and 
    parameters the Wizard is defined, it will be a list of tuples such as:
    _targets = [(DefImportMicrographs, ['voltage', sphericalAberration']),
                (DefCTFMicrographs, ['lowRes', 'highRes'])]
    The _environmets will serve to define when this wizard can be used.
    For example>
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    """
    _targets = []
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def show(self, form, *params):
        """ This will show up the wizard to select parameters.
        Params:
            form: the protocol form, given access to to all parameters.
                Some times the same wizard will modifify several elements
                in the form.
            *params: a list of params to modify, sometimes the wizard can 
                be generic and can be used for different parameters in the
                same form.
        """
        pass
    
    def getView(self):
        """ This method should return the string value of the view in web
        that will respond to this wizard. This method only should be implemented
        in those wizards that have WEB_DJANGO environment defined. 
        """
        return None
    
