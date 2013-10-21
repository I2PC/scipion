# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This module contains templates files from Spider scripts that will 
generated a file that can be executed from Spider, replacing specific
variables.
"""

from os.path import join, dirname, abspath, exists

TEMPLATE_DIR = abspath(dirname(__file__))

def createSpiderScript(templateFile, scriptFile, paramsDict):
    """ This function will read a template, substitute the params with values
    and write the resulting script. """
    f = open(templateFile, 'r')
    
    for line in f:
        line = line.strip()
        if not line.startswith(';'):
            line = line % paramsDict
        print line
        
        
def getTemplate(templateName):
    """ Return the path to the template file given its name. """
    templateFile = join(TEMPLATE_DIR, templateName)
    if not exists(templateFile):
        raise Exception("getTemplate: template '%s' was not found in templates directory" % templateName)
    
    return templateFile
    