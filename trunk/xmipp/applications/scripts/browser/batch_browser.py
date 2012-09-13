#!/usr/bin/env xmipp_python
'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

from protlib_xmipp import XmippScript
from protlib_gui_ext import showBrowseDialog

class ScriptBrowser(XmippScript):
    def __init__(self):
        XmippScript.__init__(self, True)
        
    def defineParams(self):
        self.addUsageLine('Browse some directory')
        ## params
        self.addParamsLine('[ -i <directory=".">]          : Initial directory to start browsing')
        self.addParamsLine('   alias -d;')
        self.addParamsLine('[ -f <pattern="">]             : Filter results files')
        self.addParamsLine('   alias --filter;')
        
        self.addExampleLine('Just open the browser in current directory', False)
        self.addExampleLine('xmipp_browser')
        self.addExampleLine('Only show .xmd and .stk files:', False)
        self.addExampleLine('xmipp_browser -f ".stk .xmd"')        
            
    def run(self):
        path = self.getParam('-i')
        filter = self.getParam('--filter')
        showBrowseDialog(path, main=True, seltype="none", filter=filter)

if __name__ == '__main__':
    ScriptBrowser().tryRun()
