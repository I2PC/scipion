#!/usr/bin/env xmipp_python
"""/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *              J. M. de la Rosa Trevin
 *
 * Universidad Autonoma de Madrid
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
"""
import os
from protlib_xmipp import XmippScript
from protlib_emx import ParticlePickingConverter

class ScriptEmxConverter(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Plot some values from metadata.')
        ## params
        self.addParamsLine(' -i <filename=/dev/stdin>          : Input file')
        self.addParamsLine('   alias --input_filename;')
        self.addParamsLine(' -o <filename=/dev/stdout>         : Output file')
        self.addParamsLine('   alias --output_filename;')
        self.addParamsLine(' -t <mode=coordinates>             : Conversion type')
        self.addParamsLine('     where <mode> coordinates alignment class ctf')
        self.addParamsLine('   alias --conversion_type;')

        ## examples
        self.addExampleLine('Simple plot of label "sigmaNoise" from metadata', False)
        self.addExampleLine('xmipp_metadata_plot -i results.xmd -y sigmaNoise')
        self.addExampleLine('Additionally take values for X from other label and set title', False)
        
    def run(self):   
        inputFn = self.getParam('-i')
        outputFn = self.getParam('-o')        
        type = self.getParam('-t')
        if type == 'coordinates':
            ParticlePickingConverter(inputFn, outputFn).run()  


if __name__ == '__main__':
    ScriptEmxConverter().tryRun()

