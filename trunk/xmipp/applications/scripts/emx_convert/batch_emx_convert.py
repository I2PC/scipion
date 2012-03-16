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
import sys
from protlib_xmipp import XmippScript
from protlib_emx import ParticlePickingConverter
from protlib_emx import ParticleAlignmentConverter
#from protlib_emx import ParticleClassConverter
from protlib_emx import CtfMicrographConverter

class ScriptEmxConverter(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Conversion utilities for Electron Microscopy exchange.')
        ## params
        self.addParamsLine(' -i <filename>                   : Input file')
        self.addParamsLine('   alias --input_filename;')
        self.addParamsLine('[ -o <filename="/dev/stdout"> ]    : Output file')
        self.addParamsLine('   alias --output_filename;')
        self.addParamsLine('[ -t <mode=coordinates>]       : Conversion type')
        self.addParamsLine('     where <mode> coordinates alignment class ctfMicrograph ctfParticle')
        self.addParamsLine('   alias --conversion_type;')

        ## examples
        self.addExampleLine('Convert coordinates file from EMX to Xmipp', False)
        self.addExampleLine('xmipp_emx_convert -i particlePicking.emx -o particlePicking.xmd -t coordinates')
        
    def run(self):   
        inputFn = self.getParam('-i')
	##next line will make life easier for java programmers
	inputFn = inputFn.replace('//','/')
        outputFn = self.getParam('-o')        
        convType = self.getParam('-t')
        if convType == 'coordinates':
            ParticlePickingConverter(inputFn, outputFn).run()  
#        elif convType == 'alignment':
#            ParticleAlignmentConverter(inputFn, outputFn).run()  
#        elif convType == 'class':
#            ParticleClassConverter(inputFn, outputFn).run()  
        elif convType == 'ctfMicrograph':
            CtfMicrographConverter(inputFn, outputFn).run()  
#        elif convType == 'ctfParticle':
#            CtfParticleConverter(inputFn, outputFn).run()  
        else:
            print >> sys.stderr, "ERROR: Wrong mode: ", type
            exit(0)

if __name__ == '__main__':
    ScriptEmxConverter().tryRun()

