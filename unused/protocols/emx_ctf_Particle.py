#!/usr/bin/env xmipp_python
'''
/***************************************************************************
 * Authors:     Roberto Marabini Ruiz
 *
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
import emx_struct
"""
1) read star file in star object
2) parse star object an assign to data structure EMX/XMIPP
3) convert EMX to XMIPP (or vice versa)
4) parse XMIPP/EMX data structure to star object 
5) save star file
"""
from emx_struct import CtfMicrographStructEmx,\
                       CtfMicrographStructXmd,\
                       BlockNamesEMX,\
                       BlockNamesXMD,\
                       prefix_micrograph
from protlib_emx import EmxBase
import CifFile
import StarFile

###########################################################################
##   Class Related to CTF conversion: defocus per particle
###########################################################################
class CtfParticleConverter(EmxBase):    
        
    def __init__(self, inputFn, outputFn):
        self.inputFileName  = inputFn
        self.outputFileName = outputFn       
    
    def run(self):
        print """
********************************************************************
Not implemented yet, EMX approach does not follows XMIPP philosophy
This will not be implemented until XMIPP 3.0 has been released
*********************************************************************
"""
