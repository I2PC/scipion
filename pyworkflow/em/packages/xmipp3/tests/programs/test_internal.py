# **************************************************************************
# *
# * Authors:     Josue GOmez Blanco (jgomez@cnb.csic.es)
# *
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
import os
import subprocess

from base import BaseTest
import pyworkflow.em.packages.xmipp3 as xmipp3


class TestXmipp(BaseTest):
    """ Base class to launch both Alignment and Reconstruction tests."""

    def discoverXmippTest(self):
        """ Discovers all basic xmipp tests and prints or executes them
            depending on newItemCallback argument
        """
    
        # Build base directory.
        path = "software/em/xmipp/applications/tests/"
        
        cmdList= []
        # Get all folders in current path.
        for folder in os.listdir(path):
        
            # Check if current folder starts with "test_"
            if folder.startswith("test_"):
                # Build command to execute/print.
                command = 'xmipp_%s' % folder
                cmdList.append(command)
        return cmdList
    
    def testXmippCode(self):
        for cmd in self.discoverXmippTest():
            runXmippProgram(cmd)


def runXmippProgram(cmd):
    print ">>>", cmd
    p = subprocess.Popen(cmd, shell=True, env=xmipp3.getEnviron())
    return p.wait()
