#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco     (jgomez@cnb.csic.es)
# *              I. Foche Perez         (ifoche@cnb.csic.es)
# *              J. Burguet Castell     (jburguet@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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

"""
Execute scons in order to compile Xmipp code.
"""

import sys
import os

import subprocess


SCIPION_HOME = os.environ['SCIPION_HOME']
LOGFILE = os.path.join(SCIPION_HOME, 'software', 'log', 'scons.log')



def build(args):
    """ Run scons with our python (and with the passed args) """
    # Call SCons and show the output on the screen and logfile.
    proc = subprocess.Popen(['python', 'software/bin/scons'] + args,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            cwd=SCIPION_HOME)
    with open(LOGFILE, 'a') as logFile:
        while True:
            ret = proc.poll()
            if ret is not None:
                for line in proc.stdout:  # output the remaining lines
                    sys.stdout.write(line)
                    logFile.write(line)
                    logFile.flush()
                return ret
            line = proc.stdout.readline()
            sys.stdout.write(line)
            logFile.write(line)
            logFile.flush()



if __name__ == '__main__':
    sys.exit(build(sys.argv[1:]))
