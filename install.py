#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco     (jgomez@cnb.csic.es)
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
This script will generate the pw.bashrc and pw.cshrc file to include
"""
import sys
from os.path import abspath, dirname, join
FULLPATH = dirname(abspath(__file__))
BASHRC = '.bashrc'
CSHRC = '.cshrc'

def config(shell, varList, shellVar):  
    shellFile = 'pw' + shell
    print " - Creating file: ", shellFile
    f = open(shell, 'w+')
    for name, value in varList:
    	line = shellVar(name, value)
	print >> f, line
    f.close()
    print " - Include: \n      source %(shellFile)s in your %(shell)s file" % locals()


def bashVar(name, value):
    return "export %(name)s=%(value)s" % locals()

def tcshVar(name, value):
    return "setenv %(name)s %(value)s" % locals()    
    
if __name__ == '__main__':
    print "Installing Scipion in : ", FULLPATH
    print " - Creating file: ", BASHRC
    f = open(join(FULLPATH, BASHRC), 'w+')
    
    template = """
export SCIPION_HOME=$HOME/Scipion
export PW_HOME=%(FULLPATH)s

export PYTHONPATH=$PW_HOME:$PYTHONPATH
export PATH=$PW_HOME/pyworkflow/apps:$PATH

# For XMIPP
export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$PYTHONPATH
"""
    f.write(template % locals())
    f.close()
    
    print " - Include: \n      source %s \n   in your .bashrc file" % BASHRC
    VARS = [('SCIPION_HOME', '$HOME/Scipion'),
            ('PW_HOME', FULLPATH),
	    ('PYTHONPATH', '$PW_HOME:$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$PYTHONPATH'),
            ('PATH', '$PW_HOME/pyworkflow/apps:$PATH')]
	    
    config(BASHRC, VARS, bashVar)    
    config(CSHRC, VARS, tcshVar)