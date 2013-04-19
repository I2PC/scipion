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
This module is responsible for launching protocol executions.
"""
import sys
from os.path import abspath, dirname
FULLPATH = abspath(__file__)
sys.path.append(dirname(dirname(dirname(FULLPATH))))
from pyworkflow.manager import Manager



if __name__ == '__main__':
    if len(sys.argv) > 2:
        projName = sys.argv[1]
        protId = int(sys.argv[2])
        
        print "="*100
        print "projName: ", projName
        print "protId: ", protId
        
        
        manager = Manager()
        project = manager.createProject(projName) # Now it will be loaded if exists
        protocol = project.mapper.selectById(protId)
        if protocol is None:
            print "Not protocol found"
        project.launchProtocol(protocol)
        #protocol.run()
        #protocol.printAll()
    else:
        print "usage: %s projectName protocolID" % sys.argv[0]
