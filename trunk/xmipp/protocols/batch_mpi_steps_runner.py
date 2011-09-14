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

import sys
from protlib_base import XmippProject, getProtocolFromModule
from protlib_sql import runStepGaps
from protlib_utils import printLog

if __name__ == '__main__':
    script  = sys.argv[1]
    project = XmippProject()
    project.load()
    protocol = getProtocolFromModule(script, project)
    protocol.runSetup(isMainLoop=False)
    db = protocol.Db
    #Create database and other output files
    try:
        runStepGaps(db)
    except Exception, e:
        printLog("Stopping MPI process because of error %s"%e, db.Log, out=True, err=True, isError=True)
        db.updateRunState(db.RUN_FAILED)
        exit(1)
    
    
