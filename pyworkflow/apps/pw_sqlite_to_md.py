#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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
Object browser
"""

import sys
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.em.data import SetOfParticles
from pyworkflow.mapper.sqlite import SqliteFlatDb
from pyworkflow.em import *






USAGE = "usage: pw_sqlite_to_md.py pathToSqlite pathToMd"

if __name__ == '__main__':


    if len(sys.argv) == 4:
        pathToSqlite = sys.argv[1]
        preffix = sys.argv[2]
        pathToMd = sys.argv[3]

        db = SqliteFlatDb(dbName=pathToSqlite, tablePrefix=preffix)
        setClassName = db.getProperty('self') # get the set class name
        setObj = emObjectsDict[setClassName](filename=pathToSqlite, prefix=preffix)
        print setObj.getClass()
        writeSetOfParticles(set, pathToMd)


    else:
        print USAGE

