#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     R. Marabini(roberto@cnb.csic.es)
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
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************


import os
import sqlite3 as lite
from matplotlib import pyplot
from time import sleep
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Plot cpu and mem usage. Usage: scipion run ./plotter.py --cpu --mem')
    parser.add_argument('--cpu', '-c',
          action='store_true',
          help='cpu flag' )
    parser.add_argument('--mem', '-m',
          action='store_true',
          help='mem flag' )
    args = parser.parse_args()
    doCPU = args.cpu
    doMEM = args.mem
    if doCPU == False and doMEM == False:
         doMEM=True

    # Get arguments. This will be implement later
    baseFn = os.path.join(os.environ.get('SCIPION_USER_DATA'),'tmp','log.sqlite')
    tableName = 'log'
    conn = lite.connect(baseFn, isolation_level=None)
    cur = conn.cursor()
    #conn.text_factory = str

    #I guess this can be done in a single call
    cur.execute("select julianday(timestamp)  from %s where id=1"%tableName )
    initTime = cur.fetchone()[0]
    cur.execute("select (julianday(timestamp) - %f)*24  from %s"%(initTime,tableName) )
    id=[r[0] for r in cur.fetchall()]
    cur.execute("select cpu  from %s"%tableName )
    cpu=[r[0] for r in cur.fetchall()]
    cur.execute("select mem from %s"%tableName )
    mem=[r[0] for r in cur.fetchall()]

    if doCPU:
       pyplot.plot(id,cpu,'-b', label='cpu')
    if doMEM:
       pyplot.plot(id,mem,'-r', label='mem')

    x1,x2,y1,y2 = pyplot.axis()
    pyplot.axis((0.,x2,y1,y2))
    pyplot.title("System Monitor")
    pyplot.xlabel("time (hours)")
    pyplot.ylabel("percentage")
    pyplot.legend()
    pyplot.show()



if __name__ == "__main__":
    main()
