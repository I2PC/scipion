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
import sys
import time
import psutil
import argparse

def main():
    parser = argparse.ArgumentParser(description='store memory and cpu usage, plot it using monitor.py program. Usage scipion run monitor.py')

    parser.add_argument('--interval', action="store", type=int, default=60, help="record during these minutes")
    parser.add_argument('--sleepSec', action="store", type=int, default=2, help="record each these seconds")
    args = parser.parse_args()
    # Get arguments. This will be implement later
    baseFn = os.path.join(os.environ.get('SCIPION_USER_DATA'),'tmp','log.sqlite')
    tableName = 'log'
    print "INFO: log data stored at sqlite file %s in table %s"%(baseFn,tableName)
    conn = lite.connect(baseFn, isolation_level=None)
    cur = conn.cursor()    
    createTable(conn,baseFn,tableName)
    interval = args.interval # logging during these minutes
    sleepSec = args.sleepSec # wait this seconds between logs
    loopPsUtils(cur,tableName,interval,sleepSec)

def createTable(cur,baseFn,tableName):

    cur.execute("DROP TABLE IF EXISTS %s"%tableName)
    sql = """CREATE TABLE %s(
                            id INTEGER PRIMARY KEY AUTOINCREMENT,
                            timestamp DATE DEFAULT (datetime('now','localtime')),
                            mem FLOAT,
                            cpu FLOAT,
                            openfile INT)"""%tableName
    cur.execute(sql)

def loopPsUtils(cur, tableName, interval,sleepSec):
     timeout = time.time() + 60.*interval   # interval minutes from now
     #ignore first meassure because is very unrealible
     psutil.cpu_percent(True)
     psutil.virtual_memory()
     while True:
         if time.time() > timeout:
            break
         #, percpu=True
         # non-blocking (percentage since last call)
         cpu = psutil.cpu_percent(interval=0)
         mem = psutil.virtual_memory()
         #vmem(total=8240947200L, available=7441436672L, percent=9.7, used=1939496960L, 
         #free=6301450240L, active=727162880, inactive=966086656, buffers=123904000L, cached=1016082432)
         memPercent = mem.percent
         sql = "INSERT INTO %s(mem,cpu) VALUES(%f,%f);"%(tableName,memPercent,cpu)
         cur.execute(sql)
         time.sleep(sleepSec)

if __name__ == "__main__":
    main()
