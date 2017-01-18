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
# *  e-mail address 'scipion@cnb.csic.es'
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
    parser.add_argument('--sleepSec', action="store", type=int, default=20, help="record each these seconds")
    parser.add_argument('--doDisk', action='store_true', default=False, help="record disk IO information")

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
    doDisk = args.doDisk
    loopPsUtils(cur,tableName,interval,sleepSec,doDisk)

def createTable(cur,baseFn,tableName):

    cur.execute("DROP TABLE IF EXISTS %s"%tableName)
    sql = """CREATE TABLE %s(
                            id INTEGER PRIMARY KEY AUTOINCREMENT,
                            timestamp DATE DEFAULT (datetime('now','localtime')),
                            mem FLOAT,
                            cpu FLOAT,
                            disks_read_per_sec FLOAT,
                            disks_write_per_sec FLOAT,
                            disks_read_time_sec FLOAT,
                            disks_write_time_sec FLOAT,
                            openfile INT,
                            swap FLOAT)""" % tableName
    cur.execute(sql)

def bytes2human(n):
    # http://code.activestate.com/recipes/578019
    # >>> bytes2human(10000)
    # '9.8K'
    # >>> bytes2human(100001221)
    # '95.4M'
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n


def loopPsUtils(cur, tableName, interval,sleepSec,doDisk):
     timeout = time.time() + 60.*interval   # interval minutes from now
     #ignore first meassure because is very unrealible
     psutil.cpu_percent(True)
     psutil.virtual_memory()
     disks_before = psutil.disk_io_counters(perdisk=False)    
     

     while True:
         if time.time() > timeout:
            break
         #, percpu=True
         # non-blocking (percentage since last call)
         cpu = psutil.cpu_percent(interval=0)
         mem = psutil.virtual_memory()
         swap = psutil.swap_memory()
         
         if doDisk:
            disks_after = psutil.disk_io_counters(perdisk=False)
            disks_read_per_sec   = (disks_after.read_bytes  - disks_before.read_bytes)/(sleepSec * 1024.*1024.) 
            disks_write_per_sec  = (disks_after.write_bytes - disks_before.write_bytes)/(sleepSec * 1024.*1024.) 
            disks_read_time_sec  = (disks_after.read_time   - disks_before.read_time)/(sleepSec*1000.)
            disks_write_time_sec = (disks_after.write_time  - disks_before.write_time)/(sleepSec*1000.)

            #print "read %fM, write=%fM, read_time=%fs,  write_time =%fs"%((disks_read_per_sec),\
            #                                                            (disks_write_per_sec),\
            #                                                            (disks_read_time_sec),\
            #                                                            (disks_write_time_sec) )
            disks_before = disks_after
         
         #vmem(total=8240947200L, available=7441436672L, percent=9.7, used=1939496960L, 
         #free=6301450240L, active=727162880, inactive=966086656, buffers=123904000L, cached=1016082432)


            sql = """INSERT INTO %s(mem,cpu,swap,
                                    disks_read_per_sec,
                                    disks_write_per_sec,
                                    disks_read_time_sec,
                                    disks_write_time_sec) VALUES(%f,%f,%f,%f,%f,%f,%f);""" % (tableName, mem.percent, cpu, swap.percent,disks_read_per_sec,disks_write_per_sec,disks_read_time_sec,disks_write_time_sec )
         else:
            sql = """INSERT INTO %s(mem,cpu,swap) VALUES(%f,%f,%f);""" % (tableName, mem.percent, cpu, swap.percent)
         cur.execute(sql)
         time.sleep(sleepSec)

if __name__ == "__main__":
    main()
