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
import sys
import time
import argparse
import sqlite3 as lite
from matplotlib import pyplot



def getData():
    """ Retrieve several lists with values to plot. """
    # Get arguments. This will be implement later
    baseFn = os.path.join(os.environ.get('SCIPION_USER_DATA'),'tmp','log.sqlite')
    tableName = 'log'
    conn = lite.connect(baseFn, isolation_level=None)
    cur = conn.cursor()
    #conn.text_factory = str

    #I guess this can be done in a single call
    cur.execute("select julianday(timestamp)  from %s where id=1"%tableName )
    initTime = cur.fetchone()[0]
    cur.execute("select timestamp  from %s where id=1"%tableName )
    initTimeTitle = cur.fetchone()[0]
    cur.execute("select (julianday(timestamp) - %f)*24  from %s"%(initTime,tableName) )
    idValues = [r[0] for r in cur.fetchall()] 
    
       
    def get(name):
        cur.execute("select %s from %s" % (name, tableName))
        return [r[0] for r in cur.fetchall()]
    
    data = {'initTime': initTime,
            'initTimeTitle': initTimeTitle,
            'idValues': idValues,
            'cpu': get('cpu'),
            'mem': get('mem'),
            'swap': get('swap'),
            'disks_read_per_sec': get('disks_read_per_sec'),
            'disks_write_per_sec': get('disks_write_per_sec'),
            'disks_read_time_sec': get('disks_read_time_sec'),
            'disks_write_time_sec': get('disks_write_time_sec'),
            }
    
    conn.close()
    return data


    
def main():
    parser = argparse.ArgumentParser(description='Plot cpu and mem usage. Usage: scipion run ./plotter.py [--cpu] [--mem] [--swap]')
    parser.add_argument('--cpu', '-c',
          action='store_true',
          help='Plot cpu usage percentage stored by monitor.' )
    parser.add_argument('--mem', '-m',
          action='store_true',
          help='Plot virtual memory stored by monitor' )
    parser.add_argument('--swap', '-s',
                        action='store_true',
                        help='Plot swap memory stored by monitor.')
    parser.add_argument('--disks_read_per_sec', '--drb',
                        action='store_true',
                        help='Report disk usage read bytes seg')
    parser.add_argument('--disks_write_per_sec', '--dwb',
                        action='store_true',
                        help='Report disk usage write byte seg')
    parser.add_argument('--disks_read_time_sec', '--drt',
                        action='store_true',
                        help='Report disk usage read seg (do not think is working)')
    parser.add_argument('--disks_write_time_sec', '--dwt',
                        action='store_true',
                        help='Report disk usage write seg (do not think is working)')
    #TODO: I have remove the repaint ability in plotter since it make imposible to maximize the plot
    # and it is not very useful. If you want to recover it uncomment all the lines starting with #! 
    #and comment the lines ending with #!  
    #!parser.add_argument('--sleepSec', action="store", type=int, default=20, help="record each these seconds")
    
    args = parser.parse_args()   
    
    # Plot memory by default if not option is selected
    if not (args.cpu or args.mem or args.swap or args.disks_read_per_sec or args.disks_write_per_sec or args.disks_read_time_sec or args.disks_write_time_sec):
         args.mem = True

    # Force to start the x values at 0
    #x1,x2,y1,y2 = pyplot.axis()
    #pyplot.axis((0.,x2,y1,y2))
    
    # Show plots
    
    pyplot.xlabel("time (hours)")
    pyplot.ylabel("percentage (or Mb for IO)")
    
    
    #pyplot.ion()
    
    lines = {}
    
    
    while True:
        data = getData()
        idValues = data['idValues']
        
        def plot(name):
            if getattr(args, name):
                
                values = data[name]
                
                if name in lines:
                    line = lines[name]
                    line.set_xdata(idValues)
                    line.set_ydata(values)
                else:
                    line, = pyplot.plot(idValues, values, label=name)#,'-b', label=name)
                    lines[name] = line
        plot('cpu')
        plot('mem')
        plot('swap')
        plot('disks_read_per_sec') 
        plot('disks_write_per_sec') 
        plot('disks_read_time_sec') 
        plot('disks_write_time_sec') 
        pyplot.title("System Monitor (%s)" % data['initTimeTitle'])
        pyplot.legend()
        x1,x2,y1,y2 = pyplot.axis()
        pyplot.axis((0.,x2,y1,y2))
        #!pyplot.draw()
        pyplot.show()#!
        #!time.sleep(args.sleepSec)
        break        



if __name__ == "__main__":
    main()
