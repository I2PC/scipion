#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:    Airen Zaldivar         (airenzp@gmail.com) 
#               J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import sys
import matplotlib.pyplot as plt


if __name__ == '__main__':
    #TODO: REMOVE THIS AFTER DEBUGGING
    print "ARGS: "
    for i, arg in enumerate(sys.argv):
        print "%02d: %s" % (i, arg)

    dbName = sys.argv[1]
    dbPreffix = sys.argv[2]
    columns = sys.argv[3].split()
    colors = sys.argv[4].split()
    lines = sys.argv[5].split()
    markers = sys.argv[6].split()
    xcolumn = sys.argv[7]
    ylabel = sys.argv[8]
    xlabel = sys.argv[9]
    title = sys.argv[10]
    bins = sys.argv[11]



    from pyworkflow.mapper.sqlite import SqliteFlatDb
    db = SqliteFlatDb(dbName=dbName, tablePrefix=dbPreffix)
    setClassName = db.getProperty('self') # get the set class name
    from pyworkflow.em import getObjects
    setObj = getObjects()[setClassName](filename=dbName, prefix=dbPreffix)



    if xcolumn:
        xvalues = []

        for obj in setObj:
            if hasattr(obj, xcolumn):
                value = getattr(obj, xcolumn)
            elif xcolumn == 'id':
                id = int(obj.getObjId())
                print id
                xvalues.append(id)

    else:
        xvalues = range(0, setObj.getSize())

    i = 0
    for column in columns:
        yvalues = []

        for obj in setObj.iterItems(orderBy=column):
            if hasattr(obj, column):

                value = getattr(obj, column)
                yvalues.append(value.get())

        if bins:
            plt.hist(yvalues, bins=int(bins), color=colors[i])
        else:
            plt.plot(xvalues, yvalues, colors[i])
        i += 1

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

