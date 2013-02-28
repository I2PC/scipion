#!/usr/bin/env xmipp_python


from pyworkflow.mapper.sqlite import SqliteMapper
from pyworkflow.emx import *
from pyworkflow.object import Array

from random import randint
import os

dbName = 'db.sqlite'
size = 2#0#00000

print "Removing the db..."
os.system('rm %s' % dbName)
mapper = SqliteMapper(dbName)
#c = mapper.get(2)

#c = Coordinate(10, 5)

#mapper.store(a)
#mapper.store(c)
l = Array(size)
units = ['px', 'mm']
#
print "Populating list..."
for i in range(size):
    x = randint(0, 999)
    y = randint(0, 999)
    c = BoxSize(X=x, Y=y, Z=y)
    c.setUnits(units[i])
    print "iter ", i
    print "vector: ", c
    #    c = Integer(x)#EmxVector(Integer, X=x, Y=y)

    l[i] = c


mapper.store(l)

    
#for k, v in l.__dict__.iteritems():
#    if k.startswith('item'):
#        print k, v
        
#print "len(l): ", len(l)
print "Storing..."
#mapper.store(l)

mapper.commit()
