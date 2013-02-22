#!/usr/bin/env xmipp_python

from pyworkflow.object import *
from pyworkflow.mapper.sqlite import SqliteMapper
from random import randint
import os

dbName = 'db.sqlite'
size = 1000000

print "Removing the db..."
os.system('rm %s' % dbName)
mapper = SqliteMapper(dbName)
#c = mapper.get(2)

#c = Coordinate(10, 5)

#mapper.store(a)
#mapper.store(c)
#l = Array(size)
#
#print "Populating list..."
#for i in range(size):
#    x = randint(0, 999)
#    y = randint(0, 999)
#    c = Coordinate(x, y)
#    l[i] = c

m = Micrograph()
m.Particles[0] = Coordinate(10, 10)
m.Particles[1] = Coordinate(1,3)
m.Path.set('/home/jose/micrograph.xmp')
m.N.set(2)

mapper.store(m)

    
#for k, v in l.__dict__.iteritems():
#    if k.startswith('item'):
#        print k, v
        
#print "len(l): ", len(l)
print "Storing..."
#mapper.store(l)

mapper.commit()