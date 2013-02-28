#!/usr/bin/env xmipp_python

import sys
from pyworkflow.emx import *
from pyworkflow.object import *
from pyworkflow.mapper.sqlite import SqliteMapper

mapper = SqliteMapper('db.sqlite', globals())

#c = mapper.get(2)
#print c.name
#print c.x.name
#print c.y.name
#c.x.set(100)
#c.y.set(200)

#mapper.store(c)

#mapper.commit()

l = mapper.select(classname='Array')[0]

#n = len(l)
#for i in range(n):
#    print l[i]
    
#print 'total: ', n

#c = l[n-1]
#c.x.set(999)
#c.y.set(999)

print l[0]
print l[1]

#mapper.store(c)
#mapper.commit()


#for k, v in l.__dict__.iteritems():
#    if k.startswith('item'):
#        print k, v

