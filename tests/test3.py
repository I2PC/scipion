#!/usr/bin/env xmipp_python

from pyworkflow.object import *
from pyworkflow.emx import *
from pyworkflow.mapper.xmlmapper import XmlMapper
from random import randint
import os

dbName = 'db.xml'
size = 1000000

print "Removing the xml..."
os.system('rm %s' % dbName)

rootName = 'EMX'
header = '''
  ##########################################################################
  #               EMX Exchange file 
  #               Produced by the prolib_emx module
  # 
  #  This is a EMX file.
  #
  #  Information on this file format is available at 
  #  http://i2pc.cnb.csic.es/emx
  ##########################################################################
  #  One of the best ways you can help us to improve this software
  #  is to let us know about any problems you find with it.
  #  Please report bugs to: emx@cnb.csic.es
  ##########################################################################
  '''
version = 1.0
mapper = XmlMapper(dbName, rootName, version, header)

m1 = Micrograph(id={'filename': 'mic', 'index': '1'})
m1.acceleratingVoltage.set(100)
m1.defocusU.set(1000.)
m1.pixelSpacing.X.set(5.6)
m1.pixelSpacing.Y.set(5.7)

m2 = Micrograph(id={'filename': 'mic', 'index': '2'})
m2.acceleratingVoltage.set(200)
m2.activeFlag.set(None)
m2.defocusUAngle.set(135)
#m.Particles[0] = EmxCoordinate(10, 10)
#m.Parti3cles[1] = EmxCoordinate(1, 3)
p1 = Particle(id={'filename': 'parti', 'index': '1'})
p1.boxSize.X.set(1)
p1.boxSize.Y.set(3)
p1.pixelSpacing.X.set(55.6)
p1.pixelSpacing.Y.set(55.7)
p1.defocusU.set(1000.)
p1.micrograph.set(m1)

#print "eeeeeeeeeeeeeeeee", p1.micrographFK.id

mapper.insert(m1)
mapper.insert(m2)
mapper.insert(p1)
mapper.commit()