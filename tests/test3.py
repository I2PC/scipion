#!/usr/bin/env xmipp_python

from pyworkflow.object import *
from pyworkflow.emx import *
from pyworkflow.mapper.xmlmapper import XmlMapper
from random import randint
import os

fileName = 'db.xml'
print "Removing the xml..."
os.system('rm %s' % fileName)


m1 = micrograph(id={'filename': 'mic', 'index': '1'})
m1.acceleratingVoltage.set(100)
m1.defocusU.set(1000.)
m1.pixelSpacing.X.set(5.6)
m1.pixelSpacing.Y.set(5.7)

m2 = micrograph(id={'filename': 'mic', 'index': '2'})
m2.acceleratingVoltage.set(200)
m2.activeFlag.set(None)
m2.defocusUAngle.set(135)

p1 = particle(id={'filename': 'parti', 'index': '1'})
p1.boxSize.X.set(1)
p1.boxSize.Y.set(3)
p1.pixelSpacing.X.set(55.6)
p1.pixelSpacing.Y.set(55.7)
p1.defocusU.set(1000.)
#p1.setTransformationMatrix.t11.set(11.)
mat = {'t11':11.1,'t12':12.1,'t13':13.1,'t14':14.1,
       't21':21.1,'t24':24.1,
       't31':31.1,'t32':32.1,'t33':33.1,'t34':34.1,
       }
p1.transformationMatrix.setValue(**mat)
p1.setMicrograph(m1)


emxData   = EmxData()
emxData.addObject(m1)
emxData.addObject(m2)
emxData.addObject(p1)
mapper    = XmlMapper(emxData)
mapper.emxDataToXML()
mapper.write(fileName)