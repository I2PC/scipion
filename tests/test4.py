#!/usr/bin/env xmipp_python

from pyworkflow.object import *
from pyworkflow.emx import *
from pyworkflow.mapper.xmlmapper import XmlMapper
from random import randint
import os

fileName = 'db.xml'
emxData   = EmxData()
mapper    = XmlMapper(emxData)
mapper.read(fileName)
mapper.convertToEmxData(emxData)
#testing
print emxData

