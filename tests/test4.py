#!/usr/bin/env xmipp_python

from pyworkflow.object import *
from pyworkflow.emx import *
from pyworkflow.mapper.xmlmapper import XmlMapper
from random import randint
import os

dbName = 'db.xml'

mapper = XmlMapper(dbName)

objList = mapper.select()

for o in objList:
    print o

