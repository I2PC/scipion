'''
Created on May 8, 2013

@author: antonio
'''
import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from em.packages.xmipp3.data import XmippSetOfMicrographs

projName = sys.argv[1]
manager = Manager()
proj = manager.createProject(projName)  # Now it will be loaded if exists


result = proj.mapper.selectByClass('XmippSetOfMicrographs')
if len(result):    
    for xmippSetOfMicrographs in result:
        print ("XmippSetOfMicrographs files: " + str(xmippSetOfMicrographs.getFiles()))
else:
    print "Not XmippSetOfMicrographs found"

result = proj.mapper.selectByClass('XmippSetOfCoordinates')
if len(result):    
    for xmippSetOfCoordinates in result:
        print ("XmippSetOfCoordinates files: " + str(xmippSetOfCoordinates.getFiles()))
else:
    print "Not XmippSetOfCoordinates found"

result = proj.mapper.selectByClass('EmanSetOfCoordinates')
if len(result):    
    for emanSetOfCoordinates in result:
        print ("EmanSetOfCoordinates files: " + str(emanSetOfCoordinates.getFiles()))
else:
    print "Not EmanSetOfCoordinates found"
    
result = proj.mapper.selectByClass('XmippCTFModel')
if len(result):    
    for xmippCTFModel in result:
        print ("XmippCTFModel files: " + str(xmippCTFModel.getFiles()))
else:
    print "Not XmippCTFModel found"