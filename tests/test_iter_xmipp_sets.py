'''
Created on May 8, 2013

@author: antonio
'''
import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *
from em.packages.xmipp3.data import XmippSetOfMicrographs, XmippCoordinate

projName = sys.argv[1]
manager = Manager()
proj = manager.createProject(projName)  # Now it will be loaded if exists

print ('*************************************************************************************************************************')
result = proj.mapper.selectByClass('XmippSetOfMicrographs')
if len(result):    
    for xmippSetOfMicrographs in result:
        for xmippMicrograph in xmippSetOfMicrographs:
            print ('Micrograph: ' + xmippMicrograph.getFileName())
else:
    print "Not XmippSetOfMicrographs found"
print ('*************************************************************************************************************************')

result = proj.mapper.selectByClass('XmippSetOfCoordinates')
if len(result):    
    for xmippSetOfCoordinates in result:
        for xmippCoordinate in xmippSetOfCoordinates.iterCoordinates():
            print ("Coordinate: " + str(xmippCoordinate.getPosition()))
else:
    print "Not XmippSetOfCoordinates found"
print ('*************************************************************************************************************************')

print ('*************************************************************************************************************************')
result = proj.mapper.selectByClass('XmippSetOfImages')
if len(result):    
    for xmippSetOfImages in result:
        for xmippImage in xmippSetOfImages:
            print ('Image: ' + xmippImage.getFileName())
else:
    print "Not XmippSetOfMicrographs found"
print ('*************************************************************************************************************************')