'''
Created on May 8, 2013

@author: antonio
'''
import sys

from pyworkflow.manager import Manager
from pyworkflow.project import *
from pyworkflow.em import *

projName = sys.argv[1]
manager = Manager()
proj = manager.createProject(projName)  # Now it will be loaded if exists

print ('*************************************************************************************************************************')
result = proj.mapper.selectByClass('EmanSetOfCoordinates')
if len(result):    
    for emanSetOfCoordinates in result:
        for emanCoordinate in emanSetOfCoordinates.iterCoordinates():
            print ("Coordinate: " + str(emanCoordinate.getPosition()))
else:
    print "Not EmanSetOfCoordinates found"
print ('*************************************************************************************************************************')