#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
import os

class ProtDummy(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.dummy.key, scriptname, project)
        self.Import = 'from protocol_dummy import *'

    def validate(self):
        return []
        #return ["Protocol not implemented yet..."]
    
    def summary(self):
        return ["This is a test summary",
                "Need a real summary here..."
                ]
        
    def defineSteps(self):
        print '*********************************************************************'
        self.Db.insertStep('runJob', programname="xmipp_apropos", params="-k fourier")
        tmpFile = os.path.join(self.TmpDir, "kk.txt")
        self.Db.insertStep('createTempFile', filename=tmpFile)
        

def createTempFile(log, filename):
    f = open(filename, 'w')
    f.write("Esto es una prueba")
    f.close()

