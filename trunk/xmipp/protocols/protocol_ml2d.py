#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from protlib_base import XmippProtocol, protocolMain
from protlib_utils import printLog   
from config_protocols import protDict

class ProtML2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.key, scriptname, project)
        self.Import = 'from protocol_ml2d import *'

    def validate(self):
        return []
        #return ["Protocol not implemented yet..."]
    
    def summary(self):
        return ["This is a test summary",
                "Need a real summary here..."
                ]
        
    def defineActions(self):
        print '*********************************************************************'
        progId = "ml"
        if (self.DoMlf):
            progId += "f"  
        
        program = "xmipp_%s_align2d" % progId
#        action = "Executing"
#        if (restart):
#            action = "Restarting"
#            
#        print '*  %s %s program :' % (action, program)
        restart = False
        if (restart):
            pass 
            #Not yet implemented
            #params= ' --restart ' + utils_xmipp.composeFileName('ml2d_it', RestartIter,'log')
        else: 
            params = ' -i %s --oroot %s/%s2d' % (self.InSelFile, self.WorkingDir, progId)
            # Number of references will be ignored if -ref is passed as expert option
            if self.ExtraParams.find("--ref") == -1:
                params += ' --nref %i' % self.NumberOfReferences
            params += ' ' + self.ExtraParams
            if (self.DoFast and not self.DoMlf):
                params += ' --fast'
            if (self.DoNorm):
                params += ' --norm'
            if (self.DoMirror):
                params += ' --mirror'
            if (self.NumberOfThreads > 1  and not self.DoMlf):
                params += ' --thr %i' % self.NumberOfThreads
            if (self.DoMlf):
                if (self.DoCorrectAmplitudes):
                    params += ' --ctfdat %s' % self.InCtfDatFile
                else:
                    params += ' --no_ctf --pixel_size %f' % self.PixelSize
                if (not self.ImagesArePhaseFlipped):
                    params += ' --not_phase_flipped'
                if (self.HighResLimit > 0):
                    params += ' --high %f' % self.HighResLimit
                    
        self.Db.insertAction('runJob', 
                             programname=program, 
                             params=params,
                             NumberOfMpiProcesses = self.NumberOfMpiProcesses,
                             NumberOfThreads = self.NumberOfThreads,
                             SystemFlavour = self.project.SystemFlavour)

