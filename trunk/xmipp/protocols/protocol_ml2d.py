#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

import os, sys, shutil
from xmipp import MetaData
from protlib_base import *

class ProtML2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.key, scriptname, project)
        self.Import = 'from xmipp_protocol_ml2d import *'

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
        if (DoMlf):
            progId += "f"  
        
        program = "xmipp_%s_align2d" % progId
#        action = "Executing"
#        if (restart):
#            action = "Restarting"
#            
#        print '*  %s %s program :' % (action, program)
        restart = False
        if (restart):
            params= ' --restart ' + utils_xmipp.composeFileName('ml2d_it', RestartIter,'log')
        else: 
            params = ' -i %s --oroot %s2d' % (InSelFile, progId)
            # Number of references will be ignored if -ref is passed as expert option
            if ExtraParamsMLalign2D.find("--ref") == -1:
                params += ' --nref %i' % NumberOfReferences
            params += ' ' + ExtraParamsMLalign2D
            if (DoFast and not DoMlf):
                params += ' --fast'
            if (DoNorm):
                params += ' --norm'
            if (DoMirror):
                params += ' --mirror'
            if (NumberOfThreads > 1  and not DoMlf):
                params += ' --thr %i' % NumberOfThreads
            if (DoMlf):
                if (DoCorrectAmplitudes):
                    params += ' --ctfdat my.ctfdat'
                else:
                    params += ' --no_ctf -pixel_size %f' + PixelSize
                if (not self.ImagesArePhaseFlipped):
                    params += ' --not_phase_flipped'
                if (self.HighResLimit > 0):
                    params += ' --high %f' + HighResLimit
                    
        self.Db.insertAction('launchML', program=program, params=params)
    
def launchML(log, program, params):
    #launchJob(program, params, self.log, DoParallel, NumberOfMpiProcesses, NumberOfThreads, SystemFlavour)
    print "Running program: '%s %s'" % (program, params)
    log.info("Running program: '%s %s'" % (program, params))

