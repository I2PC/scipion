#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles, according to:
#
# Example use:
# ./xmipp_protocol_ml2d.py
#
#  Author:  Sjors Scheres, January 2008
# Updated:  J. M. de la Rosa Trevin July 2011
#

import os, sys, shutil
from xmipp import MetaData
from protlib_base import *

class ProtML2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.key, scriptname, project)
        self.Import = 'from xmipp_protocol_ml2d import *'
#    def runSetup(self):
#        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
#        sys.path.append(scriptdir) # add default search path
#        # store script name
#        protocolName = sys.argv[0]
#        #assume restart
#        restart = True
#        # This is not a restart
#        if (RestartIter < 1):
#            # Delete working directory if it exists, make a new one
#            if (DoDeleteWorkingDir and DoML2D): 
#                if (self.WorkingDir==""):
#                   raise RuntimeError,"No working directory given"
#                if os.path.exists(WorkingDir):
#                    shutil.rmtree(WorkingDir)
#            if not os.path.exists(WorkingDir):
#                os.makedirs(WorkingDir)
#
#
#            # Create a selfile with absolute pathname in the WorkingDir
##            mysel = MetaData(InSelFile);
##            mysel.makeAbsPath();
##            InSelFile = os.path.abspath(os.path.join(WorkingDir, InSelFile))
##            mysel.write(InSelFile)
#
#            if (DoMlf and DoCorrectAmplitudes):
#                # Copy CTFdat to the workingdir as well
#                shutil.copy(InCtfDatFile, os.path.join(WorkingDir, 'my.ctfdat'))
#
#            # Backup script
#            log.make_backup_of_script_file(protocolName, os.path.abspath(WorkingDir))
#            restart = False
#
#        # Store current directory before moving
#        currentDir = os.getcwd()            
#        # Execute protocol in the working directory
#        os.chdir(WorkingDir)
#        self.run(restart)
#        # Return to parent dir
#        os.chdir(currentDir)

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


if __name__ == '__main__':
    protocolMain(ProtML2D)
