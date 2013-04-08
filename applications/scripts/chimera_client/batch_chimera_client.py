#!/usr/bin/env xmipp_python

from protlib_chimera import XmippChimeraClient, XmippProjectionExplorer
from optparse import OptionParser
from os import system
from os.path import exists
from protlib_filesystem import getXmippPath
import sys
from protlib_xmipp import XmippScript


class ScriptChimeraClient(XmippScript):
    def __init__(self):
        XmippScript.__init__(self, True)
        
    def defineParams(self):
        self.addUsageLine('Chimera client for visualization and projection of xmipp volumes')
        ## params
        self.addParamsLine('[ --input <input>]          : Volume to visualize')
        self.addParamsLine('   alias -i;')
        self.addParamsLine('[ --mode <mode="viewer">]             : Sets visualization mode')
        self.addParamsLine('   alias -m;')
        self.addParamsLine('[ --angulardist <angulardist="angulardist.xmd">]     : Volume angular distribution to visualize')
        self.addParamsLine('   alias -a;')
        
        self.addExampleLine('Open xmipp chimera client in projector mode:', False)
        self.addExampleLine('xmipp_chimera_client -i hand.vol --mode projector', False)
            
    def run(self):
    	
    	volfile = self.getParam('-i')
    	mode = self.getParam('-m')
    	angulardistfile = self.getParam('-a')
		
        if volfile is None or not(exists(volfile)):
        	raise ValueError(volfile)
		if options.mode is None:
			mode = 'viewer'
		if not angulardist is None:
			if not(exists(angulardistfile)):
				raise ValueError(angulardistfile)

        serverfile = getXmippPath('libraries/bindings/chimera/xmipp_chimera_server.py')
        system("chimera %s  &" % serverfile)
        print mode
        if mode == 'projector':
			XmippProjectionExplorer(volfile, angulardistfile)
			print 'created projection explorer'
        elif mode == 'viewer':
			client = XmippChimeraClient(volfile, angulardistfile)
			client.listen()
			print 'created chimera client'
if __name__ == '__main__':
    ScriptChimeraClient().tryRun()
    
    

