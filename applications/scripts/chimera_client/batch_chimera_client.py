#!/usr/bin/env xmipp_python

from protlib_chimera import XmippChimeraClient, XmippProjectionExplorer
from optparse import OptionParser
from os import system
from os.path import exists
from protlib_filesystem import getXmippPath
import sys
from protlib_xmipp import XmippScript
from xmipp import *


class ScriptChimeraClient(XmippScript):
    def __init__(self):
        XmippScript.__init__(self, False)
        
    def defineParams(self):
        self.addUsageLine('Chimera client for visualization and projection of xmipp volumes')
        ## params
        self.addParamsLine('[ --input <input>]          : Volume to visualize')
        self.addParamsLine('   alias -i;')
        self.addParamsLine('[ --mode <mode=viewer>]             : Sets visualization mode')

        self.addParamsLine(' where <mode>')
        self.addParamsLine('  viewer         : Allows volume visualization')
        self.addParamsLine('  projector <size="default"> <padding_factor="1">  <max_freq="0.5"> <spline_degree="BSPLINE2">   : Allows volume visualization and projection. spline_degree can be: NEAREST, LINEAR, BSPLINE2, BSPLINE3 and BSPLINE4.')
        self.addParamsLine('   alias -m;')
        
        
        self.addParamsLine('[ --angulardist <angulardist=none>  <color=red> <spheres_distance="default"> <spheres_maxradius="default">]     : Volume angular distribution to visualize')
        self.addParamsLine('   alias -a;')
        
        self.addExampleLine('Opens xmipp chimera client in projector mode:', False)
        self.addExampleLine('xmipp_chimera_client -i hand.vol --mode projector', False)
            
    def run(self):
    	
    	volfile = self.getParam('-i')
        if not exists(volfile):
            print "ERROR: File " + volfile + "does not exist\n"
            sys.exit(1)

    	mode = self.getParam('-m')
    	angulardistfile = self.getParam('-a')
        if angulardistfile == 'none':
            angulardistfile = None
        else:
            if not exists(angulardistfile):#either file does not exists or has a blockname
                if not(existsBlockInMetaDataFile(angulardistfile)):#check blockname
                    print "ERROR: File " + angulardistfile + "does not exist or block is missing\n"
                    sys.exit(1)
        
        spheres_color = self.getParam('-a', 1)
        spheres_distance = self.getParam('-a', 2)
        spheres_maxradius = self.getParam('-a', 3)
        
        
        isprojector = (mode == 'projector')
        if isprojector:
            size = self.getParam('-m', 1)
            padding_factor = self.getDoubleParam('-m', 2)
            max_freq = self.getDoubleParam('-m', 3)
            
            spline_degree_str = self.getParam('-m', 4)
            if spline_degree_str.lower() == 'NEAREST'.lower():
                spline_degree = NEAREST
            elif spline_degree_str.lower() == 'LINEAR'.lower():
                spline_degree = LINEAR
            elif spline_degree_str.lower() == 'BSPLINE2'.lower():
                spline_degree = BSPLINE2
            elif spline_degree_str.lower() == 'BSPLINE3'.lower():
                spline_degree = BSPLINE3
            elif spline_degree_str.lower() == 'BSPLINE4'.lower():
                spline_degree = BSPLINE4
#            print spline_degree
		              
       
		
        port = self.getFreePort()
        if not port:
            print "ERROR: Port is not available\n"
            sys.exit(1)
        serverfile = getXmippPath('libraries/bindings/chimera/xmipp_chimera_server.py')
        command = "export XMIPP_CHIMERA_PORT=%d; chimera %s  &" % (port,serverfile)
        system(command)
        if isprojector:
			XmippProjectionExplorer(volfile, port, angulardistfile, spheres_color, spheres_distance, spheres_maxradius, size, padding_factor, max_freq, spline_degree)
#			print 'created projection explorer'
        elif mode == 'viewer':
			client = XmippChimeraClient(volfile, port, angulardistfile, spheres_color, spheres_distance, spheres_maxradius)
			client.listen()
#			print 'created chimera client'

    
    def getFreePort(self,basePort=0,host=''):
        import socket
        port=0
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.bind((host, basePort))
            ipaddr, port = s.getsockname()
            s.close()
        except Exception, e:
            print e
            return 0
        return port
    
if __name__ == '__main__':
    ScriptChimeraClient().tryRun()
    
    

