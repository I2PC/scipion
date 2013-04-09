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
        XmippScript.__init__(self, True)
        
    def defineParams(self):
        self.addUsageLine('Chimera client for visualization and projection of xmipp volumes')
        ## params
        self.addParamsLine('[ --input <input>]          : Volume to visualize')
        self.addParamsLine('   alias -i;')
        self.addParamsLine('[ --mode <mode=viewer>]             : Sets visualization mode')

        self.addParamsLine(' where <mode>')
        self.addParamsLine('  viewer         : Allows volume visualization')
        self.addParamsLine('  projector <padding_factor="1">  <max_freq="0.5"> <spline_degree="BSPLINE2">   : Allows volume visualization and projection. spline_degree can be: NEAREST, LINEAR, BSPLINE2, BSPLINE3 and BSPLINE4.')
        self.addParamsLine('   alias -m;')
        
        
        self.addParamsLine('[ --angulardist <angulardist>  <color=red> <spheres_distance="default"> <spheres_maxradius="default">]     : Volume angular distribution to visualize')
        self.addParamsLine('   alias -a;')
        
        self.addExampleLine('Opens xmipp chimera client in projector mode:', False)
        self.addExampleLine('xmipp_chimera_client -i hand.vol --mode projector', False)
            
    def run(self):
    	
    	volfile = self.getParam('-i')
    	mode = self.getParam('-m')
    	angulardistfile = self.getParam('-a')
        
        spheres_color = self.getParam('-a', 1)
        spheres_distance = self.getParam('-a', 2)
        spheres_maxradius = self.getParam('-a', 3)
        
        
        isprojector = (mode == 'projector')
        if isprojector:
            padding_factor = self.getDoubleParam('-m', 1)
            max_freq = self.getDoubleParam('-m', 2)
            
            spline_degree_str = self.getParam('-m', 3)
            if spline_degree_str == 'NEAREST':
                spline_degree = NEAREST
            elif spline_degree_str == 'LINEAR':
                spline_degree = LINEAR
            elif spline_degree_str == 'BSPLINE2':
                spline_degree = BSPLINE2
            elif spline_degree_str == 'BSPLINE3':
                spline_degree = BSPLINE3
            elif spline_degree_str == 'BSPLINE4':
                spline_degree = BSPLINE4
		              
       
		

        serverfile = getXmippPath('libraries/bindings/chimera/xmipp_chimera_server.py')
        system("chimera %s  &" % serverfile)
        if isprojector:
			XmippProjectionExplorer(volfile, angulardistfile, spheres_color, spheres_distance, spheres_maxradius, padding_factor, max_freq, spline_degree)
			print 'created projection explorer'
        elif mode == 'viewer':
			client = XmippChimeraClient(volfile, angulardistfile, spheres_color, spheres_distance, spheres_maxradius)
			client.listen()
			print 'created chimera client'
if __name__ == '__main__':
    ScriptChimeraClient().tryRun()
    
    

