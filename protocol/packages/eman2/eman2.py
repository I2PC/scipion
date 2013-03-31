'''
Created on Mar 24, 2013

@author: airen
'''
from os import system
from pyworkflow.protocol import ProtParticlePicking
from optparse import OptionParser

class Eman2ProtParticlePicking(ProtParticlePicking):
    
    
    def __init__(self, micpath, boxsize, iswrite,**args):
        ProtParticlePicking.__init__(self, **args)
        print 'init'
        self.micpath = micpath
        self.boxsize = boxsize
        self.iswrite = iswrite
        self._run()
    
    def _run(self):
        print 'run'
        boxerargs = '%(micpath)s ' %self.micpath
        if not self.boxsize is None:
            boxerargs += ' --boxsize %(boxsize)s' %self.boxsize
        if iswrite:
            boxerargs += '--write_ptcls'
        
        system('e2boxer.py %s' %boxerags)
        
    
    
def parseInputAndRun():
        
#        try:
            print 'parseInputAndRun'
#            usage = "usage: eman2_prot_particle_picking [options] Example: eman2_prot_particle_picking *.mrc --boxsize 100 --iswrite true --workingdir ."
#            parser = OptionParser(usage)
#            parser.add_option("-mp", "--micpath", dest="micpath", default="/dev/stdin", type="string", help="micrographs path")
#            parser.add_option("-bs", "--boxsize", dest="boxsize", default="/dev/stdin", type="string", help="particle box size")
#            parser.add_option("-wp", "--write_ptcls", dest="iswrite", default="/dev/stdin", type="string", help="save particles to current dir")
#            parser.add_option("-wd", "--working_dir", dest="workingdir", default="/dev/stdin", type="string", help="working dir")
#            (options, args) = parser.parse_args()
#            eman2 = Eman2ProtParticlePicking(options.micpath, options.boxsize, options.iswrite, options.workingdir)
            eman2 = Eman2ProtParticlePicking('home/airen/xprojects/eman2/*.mrc', 100, True, workingDir='/home/airen/xprojects/eman2')
#        except:
#            print usage
#            print 'error'
            
parseInputAndRun()      