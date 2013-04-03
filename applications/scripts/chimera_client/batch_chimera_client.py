#!/usr/bin/env xmipp_python

from protlib_chimera import XmippChimeraClient, XmippProjectionExplorer
from optparse import OptionParser
from os import system
from os.path import exists
from protlib_filesystem import getXmippPath

class BatchXmippChimeraClient:

	def __init__(self):
		self.parseInput()
		serverfile = getXmippPath('libraries/bindings/chimera/xmipp_chimera_server.py')
		system("chimera %s  &" % serverfile)
		#print 'running client'
		print self.mode
		if self.mode == 'projector':
			XmippProjectionExplorer(self.volfile)
		elif self.mode == 'viewer':
			client = XmippChimeraClient(self.volfile)
			client.listen()
	def parseInput(self):
		
		try:
			self.usage = "usage: %prog [options] Example: %prog -i hand.vol --mode projector"
			self.parser = OptionParser(self.usage)
			self.parser.add_option("-i", "--input", dest="volfile", default="/dev/stdin", type="string", help="Volume to display")
			self.parser.add_option("-m", "--mode", dest="mode", default="/dev/stdin", type="string", help="client mode: viewer for visualization and projection for projection explorer")
			(options, args) = self.parser.parse_args()
			if options.volfile == '/dev/stdin' or not(exists(options.volfile)):#simple validation
				raise ValueError(options.volfile)
			
			self.volfile = options.volfile
			if options.mode == '/dev/stdin':
				self.mode = 'viewer'
			else:
				self.mode = options.mode
		except:
			print self.usage
			exit()
			
  



if __name__ == '__main__':
	BatchXmippChimeraClient()