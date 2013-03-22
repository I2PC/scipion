#!/usr/bin/env xmipp_python

from protlib_chimera import XmippProjectionExplorer
from optparse import OptionParser
from os import system
from os.path import exists
from protlib_filesystem import getXmippPath

class BatchXmippChimeraClient:

	def __init__(self):
		self.volfile = self.parseInput()
		serverfile = getXmippPath('libraries/bindings/chimera/xmipp_chimera_server.py')
		system("chimera %s&" % serverfile)
		#print 'running client'
		XmippProjectionExplorer(self.volfile)


	def parseInput(self):
		
		try:
			self.usage = "usage: %prog [options] Example: %prog -i hand.vol"
			self.parser = OptionParser(self.usage)
			self.parser.add_option("-i", "--input", dest="volfile", default="/dev/stdin", type="string", help="Volume to display")
			(options, args) = self.parser.parse_args()
			if options.volfile == '/dev/stdin' or not(exists(options.volfile)):#simple validation
				raise ValueError(options.volfile)
			return (options.volfile)
		except:
			print self.usage
			exit()
			
  



if __name__ == '__main__':
	BatchXmippChimeraClient()