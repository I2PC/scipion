#!/usr/bin/env python

import os
from subprocess import Popen, PIPE
from protlib_xmipp import XmippScript
from protlib_filesystem import getXmippPath
from protlib_utils import blueStr, redStr
from protlib_sql import ProgramDb

class ScriptApropos(XmippScript):
	def __init__(self):
		XmippScript.__init__(self)
		self.weights = {'name': 5, 'keywords': 3, 'usage': 1}
		self.results = []
		
	def defineParams(self):
		self.addUsageLine("Search for Xmipp programs that are related to some keywords.")
		self.addUsageLine("Useful for when does not remeber a program name. gaussian noise header")
		## params
		self.addParamsLine(" -i <...>          : Keyword list to search programs matching")
		self.addParamsLine("   alias --input;")
		self.addParamsLine("or -u		       : Update the database with programs info")
		self.addParamsLine("   alias --update;")
		self.addParamsLine("or -l <list=programs>		          : List all Xmipp programs or categories")
		self.addParamsLine("    where <list> programs categories both")		
		self.addParamsLine("   alias --list;")	
		## examples
		self.addExampleLine("Search for program containing the keyword 'header'", False)
		self.addExampleLine("   xmipp_apropos -i header")
		self.addExampleLine("Search for program containing the keywords 'noise' and 'gaussian'", False)
		self.addExampleLine("   xmipp_apropos noise gaussian")
		self.addExampleLine("List all xmipp programs", False)
		self.addExampleLine("   xmipp_apropos --list")
			
	def readParams(self):
		self.keywords = []
		self.dbName = os.path.join(getXmippPath(), 'programs.sqlite')
		if self.checkParam('-i'):
			self.keywords = self.getListParam('-i')
			
	
	def getRank(self, program):
		rank = 0
		for k in self.keywords:
			for wkey, wvalue in self.weights.iteritems():
				if program[wkey].find(k) != -1:
					rank += wvalue
		return rank
	
	def run(self):
		if self.checkParam('-u'):
			db = createDb(self.dbName)
		else:
			if not os.path.exists(self.dbName):
				db = createDb(self.dbName)
			else:
				db = ProgramDb(self.dbName)

			onlyList = self.checkParam('--list')
			
			if onlyList and self.getParam('--list') != 'programs':
				categories = db.selectCategories()
				doBoth = self.getParam('--list') == 'both'
				for c in categories:
					print blueStr(c['name'])
					if doBoth:
						programs = db.selectPrograms(c)
						#self.maxlen = 50
						self.printPrograms(programs)
			else:
				results = []
				programs = db.selectPrograms()
				# Calculate ranks
				for p in programs:
					if onlyList:
						rank = 1
					else:
						rank = self.getRank(p)
					#Order by insertion sort
					if rank > 0:
						name = self.highlightStr(os.path.basename(p['name']))
						#self.maxlen = max(self.maxlen, len(name))
						pos = len(results)
						for i, e in reversed(list(enumerate(results))):
							if e['rank'] < rank:
								pos = i
							else: break

						results.insert(pos, {'rank': rank, 'name': name, 'usage': p['usage']})
				# Print results
				self.printPrograms(results)
		
	def printPrograms(self, programs):
		maxlen = 50
		for p in programs:
			name = p['name']
			desc = p['usage']
			if len(desc) > 0:
				desc = self.highlightStr(desc.splitlines()[0])
			print name.ljust(maxlen), desc	
						 
	def highlightStr(self, str):
		for k in self.keywords:
			str = str.replace(k, redStr(k))
		return str

def skipProgram(programName):
	if programName in ['xmipp_sqlite3',
					'xmipp_classify_CL2D_core_analysis',
					'xmipp_transform_threshold']:
		return True
	for p in ['xmipp_test', 'xmipp_template', 'mpi_']:
		if programName.find(p) != -1:
			return True
	return False

#Some helper functions
def createDb(dbName):
	db = ProgramDb(dbName)
	print 'Created db with name: %(dbName)s' % locals()
	db.create()
	#Create categories dictionary to classify programs
	#looking in program name and category prefixes
	categories = db.selectCategories()
	categoryDict = {}
	for c in categories:
		prefixes = c['prefixes'].split()
		for p in prefixes:
			categoryDict[p] = c
			
	import glob
	programs = glob.glob(os.path.join(getXmippPath('bin'), 'xmipp_*'))
	for p in programs:
		try:
			p = os.path.basename(p)
			if not skipProgram(p):
				if p.find('_mpi') != -1:
					p = "mpirun -np 2 %s" % p					
				cmd = [p, "--xmipp_write_definition"]
				ps = Popen(cmd, stdout=PIPE, stderr=PIPE)
				stdoutdata, stderrdata = ps.communicate()
				if stderrdata != '':
					raise Exception(stderrdata)
				for prefix, category in categoryDict.iteritems():
					if prefix in p:
						db.updateProgramCategory(p, category)
						break
		except Exception, e:
			print "PROGRAM: " + p
			print redStr("ERROR: " + str(e)) 
	return db
		
if __name__ == '__main__':
	ScriptApropos().tryRun()

