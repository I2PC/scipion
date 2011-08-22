#!/usr/bin/env python

import os

from protlib_xmipp import XmippScript, createProgramsDb, skipProgram, ProgramKeywordsRank, getProgramsDbName
from protlib_filesystem import getXmippPath
from protlib_utils import blueStr, redStr


class ScriptApropos(XmippScript):
	def __init__(self):
		XmippScript.__init__(self)
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
		if self.checkParam('-i'):
			self.keywords = self.getListParam('-i')
		self.progRank = ProgramKeywordsRank(self.keywords)
	
	def run(self):
		if self.checkParam('-u'):
			db = createProgramsDb()
		else:
			dbName = getProgramsDbName()
			if not os.path.exists(dbName):
				db = createDb()
			else:
				from protlib_sql import ProgramDb
				db = ProgramDb(dbName)

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
					rank = self.progRank.getRank(p)
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

if __name__ == '__main__':
	ScriptApropos().tryRun()

