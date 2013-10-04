#!/usr/bin/env xmipp_python

import os

from protlib_xmipp import XmippScript, createProgramsDb, createProgramsAutocomplete, ProgramKeywordsRank, getProgramsDbName,\
    getXmippLabels, blueStr, redStr
from protlib_filesystem import getXmippPath


class ScriptApropos(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        self.results = []
        
    def defineParams(self):
        self.addUsageLine("Search for Xmipp programs that are related to some keywords.")
        self.addUsageLine("Useful for when does not remember a program name.")
        ## params
        self.addParamsLine(" -i <...>          : Keyword list to search programs matching")
        self.addParamsLine("   alias --input;")
        self.addParamsLine("or -u               : Update the database with programs info")
        self.addParamsLine("   alias --update;")
        self.addParamsLine("or -l           : List all Xmipp programs or categories")
        self.addParamsLine("   alias --list;")    
        self.addParamsLine(" [-t <type=programs>]      : Type of operations")
        self.addParamsLine("    where <type> programs categories both labels")
        self.addParamsLine("   alias --type;")    
        ## examples
        self.addExampleLine("Search for program containing the keyword 'header'", False)
        self.addExampleLine("   xmipp_apropos -i header")
        self.addExampleLine("Search for program containing the keywords 'noise' and 'gaussian'", False)
        self.addExampleLine("   xmipp_apropos noise gaussian")
        self.addExampleLine("List all xmipp programs", False)
        self.addExampleLine("   xmipp_apropos --list")
        self.addExampleLine("List all xmipp metadata labels", False)
        self.addExampleLine("   xmipp_apropos --list -t labels")
            
    def readParams(self):
        self.keywords = []
        if self.checkParam('-i'):
            self.keywords = [k.lower() for k in self.getListParam('-i')]
        self.progRank = ProgramKeywordsRank(self.keywords)
        self.type = self.getParam('--type')
        
    def hasKeywords(self, label):
        for key in self.keywords:
            for v in label.values():
                if key in v.lower():
                    return True
        return False 
            
    def run(self):
        if self.checkParam('-u'):
            db = createProgramsDb()
            createProgramsAutocomplete()
        elif self.type == 'labels':
            labels = getXmippLabels()
            if self.checkParam("-i"):
                labels = [l for l in labels if self.hasKeywords(l)]
            # Check which labels are ported to python
            xm = __import__('xmipp') # load xmipp module object
            
            for l in labels:
                try:
                    getattr(xm, l['enum'])
                except Exception, ex:
                    l['comment'] = l['enum'] + ' (MISSING IN PYTHON BINDING)'
                print '{name:<30} {type:<30} {enum:<30} {comment}'.format(**l)
        else:
            dbName = getProgramsDbName()
            if not os.path.exists(dbName):
                db = createDb()
            else:
                from protlib_sql import ProgramDb
                db = ProgramDb(dbName)

            onlyList = self.checkParam('--list')
            
            if onlyList and self.type != 'programs':
                categories = db.selectCategories()
                doBoth = self.type == 'both'
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

