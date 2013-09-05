import unittest
    
    
class TestWorkflow(unittest.TestCase):
    
    protDict = {}
        
    def validateFiles(self, key, prot):
        """ Validate if the produced files are the expected ones.
        Params:
            prot: the protocol to validate. 
            filesSet: the known files that should be produced (set)
        """
        #FIXME: This is disabled until changes of CTF and Relations
        return True
        self.protDict[key] = prot
        filesSet = self.getProtocolFiles(key)
        self.assertEqual(prot.getFiles(), filesSet)
        
    def printSet(self, msg, s):
        print "============= %s ==========" % msg
        for i in s:
            print i
            
    def getProtocolFiles(self, key):
        fileList = self.GOLD_FILES[key]
        fileSet = set([self.__replaceFilename(f) for f in fileList])
        
        return fileSet
    
    def __replaceFilename(self, filename):
        """ Convert list to set and replace the key
        in the filename by the protocol working dir. 
        """
        for k, v in self.protDict.iteritems():
            if filename.startswith(k):
                return filename.replace(k, v.getWorkingDir())
        return filename
                

