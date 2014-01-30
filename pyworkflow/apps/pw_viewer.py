#!/usr/bin/env python


from protlib_xmipp import ScriptShowJ

class ScriptScipionViewer(ScriptShowJ):
    def __init__(self):
        ScriptShowJ.__init__(self, 'xmipp.viewer.scipion.ScipionViewer')
        
    def defineOtherParams(self):
        ScriptShowJ.defineOtherParams(self)
        self.addParamsLine('  --command  <button> <script> <projectid> <protid> <dbpath>  : Command button will be included on the GUI to execute script associated');
       
        
    def readOtherParams(self):
         ScriptShowJ.readOtherParams(self)
         cmdbutton = self.getParam('--command', 0)
         cmdscript = self.getParam('--command', 1)
         projectid = self.getParam('--command', 2)
         protid = self.getParam('--command', 3)
         dbpath = self.getParam('--command', 4)
         self.args += " --command %(cmdbutton)s %(cmdscript)s %(projectid) %(protid) %(dbpath)" %locals()
        
if __name__ == '__main__':
    ScriptScipionViewer().tryRun()