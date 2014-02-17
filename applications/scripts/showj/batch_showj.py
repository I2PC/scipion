#!/usr/bin/env xmipp_python


from protlib_xmipp import ScriptShowJ
        if self.checkParam('--label_relion') or self.getParam('-i').endswith('.star'):
            from protlib_import import relionLabelString
            os.environ['XMIPP_EXTRA_ALIASES'] = relionLabelString()
        
if __name__ == '__main__':
    ScriptShowJ().tryRun()

