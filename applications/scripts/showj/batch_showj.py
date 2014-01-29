#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import ScriptAppIJ
from xmipp import Program, FileName
import xmipp

class ScriptShowJ(ScriptAppIJ):
    def __init__(self):
        ScriptAppIJ.__init__(self, 'xmipp.viewer.Viewer')
        
    def defineOtherParams(self):
        self.addParamsLine('  [--mode <mode_value=image>]           : List of params ')
        self.addParamsLine('     where <mode_value> image gallery metadata rotspectra')
        self.addParamsLine('         alias -d;')
        self.addParamsLine('  [--poll]                            : Keeps checking for changes on input files  (for image mode only!)')
        self.addParamsLine('         alias -p;')
        self.addParamsLine('  [--render <label=first>]    : Activates images rendering (for metadata mode only)')
        self.addParamsLine('                               : you can pass which label to render, by default the first one that can be visualized')
        self.addParamsLine('         alias -e;')
        self.addParamsLine('  [--rows <rows>]                            : number of rows in table')
        self.addParamsLine('         alias -r;')
        self.addParamsLine('  [--columns <columns>]                            : number of columns in table')
        self.addParamsLine('         alias -c;')
        self.addParamsLine('  [--zoom <zoom>]                            : zoom for images')
        self.addParamsLine('         alias -z;')
        self.addParamsLine('  [--view <axis="z">]                        : Viewer position (for volumes only)')
        self.addParamsLine('     where <axis> z y x z_pos y_pos x_pos')
        self.addParamsLine('  [--dont_apply_geo]                        : Does not read geometrical information(for metadata only)')
        self.addParamsLine('  [--dont_wrap]                             : Does not wrap (for metadata only)')
        self.addParamsLine('  [--debug] : debug')
        self.addParamsLine('  [--mask_toolbar] : Open mask toolbar (only valid for images)')
        self.addParamsLine('  [--label_alias <alias_string>]  : Activate some metadata label alias, for example')
        self.addParamsLine('                                  : anglePsi=aPsi;shiftX=sX;shiftY:sY')
        self.addParamsLine('  [--label_relion]                : Activates the mapping to Relion labels')
        
    def readOtherParams(self):
        params = ['--mode', '--rows', '--columns', '--zoom', '--view', '--render']
        for p in params:
            if self.checkParam(p):
                self.args += " %s %s" % (p, self.getParam(p))
        params = ['--poll', '--debug', '--dont_apply_geo', '--dont_wrap', '--mask_toolbar']
        for p in params:
            if self.checkParam(p):
                self.args += " %s" % p
                
        # Set environment var for extra label alias
        if self.checkParam('--label_alias'):
            os.environ['XMIPP_EXTRA_ALIASES'] = self.getParam('--label_alias')

        if self.checkParam('--label_relion'):
            from protlib_import import relionLabelString
            os.environ['XMIPP_EXTRA_ALIASES'] = relionLabelString()
        
if __name__ == '__main__':
    ScriptShowJ().tryRun()

