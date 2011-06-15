'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''


class XmippLog:
    '''This class is a simple wrapper around the most rich Python logging system
    Also providing a basic file logging for older python versions
    '''    
    def __init__(self, logname, filename):
        self.is_basic = False
        try:
            import logging
            mylog = logging.getLogger(logname)
            hdlr = logging.FileHandler(filename)
            formatter = logging.Formatter('(%(asctime)s) %(levelname)s (%(lineno)4d) %(message)s')
            hdlr.setFormatter(formatter)
            mylog.addHandler(hdlr) 
            mylog.setLevel(logging.INFO)
            self._log = mylog
        except ImportError:
            self._logfile = open(filename, 'a')
            self.is_basic = True
        # append a line with user, machine and date
        import os, socket
        myusername = str(os.environ.get('USERNAME'))
        myhost = str(socket.gethostname())
        mypwd = str(os.environ.get('PWD'))
        event = "\n"
        event += "NEW LOG SESSION\n" 
        event += "===============\n" 
        event += myusername + '@' 
        event += myhost + ':' 
        event += mypwd 
        self.info(event)
        
    def debug(self, message):
        if self.is_basic:
            self.info("DEBUG: " + message)
        else:
            self._log.debug(message)

    def info(self, message):
        if self.is_basic:
            import time
            self.fh_log.write("%s %s\n" % (message, time.asctime(time.localtime(time.time()))))
            self.fh_log.flush()
        else:
            self._log.info(message)
            
    def cat(self, filename):
        '''Cat a file contents into a log'''
        fh = open(filename, 'r')
        self.info("# Content of " + filename)
        for line in fh:
            self.info("#       " + line.strip())
        fh.close()       

    def __del__(self):
        if self.is_basic:
            self.fh_log.close()


def getScriptPrefix(script):
    '''This function will extract the root of the protocol filename
    By example:
    xmipp_protocol_ml2d_00001.py -> xmipp_protocol_ml2d
    xmipp_protocol_ml2d.py       -> xmipp_protocol_ml2d
    
    script - script filename, could also contains path
    '''
    import re
    import os
    script = os.path.basename(script)
    #all protocols script should start by 'xmipp_protocol', the numbering part is optional
    s = re.match('(xmipp_protocol(?:_[a-zA-Z]\w*[a-zA-Z])+)(?:_\d+)?.py', script)
    if not s:
        raise Exception('script %s doesn\'t conform Xmipp protocol name convention' % script)
    return s.group(1)

def makeScriptBackup(script, absolute_path_to_working_dir):
    import shutil, os
    script_prefix = getScriptPrefix(script)
    script_out = os.path.join(absolute_path_to_working_dir, script_prefix + '_backup.py')
    shutil.copy(script, script_out)

