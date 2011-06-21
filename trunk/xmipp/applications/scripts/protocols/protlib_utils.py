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

import os
import sys
'''number of leading zeros in filenames'''
FILENAMENUMBERLENTGH=6

#---------------------------------------------------------------------------
# Logging utilities
#---------------------------------------------------------------------------
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
        import socket
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

#---------------------------------------------------------------------------
# Naming conventions
#---------------------------------------------------------------------------
def getScriptPrefix(script):
    '''This function will extract the root of the protocol filename
    By example:
    xmipp_protocol_ml2d_00001.py -> xmipp_protocol_ml2d
    xmipp_protocol_ml2d.py       -> xmipp_protocol_ml2d
    
    script - script filename, could also contains path
    '''
    import re
    script = os.path.basename(script)
    #all protocols script should start by 'xmipp_protocol', the numbering part is optional
    s = re.match('((?:\w*[a-zA-Z])+)(?:_\d+)?.py', script)
    if not s:
        raise Exception('script %s doesn\'t conform Xmipp protocol name convention' % script)
    return s.group(1)

#---------------------------------------------------------------------------
# Other utilities
#---------------------------------------------------------------------------
def makeScriptBackup(log, script, WorkingDir):
    '''Make a backup of the script
    This function assumes the execution from ProjectDir
    '''
    log.info("Making backup of script " + script);
    import shutil, os
    try:        
        script_prefix = getScriptPrefix(script)
        script_out = os.path.join(WorkingDir, script_prefix + '_backup.py')
        shutil.copy(script, script_out)
    except shutil.Error, e:
        printLog(log, e.message)     
    
#---------------------------------------------------------------------------
# Parsing of arguments
#---------------------------------------------------------------------------
def getComponentFromVector(__vector, _iteration):
    ''' Convert a string to a vector of parameters'''
    _vector = __vector.strip()
    listValues = getListFromVector(_vector)
    if _iteration < 0: _iteration = 0
    if _iteration < len(listValues): return listValues[_iteration]
    else:                          return listValues[len(listValues) - 1]


def getListFromVector(_vector,numberIteration=None):
    ''' Convert a string to a list '''
    import string
    intervalos = string.split(_vector)
    if len(intervalos) == 0:
        raise RuntimeError, "Empty vector"
    listValues = []
    for i in range(len(intervalos)):
        intervalo = intervalos[i]
        listaIntervalo = string.split(intervalo, 'x')
        if len(listaIntervalo) == 1:
            listValues += listaIntervalo
        elif len(listaIntervalo) == 2:
            listValues += [ listaIntervalo[1] ] * string.atoi(listaIntervalo[0])
        else:
            raise RuntimeError, "Unknown syntax: " + intervalos
    #fill with last value the iterations
    if( numberIteration):
        for i in range(len(listValues),numberIteration):
            listValues.append(listValues[-1])
        
    return listValues

def getBoolListFromVector(_vector,numberIteration=None):
    ''' Convert a string to a list of booleans'''
    listValues = getListFromVector(_vector,numberIteration)
    listValuesBool = []
    for i in range(len(listValues)):
        if listValues[i]=='0':
            listValuesBool.append(False)
        else:
            listValuesBool.append(True)
    return listValuesBool

#---------------------------------------------------------------------------
# Error handling
#---------------------------------------------------------------------------
def reportError(msg):
    '''Function to write error message to log, to screen and exit'''
    print >> sys.stderr, "ERROR: ", msg
    exit(1)
        
def printLogError(log, msg):
    '''Function to write error message to log, to screen and exit'''
    log.error(msg)
    print >> sys.stderr, "ERROR: ", msg
    exit(1)
    
def printLog(log, msg):
    '''Just print a msg and log'''
    log.info(msg)
    print msg
    
#---------------------------------------------------------------------------
# Jobs launching
#---------------------------------------------------------------------------    
# The job should be launched from the working directory!
def runJob(log, 
           programname,
           params,
           DoParallel,
           NumberOfMpiProcesses,
           NumberOfThreads,
           SystemFlavour):
    
    command = buildRunCommand(log,
               programname,
               params,
               DoParallel,
               NumberOfMpiProcesses,
               NumberOfThreads,
               SystemFlavour)
    printLog(log, "Running command: %s" % command)

    from subprocess import call
    retcode = 0
    try:
        retcode = call(command, shell=True)
        printLog(log, "Process returned with code %d" % retcode)
    except OSError, e:
        printLogError(log, "Execution failed %s" % e)

    return (command, retcode)

def buildRunCommand(
               log,
               programname,
               params,
               DoParallel,
               NumberOfMpiProcesses,
               NumberOfThreads,
               SystemFlavour):
    if not DoParallel:
        command = programname + ' ' + params
    else:
        paramsDict['prog'] = programname.replace('xmipp', 'xmipp_mpi')
        paramsDict['jobs'] = NumberOfMpiProcesses
        paramsDict['params'] = params
        
        if (SystemFlavour == 'SLURM-MPICH'): # like BSCs MareNostrum, LaPalma etc
            mpicommand = 'srun '
        elif (SystemFlavour == 'TORQUE-OPENMPI'): # like our crunchy
            mpicommand = 'mpirun -mca mpi_yield_when_idle 1 -np %(jobs)d'
            if (int(NumberOfThreads) > 1):
                mpicommand += ' --bynode'
        elif (SystemFlavour == 'SGE-OPENMPI'): # like cluster at imp.ac.at (no variable nr_cpus yet...)
            mpicommand = 'mpiexec -n  %(jobs)d' 
        elif (SystemFlavour == 'PBS'): # like in Vermeer and FinisTerrae
            paramsDict['file']  = os.environ.get('PBS_NODEFILE')
            mpicommand = 'mpirun -np  %(jobs)d -hostfile %(file)s'
        elif (SystemFlavour == 'XMIPP_MACHINEFILE'): # environment variable $XMIPP_MACHINEFILE points to machinefile
            paramsDict['file']  = os.environ.get('XMIPP_MACHINEFILE')
            mpicommand = 'mpirun -np  %(jobs)d -machinefile %(file)s'
        elif (SystemFlavour == 'HOME_MACHINEFILE'): # machinefile is called $HOME/machines.dat
            paramsDict['file'] = os.environ.get('HOME') + '/machinefile.dat'
            mpicommand = 'mpirun -np   %(jobs)d -machinefile %(file)s'
        elif (SystemFlavour == ''):
            mpicommand = 'mpirun -mca mpi_yield_when_idle 1 -np %(jobs)d'
        else:
            printLogError(log, 'Unrecognized SystemFlavour %s' % SystemFlavour)
        command = (mpicommand + ' `which %(prog)s` %(params)s') % paramsDict
    return command

def loadModule(moduleName):
    try:
        if moduleName in sys.modules:
            module = sys.modules[moduleName]
            reload(module)
        else:
            module = __import__(moduleName)
    except ImportError, e:
        reportError(str(e))
    return module
        
def loadLaunchConfig(configFile):
    moduleName = configFile.replace('.py', '')
    try:
        if moduleName in sys.modules:
            module = sys.modules[moduleName]
            reload(module)
        else:
            module = __import__(moduleName)
    except ImportError, e:
        reportError(str(e))
        
def createLaunchFile(outFilename, fileTemplate, **kargs):
    file = open(outFilename, 'w')
    file.write(fileTemplate % kargs)
    file.close()
    
def launchJob(cmdTemplate, jobFile):
    command = cmdTemplate % {'file':jobFile}
    printLog(log, "Lauching job: %s" % command)
    from subprocess import call
    retcode = 0
    try:
        retcode = call(command, shell=True)
    except OSError, e:
        reportError("Launch failed %s" % e)


#---------------------------------------------------------------------------
# Metadata stuff
#--------------------------------------------------------------------------- 
#create a metadata file with original image name, and two other 
#lines with variation over the original name
def intercalate_union_3(inFileName,outFileName, src1,targ1,src2,targ2):
   
   mD = MetaData(inFileName)
   mDout = MetaData()
   
   for id in mD:       
       idOut = mDout.addObject()
       sIn = mD.getValue(MDL_IMAGE,id)
       mDout.setValue(MDL_IMAGE, sIn, idOut)
       enabled= mD.containsLabel(MDL_ENABLED)
       
       if  (enabled):
       
            i = int(mD.getValue(MDL_ENABLED,id))
            mDout.setValue(MDL_ENABLED, i, idOut)
       
       idOut = mDout.addObject()

       ss = sIn.replace(src1,targ1)
       mDout.setValue(MDL_IMAGE, ss, idOut)
       
       if  (enabled):
           mDout.setValue(MDL_ENABLED, i, idOut)
           
       idOut = mDout.addObject()
       
       ss = sIn.replace(src2,targ2)
       mDout.setValue(MDL_IMAGE, ss, idOut)
       
       if  (enabled):
           mDout.setValue(MDL_ENABLED, i, idOut)
       
   mDout.write(outFileName)

#set rot and tilt between -180,180 and -90,90
def check_angle_range(inFileName,outFileName):

    mD    = MetaData(inFileName)
    doWrite=False
    
    for id in mD: 
        doWrite2=False
        rot = mD.getValue(MDL_ANGLEROT,id)
        tilt = mD.getValue(MDL_ANGLETILT,id)
        if tilt > 90.: 
            tilt = -(int(tilt)-180)
            rot  += 180.
            doWrite=True
            doWrite2=True
        if tilt < -90.: 
            tilt = -(int(tilt)+180)
            rot  -= 180. 
            doWrite=True
            doWrite2=True
        if (doWrite2):
            mD.setValue(MDL_ANGLEROT , rot, id)
            mD.setValue(MDL_ANGLETILT, tilt, id)
        
    if(doWrite or inFileName != outFileName):
        mD.write(outFileName)


#compute histogram
def compute_histogram(mD,bin,col,min,max):
    
   allMD = MetaData()
   outMD = MetaData()   
   _bin = (max-min)/bin
   
   for h in range(0,bin):
       outMD.removeObjects(MDQuery("*"))
       if (h==0):
           outMD.importObjects(mD, MDValueRange(col, float(min), float(_bin*(h + 1)+min)))
       if (h>0 and h<(bin-1)):
           outMD.importObjects(mD, MDValueRange(col, float(_bin * h + min), float(_bin*(h + 1)+min)))
       if (h==(bin-1)):
           outMD.importObjects(mD, MDValueRange(col, float(_bin * h + min), float(max)))
       
       _sum=float(outMD.aggregateSingle(AGGR_SUM,MDL_WEIGHT))
       outMD.addLabel(MDL_COUNT)
       outMD.setValueCol(MDL_COUNT, int(_sum+0.1))
       allMD.unionAll(outMD)
       
   return allMD

#########################
# FileName Handling
###########################
def unique_filename(file_name):
    ''' Create a unique filename (not file handler)
       this approach is unsecure but good enought for most purposes'''
    import os
    counter = 1
    file_name_parts = os.path.splitext(file_name) # returns ('/path/file', '.ext')
    while os.path.isfile(file_name):
        file_name = file_name_parts[0] + '_' + str(counter) + file_name_parts[1]
        counter += 1
    return file_name 
