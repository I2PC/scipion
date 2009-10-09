import os, sys
import socket
class XmippLog:
    def __init__(self,logName):
        self.fh_log=open(logName,'a')
    
    def debug(self,message):
        self.info("DEBUG: "+message)

    def info(self,message):
        import time
        self.fh_log.write(message+" "+\
           time.asctime(time.localtime(time.time()))+"\n")
        self.fh_log.flush()

    def __del__(self):
        self.fh_log.close()

# --------------------------------------------------------------------------
# Init log system
# --------------------------------------------------------------------------
def init_log_system(projectdir,logdir,scriptname,WorkDirectory):
    """ Set up logging information
    more info at:
    http://webmaster.iu.edu/tool_guide_info/python/lib/module-logging.html
    """

    if logdir[0]=='/':
       LogName = logdir
    else:
       LogName = projectdir + '/' + logdir
    if not LogName[-1]=='/':
       LogName+='/'
    if not os.path.exists(LogName):
        os.makedirs(LogName)
    scriptname=os.path.basename(scriptname)
    LogName += scriptname.replace('.py','')
    if not (WorkDirectory=="." or WorkDirectory=='.'):
        LogName += '_'
        LogName += os.path.basename(WorkDirectory)
    LogName += '.log'

    try:
       import logging
       mylog = logging.getLogger(scriptname)
       hdlr = logging.FileHandler(LogName)
       formatter = logging.Formatter('(%(asctime)s) %(levelname)s (%(lineno)4d) %(message)s')
       hdlr.setFormatter(formatter)
       mylog.addHandler(hdlr) 
       mylog.setLevel(logging.INFO)
    except ImportError:
       mylog = XmippLog(LogName)
       

    # append a line with user, machine and date
    myusername = str(os.environ.get('USERNAME'))
    myhost = str(socket.gethostname())
    mypwd = str(os.environ.get('PWD'))
    event = "\n"
    event += "NEW LOG SESSION\n" 
    event += "===============\n" 
    event += myusername + '@' 
    event += myhost + ':' 
    event += mypwd 
    mylog.info(event)

    return mylog
            
def make_backup_of_script_file(script_file_name,
                               absolute_path_to_working_dir):
    import shutil,os,re
    #check is temporal name, i.e. ends with _ddddd.py
    #where d is a digit 
    in_file_name=script_file_name
    script_file_name=re.sub('_\d\d\d\d\d\d.py$','.py',in_file_name)

    protocol_name=str(os.path.basename(script_file_name)).replace('.py','')
    out_file_name=absolute_path_to_working_dir +\
                   '/' + protocol_name + "_backup.py"
    shutil.copy(in_file_name,out_file_name)
