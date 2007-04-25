def init_log_system(projectdir,logdir,scriptname,WorkDirectory):
    import os, sys
    import logging
    import socket
    
    """ Set up logging information
    more info at:
    http://webmaster.iu.edu/tool_guide_info/python/lib/module-logging.html
    """
            
    LogName = projectdir + '/' + logdir + '/'
    if not os.path.exists(LogName):
        os.makedirs(LogName)
    scriptname=os.path.basename(scriptname)
    LogName += scriptname.replace('.py','')
    if not (WorkDirectory=="." or WorkDirectory=='.'):
        LogName += '_'
        LogName += WorkDirectory
    LogName += '.log'

    mylog = logging.getLogger(scriptname)
    hdlr = logging.FileHandler(LogName)
    formatter = logging.Formatter('%(levelname)s (%(lineno)d) %(message)s (%(asctime)s)')
    hdlr.setFormatter(formatter)
    mylog.addHandler(hdlr) 
    mylog.setLevel(logging.INFO)

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
    import shutil,os
    #strip_file_name
    protocol_name=str(os.path.basename(script_file_name))                        
    in_file_name=script_file_name
    out_file_name=absolute_path_to_working_dir +\
                   '/' + protocol_name + "_backup"
    shutil.copy(in_file_name,out_file_name)
