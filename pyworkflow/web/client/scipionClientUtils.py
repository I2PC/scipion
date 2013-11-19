'''
Created on Sep 30, 2013

@author: adrian
'''
import subprocess
from os.path import join, basename, dirname
from os import mkfifo, remove, rmdir
from tempfile import mkdtemp

path = {"templates": "templates",
        "scripts": "scripts"}

scripts = {"pipe": "managePipeScript",
           "transferFiles": "transferFilesScript",
           "interactiveProtocol": "interactiveProtocolScript"}

templates = {"interactiveProtocol": "interactive_protocol_menu.tpl",
             "transferFiles":"transfer_files_menu.tpl",
             "error": "error_page.tpl"}

  
def managePipe(pathAskpass, password):
    print "Managing Pipe"
    
    tmpdir = mkdtemp()
    filename = join(tmpdir, pathAskpass)

    shellCommand="bash " + join(path["scripts"],scripts["pipe"]) + " " + filename + " " + password
    print "Shell Command: ", shellCommand
    subprocess.call(shellCommand, shell=True)
    
    return filename
#    tmpdir = mkdtemp()
#    filename = join(tmpdir, pathAskpass)
#    print filename
#    try:
#        mkfifo(filename)
#    except OSError, e:
#        print "Failed to create FIFO: %s" % e
#    else:
#        fifo = open(filename, 'w')
#        # write stuff to fifo
#        print >> fifo, "hello"
#        fifo.close()
#        
#        return filename
#    return None    
      
def removePipe(filename):
    remove(basename(filename))
    rmdir(dirname(filename)) 
    
def getPassword(machine,user):
    print "Getting Password. TBD: Store password locally"
    return ""        
    

