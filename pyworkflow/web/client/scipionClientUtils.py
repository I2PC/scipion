'''
Created on Sep 30, 2013

@author: adrian
'''
import subprocess
from os.path import join


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
    shellCommand="bash " + join(path["scripts"],scripts["pipe"]) + " " + pathAskpass + " " + password
    print "Shell Command: ", shellCommand
    subprocess.call(shellCommand, shell=True)   
    
def getPassword(machine,user):
    print "Getting Password. TBD: Store password locally"
    return ""        
    

