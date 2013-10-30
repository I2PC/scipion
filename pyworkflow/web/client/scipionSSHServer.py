'''
Created on Sep 30, 2013

@author: adrian
'''
from bottle import route, run, template, request
import subprocess

@route('/openInteractiveProtocolMenu')
def openInteractiveProtocolMenu():
    commands, machine, user = readParametersFromScipion()
    password = getPassword(machine,user)
    
        
    return template('interactive_protocol_menu.tpl',
                    commands=(';').join(commands),
                    machine=machine,
                    user=user,
                    password=password)

@route('/runInteractiveProtocol', method='POST')
def runInteractiveProtocol():
    print "runningInteractiveProtocol"
    commands, machine, user, password = readParametersLocal()
    callScript(commands, machine, user, password)
    
#    return template('<b>Hello {{name}}</b>!', name="Name")
    return template('<b>Protocol Executed</b>!')


def readParametersFromScipion():
    print "readingparameter"
    
    machine=request.query.get("machine","takarras.cnb.csic.es")
    user=request.query.get("user","aquintana")
    commands=request.query.getall("commands")
    
    print "Execute ", commands, " as ", user," at ", machine,    
        
    return commands, machine, user 

def readParametersLocal():
    print "readingparameterLocal"
    machine = request.forms.get('machine')
    user = request.forms.get('user')
    password = request.forms.get('password')
    commands = request.forms.get('commands')
    return commands, machine, user, password  
    
def getPassword(machine,user):
    print "gettingpassword"
    return ""
    
def callScript(commands, machine, user, password):     
   
    print "callingscript"

    shellCommand="bash managePipe "+password
    print "shellCommand",shellCommand
    subprocess.call(shellCommand, shell=True) 
    
    shellCommand="bash sshServerScript "+machine +" "+user+" "+password+" \""+commands+"\""
    print "shellCommand",shellCommand
    subprocess.call(shellCommand, shell=True) 
    
    

run(host='localhost', port=8081)