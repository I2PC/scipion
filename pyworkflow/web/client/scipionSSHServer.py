'''
Created on Sep 30, 2013

@author: adrian
'''
from bottle import route, run, template, request
import subprocess

@route('/runInteractiveProtocol')
def runInteractiveProtocol():
    machine, user, commands, isXTerm = readParameters()
    password = getPassword(machine,user)
    callScript(machine,user,password,commands)
    
    #subprocess.call("bash sshServerScript", shell=True) 
 
    
    return template('<b>Hello {{name}}</b>!', name="Name")


def readParameters():
    machine=request.query.get("machine","pitagoras.cnb.csic.es")
    user=request.query.get("user","aquintana")
    commands=request.query.getall("commands")
    isXTerm=request.query.get("isXTerm",False)
    
    print "Execute ", commands, " with Xterm:", isXTerm, " as ", user," at ", machine,    
        
    return machine, user, commands, isXTerm
    
def getPassword(machine,user):
    print "gettingpassword"
    return "password"
    
def callScript(machine,user,pw,commands):     
    print "callingscript"
    subprocess.call("bash sshServerScript", shell=True) 

run(host='localhost', port=8081)