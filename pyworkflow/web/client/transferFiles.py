'''
Created on Sep 30, 2013

@author: adrian
'''
import subprocess
from os.path import join
from scipionClientUtils import *
from bottle import route, run, template, request

def openInteractiveProtocolMenu():
    inputParameters = readParametersFromScipion()
    password = getPassword(inputParameters['machine'], inputParameters['user'])
    
        
    forwardTemplate = template(join(path["templates"], templates["error"]), inputParameters = inputParameters)
        
    elif inputParameters['action'] == "TransferFilesMenu":
        forwardTemplate =  template(join(path["templates"],templates["transferFiles"]),
                                    commands = (';').join(inputParameters['commands']),
                                    machine = inputParameters['machine'],
                                    user = inputParameters['user'],
                                    password = password)        
    return forwardTemplate

@route('/runInteractiveProtocol', method='POST')
def runInteractiveProtocol():
    print "Running Interactive Protocol"
    inputParameters = readParametersFromLocal()
    
    callInteractiveProtocolScript(inputParameters)
    
#    return template('<b>Hello {{name}}</b>!', name="Name")
    return template('<b>Protocol Executed</b>!')

@route('/runTransferFiles', method='POST')
def runTransferFiles():
    print "Running Transfer Files"
    inputParameters = readParametersFromLocal()
    
    validateHost(inputParameters)
    existFiles(inputParameters)
    createTunnel(inputParameters)
    buildFilesList(inputParameters)
    createDir(inputParameters)
    
    callTransferFilesScript(inputParameters)

    return template('<b>Protocol Executed</b>!')

def createDir(inputParameters):
    print "Creating Dir"
    return ""

def buildFilesList(inputParameters):
    print "Building Files List"
    return ""

def createTunnel(inputParameters):
    print "Creating Tunnel (if GW parameters are provided)"
    return ""

def existFiles(inputParameters):
    print "Checking if files exist"
    return ""

def validateHost(inputParameters):
    print "Validating Host"
    print "Testing Connection"
    print "Testing Rsync"
    
    return ""

def readParametersFromScipion():
    print "Reading Parameters From Scipion..."
    
    inputParameters={}
    inputParameters['action'] = request.query.get("action") #InteractiveProtocolMenu || TransferFilesMenu || RunInteractiveProtocol || RunTransferFiles 
    
    inputParameters['machine'] = request.query.get("machine")
    inputParameters['user'] = request.query.get("user")
    inputParameters['commands'] = request.query.getall("commands")
    
    inputParameters['numberTrials'] = request.query.get("numberTrials",3)
    inputParameters['pathAskpass'] = request.query.get("pathAskpass","scipionPipe")
    
    
    if inputParameters['action'] == "TransferFilesMenu": 
        inputParameters['mode'] = request.query.get("mode","COPY")
        inputParameters['filesSource'] = request.query.getall("filesSource")
        inputParameters['sourcePath'] = request.query.get("sourcePath",None) #If None files Source will be absolute path
        inputParameters['targetPath'] = request.query.get("targetPath")
        inputParameters['canRsync'] = request.query.get("canRsync","true")
        inputParameters['availableTunnel'] = request.query.get("availableTunnel","true")
        inputParameters['portTunnel'] = request.query.get("portTunnel",None)
    
    print "Parameters: ", inputParameters
    return inputParameters 

def readParametersFromLocal():
    print "Reading Parameters From Local..."
    
    inputParameters={}
    inputParameters['action'] = request.forms.get("action") #InteractiveProtocolMenu || TransferFilesMenu || RunInteractiveProtocol || RunTransferFiles 
    
    inputParameters['machine'] = request.forms.get("machine")
    inputParameters['user'] = request.forms.get("user")
    inputParameters['commands'] = request.forms.getall("commands")
    inputParameters['password'] = request.forms.get('password')
    
    inputParameters['numberTrials'] = request.forms.get("numberTrials",3)
    inputParameters['pathAskpass'] = request.forms.get("pathAskpass","scipionPipe")
    

    if inputParameters['action'] == "RunTransferFiles": 
        inputParameters['mode'] = request.forms.get("mode","COPY")
        inputParameters['filesSource'] = request.forms.getall("filesSource")
        inputParameters['sourcePath'] = request.forms.get("sourcePath",None) #If None files Source will be absolute path
        inputParameters['targetPath'] = request.forms.get("targetPath")
        inputParameters['canRsync'] = request.forms.get("canRsync","true")
        inputParameters['availableTunnel'] = request.forms.get("availableTunnel","true")
        inputParameters['portTunnel'] = request.forms.get("portTunnel",None)
    
    print "Parameters: ", inputParameters
    return inputParameters 
    
    
def getPassword(machine,user):
    print "Getting Password. TBD: Store password locally"
    return ""
    
def callInteractiveProtocolScript(inputParameters):     
    print "Calling Interactive Protocol Script"
    managePipe(inputParameters['password'])
    
    shellCommand="bash " + join(path["scripts"],scripts["interactiveProtocol"]) + " " + \
        inputParameters['machine'] + " " + inputParameters['user'] + " " + \
        " \"" + inputParameters['commands'] + "\" " + inputParameters['pathAskpass']
    print "shellCommand ", shellCommand
    subprocess.call(shellCommand, shell=True) 

def callTransferFilesScript(inputParameters):     
    print "Calling Transfer Files Script"
    managePipe(inputParameters['password'])
    
    shellCommand="bash " + join(path["scripts"],scripts["transferFiles"]) + " " + inputParameters['commands'] + " " + \
        inputParameters['remoteuser'] + " " + inputParameters['pathAskpass'] + " " + \
        inputParameters['fileList'] + " " + inputParameters['sourcePath'] + " " + \
        inputParameters['targetPath'] + " " + inputParameters['canRsync'] + " " + \
        inputParameters['availableTunnel'] + " " + inputParameters['portTunnel'] 
    print "shellCommand ", shellCommand
    subprocess.call(shellCommand, shell=True)    
  
def managePipe(password):
    print "Managing Pipe"
    shellCommand="bash " + join(path["scripts"],scripts["pipe"]) + " " + password
    print "Shell Command: ", shellCommand
    subprocess.call(shellCommand, shell=True)       
    

run(host='localhost', port=8081)