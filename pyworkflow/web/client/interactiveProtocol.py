'''
Created on Sep 30, 2013

@author: adrian
'''
import subprocess
from os.path import join
from scipionClientUtils import *
from bottle import route, run, template, request

def openInteractiveProtocolMenu():
    print "Opening Interactive Protocol Menu"
    
    inputParameters = readParametersFromScipion()
    password = getPassword(inputParameters['machine'], inputParameters['user'])
        
    if inputParameters['action'] == "openInteractiveProtocolMenu":        
        return  template(join(path["templates"],templates["interactiveProtocol"]),
                                    commands = (';').join(inputParameters['commands']),
                                    machine = inputParameters['machine'],
                                    user = inputParameters['user'],
                                    password = password,
                                    numberTrials = inputParameters['numberTrials'],
                                    pathAskpass = inputParameters['pathAskpass'])
    else:
        return template(join(path["templates"], templates["error"]), inputParameters = inputParameters)

def runInteractiveProtocol():
    print "Running Interactive Protocol"
    inputParameters = readParametersFromLocal()
    
    callInteractiveProtocolScript(inputParameters)
    
#    return template('<b>Hello {{name}}</b>!', name="Name")
    return template('<b>Protocol Executed</b>!')



def readParametersFromScipion():
    print "Reading Parameters From Scipion..."
    
    inputParameters={}
    inputParameters['action'] = request.query.get("action") #InteractiveProtocolMenu || TransferFilesMenu || RunInteractiveProtocol || RunTransferFiles 
    
    inputParameters['machine'] = request.query.get("machine")
    inputParameters['user'] = request.query.get("user")
    inputParameters['commands'] = request.query.getall("commands")
    
    inputParameters['numberTrials'] = request.query.get("numberTrials",3)
    inputParameters['pathAskpass'] = request.query.get("pathAskpass","scipionPipe")

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

    print "Parameters: ", inputParameters
    return inputParameters 
    
def callInteractiveProtocolScript(inputParameters):     
    print "Calling Interactive Protocol Script"
    managePipe(inputParameters['pathAskpass'], inputParameters['password'])
    
    shellCommand="bash " + join(path["scripts"],scripts["interactiveProtocol"]) + " " + \
        inputParameters['machine'] + " " + inputParameters['user'] + " " + \
        " \"" + inputParameters['commands'] + "\" " + inputParameters['pathAskpass']
    print "shellCommand ", shellCommand
    subprocess.call(shellCommand, shell=True) 
    
