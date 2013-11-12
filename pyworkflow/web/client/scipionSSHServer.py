'''
Created on Sep 30, 2013

@author: adrian
'''
from bottle import route, run, template, request
import subprocess
from os.path import join

templatesPath = "templates"
scriptsPath = "scripts"

@route('/openScipionLocalMenu')
def openScipionLocalMenu():
    inputParameters = readParameters()
    password = getPassword(inputParameters['machine'], inputParameters['user'])
    
        
    forwardTemplate = template(join(templatesPath, 'error_page.tpl'), inputParameters = inputParameters)
        
    if inputParameters['action'] == "InteractiveProtocolMenu":        
        forwardTemplate =  template(join(templatesPath,'interactive_protocol_menu.tpl'),
                                    commands = (';').join(inputParameters['commands']),
                                    machine = inputParameters['machine'],
                                    user = inputParameters['user'],
                                    password = password)
    elif inputParameters['action'] == "TransferFilesMenu":
        forwardTemplate =  template(join(templatesPath,'transfer_file_menu.tpl'),
                                    commands = (';').join(inputParameters['commands']),
                                    machine = inputParameters['machine'],
                                    user = inputParameters['user'],
                                    password = password)
        
    return forwardTemplate

@route('/runInteractiveProtocol', method='POST')
def runInteractiveProtocol():
    print "Running Interactive Protocol"
    inputParameters = readParameters()
    callInteractiveProtocolScript(inputParameters)
    
#    return template('<b>Hello {{name}}</b>!', name="Name")
    return template('<b>Protocol Executed</b>!')

@route('/runTransferFiles', method='POST')
def runTransferFiles():
    print "Running Transfer Files"
    inputParameters = readParameters()
    #TOBEDONE
    callTransferFilesScript(inputParameters)
    
#    return template('<b>Hello {{name}}</b>!', name="Name")
    return template('<b>Protocol Executed</b>!')


def readParameters():
    print "Reading Parameters..."
    
    inputParameters = {'path_askpass':'scipionPipe'}
    
    
    inputParameters['action'] = request.query.get("action") #InteractiveProtocolMenu || TransferFilesMenu || RunInteractiveProtocol || RunTransferFiles 
    
    inputParameters['machine'] = request.query.get("machine")
    inputParameters['user'] = request.query.get("user")
    inputParameters['commands'] = request.query.getall("commands")
    if inputParameters['action'] == "RunInteractiveProtocol" or inputParameters['action'] == "RunTransferFiles":
        inputParameters['password'] = request.forms.get('password')
    if inputParameters['action'] == "TransferFilesMenu" or inputParameters['action'] == "RunTransferFiles": 
        inputParameters['mode'] = request.query.get("mode","COPY")
        inputParameters['fileList'] = request.query.get("fileList","fileList")
        inputParameters['sourcePath'] = request.query.get("sourcePath")
        inputParameters['targetPath'] = request.query.get("targetPath")
        inputParameters['canRsync'] = request.query.get("canRsync","true")
        inputParameters['availableTunnel'] = request.query.get("availableTunnel","true")
        inputParameters['portTunnel'] = request.query.get("portTunnel",None)
         
    
    print "Parameters: ", inputParameters
    return inputParameters 

    
def getPassword(machine,user):
    print "Getting Password. TBD: Store password locally"
    return ""
    
def callInteractiveProtocolScript(inputParameters):     
    print "Calling Interactive Protocol Script"

    shellCommand="bash " + join(scriptsPath,"managePipe") + inputParameters['password']
    print "shellCommand ", shellCommand
    subprocess.call(shellCommand, shell=True) 
    
    shellCommand="bash " + join(scriptsPath,"sshServerScript") + \
        inputParameters['machine'] + " " + inputParameters['user'] + " " + \
        " \"" + inputParameters['commands'] + "\" " + inputParameters['path_askpass']
    print "shellCommand ", shellCommand
    subprocess.call(shellCommand, shell=True) 

def callTransferFilesScript(inputParameters):     
    print "Calling Transfer Files Script"

    shellCommand="bash " + join(scriptsPath,"managePipe") + inputParameters['password']
    print "shellCommand ", shellCommand
    subprocess.call(shellCommand, shell=True) 
    
    
    #POR AKI
    shellCommand="bash " + join(scriptsPath,"rsyncScript") + inputParameters['commands'] + " " + \
        inputParameters['user'] + " \"" + inputParameters['commands'] + "\""
    print "shellCommand ", shellCommand
    subprocess.call(shellCommand, shell=True)    
  
  
    

run(host='localhost', port=8081)