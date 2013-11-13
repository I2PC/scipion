'''
Created on Sep 30, 2013

@author: adrian
'''
from bottle import route, run, template, request
import subprocess
from os.path import join

#from transferFiles import openTransferFilesMenu


@route('/openScipionLocalMenu')
def openScipionLocalMenu():
    print "Opening Scipion Local Menu"
    action = request.query.get("action") #openInteractiveProtocolMenu || openTransferFilesMenu || runInteractiveProtocol || runTransferFiles

    if action == 'openInteractiveProtocolMenu':
        from interactiveProtocol import openInteractiveProtocolMenu
        return openInteractiveProtocolMenu()
#    elif action == 'openTransferFilesMenu':
#        return openTransferFilesMenu ()

#    m = __import__ ('foo')
#    func = getattr(m,'bar')
#    func()
    

    


@route('/runAction', method='POST')
def runAction():
    print "Running Action..."
    
    action = request.forms.get("action") #openInteractiveProtocolMenu || openTransferFilesMenu || runInteractiveProtocol || runTransferFiles
    if action == 'runInteractiveProtocol':
        from interactiveProtocol import runInteractiveProtocol
        return runInteractiveProtocol()
    

    

run(host='localhost', port=8081)