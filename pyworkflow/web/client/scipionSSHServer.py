'''
Created on Sep 30, 2013

@author: adrian
'''
from bottle import route, run, template
import subprocess


@route('/hello/<name>')
def index(name='World'):
    
#    return_code = call("echo Hello World", shell=True)
#    print "return code",return_code
    machineName="pitagoras.cnb.csic.es"
    usernameString="aquintana"
    passwordString=""
    
    subprocess.call("bash sshServerScript", shell=True) 
 
    
    return template('<b>Hello {{name}}</b>!', name=name)

run(host='localhost', port=8081)