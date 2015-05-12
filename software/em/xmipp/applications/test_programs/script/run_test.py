#!/usr/bin/env python
#update subveresion
from  config import XMIPP_BASE
from  config import XMIPP_HOME
from  config import XMIPP_TEST
from  config import XMIPP_LOGS
from  config import XMIPP_BRANCH
from  config import XMIPP_MPI_LIBDIR
from  config import XMIPP_JAVA_HOME
from  config import XMIPP_JNI_CPPPATH
from  config import XMIPP_MPI_INCLUDE

import datetime, os,shutil
#create log file
d= datetime.datetime.today()
filename = XMIPP_LOGS + '/' + d.strftime("%d_%H.log") 
print "Logging file:", filename
#run svn

def addCmd(title, cmd):
    os.system('echo "%s" >> %s  2>&1' % (title, filename))
    os.system('echo "command: %s" >> %s  2>&1' % (cmd, filename))
    os.system('echo "output:" >> %s  2>&1' % filename)
    os.system('%s >> %s  2>&1' % (cmd, filename)) 

if not os.path.exists(XMIPP_LOGS):
    os.mkdir(XMIPP_LOGS) 
if  os.path.exists(XMIPP_HOME):
    shutil.rmtree(XMIPP_HOME) 
   
#Svn update
os.chdir(XMIPP_BASE)
addCmd("Git CLONE", "git clone http://git.code.sf.net/p/newxmipp/code "+XMIPP_HOME)
#configure with test on
os.chdir(XMIPP_HOME)
addCmd("Git CHECKOUT", "git checkout "+XMIPP_BRANCH)
confCmd="./install.sh gui=false -j 8"
addCmd("Configuration compilation and installation", confCmd)
#confCmd="./setup.py configure MPI_LIBDIR="+XMIPP_MPI_LIBDIR+" JAVA_HOME="+XMIPP_JAVA_HOME+" JNI_CPPPATH="+str(XMIPP_JNI_CPPPATH)+" MPI_INCLUDE="+XMIPP_MPI_INCLUDE+" gtest=yes java=yes -j 8"
#confCmd="./scons.configure MPI_LIBDIR="+XMIPP_MPI_LIBDIR+" JAVA_HOME="+XMIPP_JAVA_HOME+" JNI_CPPPATH="+str(XMIPP_JNI_CPPPATH)+" MPI_INCLUDE=/usr/lib64/mpi/gcc/openmpi/include/  MPI_LIB='mpi' gtest=yes java=yes"
#addCmd("Scons CONFIGURE", confCmd)
#In Theory I do not need to compile twice but I keep getting text file busy messages ROB
#compCmd="./setup.py compile -j 8"
#addCmd("Scons COMPILE", compCmd)
#compCmd="./setup.py run_tests -j 8"
#addCmd("Scons Run Test", compCmd)

import parse_tests_results
#parse xml files
failed = False
message = ""
subject = "XMIPP compilation "
results = parse_tests_results.parseResults()
for group in results:
    n = len(group.tests)
    if group.failed > 0:
        message += "   [  FAILED  ] Group: %s, failed %d of %d tests.\n" % (group.name, group.failed, n)
        failed = True
    else:
        message += "   [  OK      ] Group: %s, succeed %d tests.\n" % (group.name, n) 
        
if failed:
    subject += "FAILED"
else:
    subject += "OK"
subject += "(" + XMIPP_BRANCH + ")"

#print subject
#print message
#Send notification mail
import mail
from config import toaddrs
from config import fromaddr 
mail.mail(toaddrs, fromaddr, subject, message, filename)
    


