#!/usr/bin/env python
#update subveresion
from  config import XMIPP_BASE
from  config import XMIPP_HOME
from  config import XMIPP_TEST
from  config import XMIPP_LOGS

import datetime, os,shutil
#create log file
d= datetime.datetime.today()
filename = XMIPP_LOGS + '/' + d.strftime("%Y_%m_%d_%H_%M_%S.log") 
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
addCmd("Subversion CHECKOUT", "svn co http://newxmipp.svn.sourceforge.net/svnroot/newxmipp/trunk/xmipp")
#configure with test on
os.chdir(XMIPP_HOME)
confCmd="./scons.configure  QTDIR=/usr/share/qt3 MPI_LIBDIR=/usr/lib64/mpi/gcc/openmpi/lib64/ MPI_INCLUDE=/usr/lib64/mpi/gcc/openmpi/include/  MPI_LIB='mpi' gtest=yes java=yes"
addCmd("Scons CONFIGURE", confCmd)
#In Theory I do not need to compile twice but I keep gettiong text file fusy messages ROB
compCmd="./scons.compile -j 3"
addCmd("Scons COMPILE", compCmd)
compCmd="./scons.compile -j 3 run_tests"
addCmd("Scons Run Test", compCmd)
import parse_test
#parse xml files
parse_test.main(filename)
