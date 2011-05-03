#!/usr/bin/env python
#update subveresion
XMIPP_HOME='/home/roberto/xmipp_svn'
XMIPP_TEST=XMIPP_HOME + '/applications/tests'
XMIPP_LOGS=XMIPP_TEST + '/Logs'
import datetime, os
#create log file
d= datetime.datetime.today()
filename = XMIPP_LOGS + '/' + d.strftime("%Y_%m_%d_%H_%M_%S.log") 
print filename
#run svn

def addCmd(title, cmd):
  os.system('echo "%s" >> %s  2>&1' % (title, filename))
  os.system('echo "command: %s" >> %s  2>&1' % (cmd, filename))
  os.system('echo "output:" >> %s  2>&1' % filename)
  os.system('%s >> %s  2>&1' % (cmd, filename)) 

if os.path.exists(filename):
    os.remove(filename) 
os.chdir(XMIPP_HOME)
#Svn update
addCmd("Subversion UPDATE", "svn update")
#configure with test on
confCmd="./scons.configure  MPI_LIBDIR=/usr/lib64/mpi/gcc/openmpi/lib64/ MPI_INCLUDE=/usr/lib64/mpi/gcc/openmpi/include/  MPI_LIB='mpi' gtest=yes java=yes"
addCmd("Scons CONFIGURE", confCmd)
#compile
compCmd="./scons.compile -j 3"
addCmd("Scons COMPILE", compCmd)

