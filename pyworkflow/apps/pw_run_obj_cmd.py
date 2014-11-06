import sys
from pyworkflow.em.showj import *

if __name__ == '__main__':
    #TODO: REMOVE THIS AFTER DEBUGGING
    print "ARGS: "
    for i, arg in enumerate(sys.argv):
        print "%02d: %s" % (i, arg)

    cmd = sys.argv[1]
    if cmd == OBJ_CMDA:
        print 'doing cmd ' + cmd