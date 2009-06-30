#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_CL2D.py
#
# Example use:
# python visualize_cl2d.py protocol_cl2d.py
#
# This script requires that protocol_cl2d.py is in the current directory
#
# Author: Carlos Oscar Sánchez Sorzano
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize the class averages of all iterations in matrix-view?
DoShowAllIter=True
# Separately visualize class averages of the last iteration?
DoShowLastIter=False
# Open files with the number of images assigned?
DoShowFractions=True

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_CL2D_class:

    #init variables
    def __init__(self,
                 DoShowAllIter,
                 DoShowLastIter,
                 DoShowFractions,
                 ProtocolName):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log

        # Import the corresponding protocol and get WorkingDir
        shutil.copy(ProtocolName,'protocol.py')
        import protocol
        self.WorkingDir=protocol.WorkingDir
        self.DoShowFractions=DoShowFractions

        if (DoShowAllIter):
            self.show_all_iter()
        if (DoShowLastIter):
            self.show_last_iter()

        # Return to parent dir and remove protocol.py(c)
        if (os.path.exists('protocol.py')):
            os.remove('protocol.py')
        if (os.path.exists('protocol.pyc')):
            os.remove('protocol.pyc')
        
        # Finish
        self.close()

    def show_all_iter(self):
        import os,glob
        import selfile
        
        selfiles=glob.glob(self.WorkingDir+'/class_level_??_.sel')
        selfiles.sort()
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            for selfile in selfiles:
                command='xmipp_show -sel '+selfile+ ' &'
                print '* ',command
                os.system(command)
                if (self.DoShowFractions):
                    docfile=selfile.replace('.sel','.doc')
                    command='xmipp_edit -i '+docfile+ ' &'
                    print '* ',command
                    os.system(command)                    

    def show_last_iter(self):
        # Visualize class averages:
        import os,glob
        selfiles=glob.glob(self.WorkingDir+'/class_level_??_.sel')
        selfiles.sort()
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            # print logfile to screen:
            lastselfile=selfiles[-1]
            command='xmipp_show -sel '+lastselfile+ ' &'
            print '* ',command
            os.system(command)
            if (self.DoShowFractions):
                docfile=lastselfile.replace('.sel','.doc')
                command='xmipp_edit -i '+docfile+ ' &'
                print '* ',command
                os.system(command)                    
                          
    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    import sys
    ProtocolName=sys.argv[1]

    # create ML2D_class object
    visualize_CL2D=visualize_CL2D_class(DoShowAllIter,
                                        DoShowLastIter,
                                        DoShowFractions,
                                        ProtocolName)
