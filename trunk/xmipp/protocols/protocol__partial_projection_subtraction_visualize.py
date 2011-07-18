#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_projmatch.py
#
# Example use:
# python visualize_projmatch.py
#
# This script requires that protocol_projmatch.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
#show results for iterations
""" Nothing here yet
"""
DoNotFillThis=1

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#


class visualize_partial_projection_subtraction:

    #init variables
    def __init__(self,
                _ProtocolName
                ):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,logging,arg
        import visualization

        # Import the corresponding protocol, get WorkingDir and go there
        pardir=os.path.abspath(os.getcwd())
	#shutil.copy(_ProtocolName,'v_partial_protocol.py')

        #import v_partial_protocol
	xmpi_run_file      = 'readDocfileAndPairExperimentalAndReferenceImages.sh'
	WorkingDir='ProjMatch/xmipp_2.4_subtraction_crunchy'
        _Iteration_Working_Directory='/Iter_'+  str(6)
	os.chdir(WorkingDir+'/'+_Iteration_Working_Directory)

        #command = "cat " + xmpi_run_file + " | awk '{print $3 \" 1\\n\" $5 \" 1\\n\" $7 \" 1\"}'"
        fn=open('kk.sel','w')
	for y in open(xmpi_run_file).readlines():
            x =y.split(' ')
	    #s=  "%s 1\n%s 1\n%s 1\n"%(x[4],x[6],x[8].strip('\n'))
            s=  "%s 1\n%s 1\n%s 1\n"%(x[12],x[14],x[16].strip('\n'))
	    s = s.replace(';','')
            fn.write(s)
	fn.close()


    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'


#		
# Main
#      
if __name__ == '__main__':

    import sys
    ProtocolName=""#sys.argv[1]

    # create projmatch_class object
    visualize_partial_projection=visualize_partial_projection_subtraction(
                                                  ProtocolName)
    # close 
    visualize_partial_projection.close()

