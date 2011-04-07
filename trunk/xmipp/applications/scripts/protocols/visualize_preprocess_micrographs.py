#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_preprocess_micrographs.py
#
# Example use:
# ./visualize_preprocess_micrographs.py
#
# This script requires that protocol_preprocess_micrographs.py is in the current directory
#
# Author: Carlos Oscar Sorzano, April 2011
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
import os,sys
class visualize_micrographs_class:
    #init variables
    def __init__(self,WorkingDir):
        summaryFile=WorkingDir+"/micrographs.sel"
        if os.path.exists(summaryFile):
            os.system("xmipp_visualize_preprocessing_micrographj -i "+summaryFile+" &")
        else:
            import tkMessageBox
            message="There is no result yet"
            tkMessageBox.showerror("Error", message)

if __name__ == '__main__':
    WorkingDir=sys.argv[1]
    visualize=visualize_micrographs_class(WorkingDir)

