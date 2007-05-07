#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_rotspectra.py
#
# Example use:
# python visualize_rotspectra.py
#
# This script requires that protocol_rotspectra.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize aligned images
DoShowAlignedImages=False
# Visualize average image
DoShowAverage=False
# Visualize spectra
DoShowSpectra=False
# Visualize SOM made with spectra data
DoShowSOMSpectra=False
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_rotspectra_class:

    #init variables
    def __init__(self,_DoShowAlignedImages,
                      _DoShowAverage,
                      _DoShowSpectra,
                      _DoShowSOMSpectra,
                      _WorkingDir,
                      ):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self._WorkDirectory=_WorkingDir
        os.chdir(self._WorkDirectory)
        # import corresponding protocol
        sys.path.append(self._WorkDirectory)
        import protocol_rotspectra_backup

        _LogDir=protocol_rotspectra_backup.LogDir
        _ProjectDir=protocol_rotspectra_backup.ProjectDir
        self._SelFileName=_ProjectDir+'/'+str(protocol_rotspectra_backup.SelFileName)
        self._SomName=protocol_rotspectra_backup.SomName
        self._SpectraName=protocol_rotspectra_backup.SpectraName
        self.mylog=log.init_log_system(_ProjectDir,
                                       _LogDir,
                                       sys.argv[0],
                                       self._WorkDirectory)

        if (_DoShowAlignedImages):
            self.visualize_AlignedImages(self._SelFileName)
        if (_DoShowAverage):
            self.visualize_Average(self._SelFileName)
        if (_DoShowSpectra):
            self.visualize_Spectra(self._SpectraName)
        if (_DoShowSOMSpectra):
            self.visualize_SOMSpectra(self._SomName,self._SpectraName)
            
        # Return to parent dir
        os.chdir(os.pardir)

    def visualize_AlignedImages(self,_SelFileName):
         import os
         command='xmipp_show -sel '+ os.path.basename(_SelFileName) + '&'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)


    def visualize_Average(self,_SelFileName):
         import os
         selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
         command='xmipp_show -img '+ selfile_without_ext + '.med.xmp' + '&'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)

    def visualize_Spectra(self,_SpectraName):
         import os
         selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
         command='xmipp_show'+ \
                  ' -spect ' + str(_SpectraName) + '&'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)

    def visualize_SOMSpectra(self,_SomName,_SpectraName):
         import os
         selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
         command='xmipp_show -spectsom ' + \
                  str(_SomName) + \
                  ' -din ' + str(_SpectraName) + '&'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
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
    WorkingDir=sys.argv[1]

    # create rotspectra_class object
    visualize_rotspectra=visualize_rotspectra_class(DoShowAlignedImages,
                                                    DoShowAverage,
                                                    DoShowSpectra,
                                                    DoShowSOMSpectra,
                                                    WorkingDir)
    # close 
    visualize_rotspectra.close()

