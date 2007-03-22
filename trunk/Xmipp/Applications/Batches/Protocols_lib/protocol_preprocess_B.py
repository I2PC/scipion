#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script (part B) for Xmipp-based pre-processing of single-particles: 
#  - CTF estimation
#  - extraction of particles
#  - phase flipping
#  - normalization
#  - sort_junk
#
# It is assumed that you have already run part A of this script,
#  and that you have picked the particles for each micrograph
# You also need the xmipp_preprocess_A.py file in the current directory
#
# Example use:
# ./xmipp_preprocess_B.py *.tif
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# All filenames on which to perform preprocessing
ALLMICROGRAPHS='*.tif'
# Project directory
PROJECTDIR="/home2/bioinfo/scheres/work/protocols"
# Downsampling factor (use the same as in part A of this script!)
DOWN=3
# {expert} Log File Directory
LogPath="Logs"
# {expert} Log File Name (.log will be added"
LogFileName="G40P"
#------------------------------------------------------------------------------------------------
# {section} Contrast Transfer Function estimation
#------------------------------------------------------------------------------------------------
# Perform CTF estimation?
DO_CTFESTIMATE=False
# Microscope voltage (in kV)
KV=200
# Spherical aberration
SA=2.26
# Sampling rate after downsampling (in Ang./pix.)
SAM=2.8
# Amplitude Contrast
AMPCONTR=0.1
# Visualize CTF-estimation?
DO_SHOWCTF=True
#------------------------------------------------------------------------------------------------
# {section} Extract particles
#------------------------------------------------------------------------------------------------
# Extract particles from micrographs?
DO_EXTRACT=True
# Name for output selfile
OUT_SELFILE="all_images.sel"
# It is recommend to use the absolute path name
# Directory for extracted images (create if unexisiting)
IMAGESDIR="images"
# (these files are assumed to be in the corresponding subdirectories!)
# Extension for the picked coordinates (Family) file
FAMILY="raw.Common.pos"
# Dimension of the particles to extract (in pix.)
SIZE=64
#------------------------------------------------------------------------------------------------
# {section} Phase correction
#------------------------------------------------------------------------------------------------
# Perform CTF-phase flipping?
DO_FLIP=True
#------------------------------------------------------------------------------------------------
# {section} Normalization
#------------------------------------------------------------------------------------------------
# Perform particle normalization?
DO_NORM=True
# Pixels outside this circle are assumed to be noise and their stddev is set to 1.
# Radius for background circle definition (in pix.)
BG_RAD=30
# Perform ramping background correction?
DO_USE_RAMP=True
# Perform black dust particles removal?
DO_REMOVE_BLACK_DUST=False
# Perform white dust particles removal?
DO_REMOVE_WHITE_DUST=False
#------------------------------------------------------------------------------------------------
# {section} Sort junk
#------------------------------------------------------------------------------------------------
# Perform particle sorting to identify junk particles?
DO_SORT=True
# Visualize sorted particles?
DO_VISUALIZESORT=True
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#  {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class preprocess_B_class:

	#init variables
	def __init__(self,ALLMICROGRAPHS,
                     PROJECTDIR,
                     DOWN,
		     DO_CTFESTIMATE,
                     KV,
                     SA,
                     SAM,
                     AMPCONTR,
                     DO_SHOWCTF,
                     DO_EXTRACT,
                     OUT_SELFILE,
                     IMAGESDIR,
                     FAMILY,
                     SIZE,
                     DO_FLIP,
                     DO_NORM,
                     BG_RAD,
                     DO_USE_RAMP,
                     DO_REMOVE_BLACK_DUST,
                     DO_REMOVE_WHITE_DUST,
                     DO_SORT,
		     LogPath='Logs',
		     LogName='logfile',
		     debug='off',
		     theSmtpLocalServer='mail.cnb.uam.es'
		     ):
	     
            import os
            import sys
            import string
            import log
            
            self.ALLMICROGRAPHS=ALLMICROGRAPHS
            self.PROJECTDIR=PROJECTDIR
            self.DOWN=DOWN
            self.DO_CTFESTIMATE=DO_CTFESTIMATE
            self.KV=KV
            self.SA=SA
            self.SAM=SAM
            self.AMPCONTR=AMPCONTR
            self.DO_SHOWCTF=DO_SHOWCTF
            self.DO_EXTRACT=DO_EXTRACT
            self.OUT_SELFILE=OUT_SELFILE
            self.IMAGESDIR=IMAGESDIR
            self.FAMILY=FAMILY
            self.SIZE=SIZE
            self.DO_FLIP=DO_FLIP
            self.DO_NORM=DO_NORM
            self.BG_RAD=BG_RAD
            self.DO_USE_RAMP=DO_USE_RAMP
            self.DO_REMOVE_BLACK_DUST=DO_REMOVE_BLACK_DUST
            self.DO_REMOVE_WHITE_DUST=DO_REMOVE_WHITE_DUST
            self.DO_SORT=DO_SORT
            self.DO_VISUALIZESORT=DO_VISUALIZESORT
            self.log = log.logpyl(LogPath,LogName,'off')                  
            
            # Execute program
            self.process_all_micrographs()

            
        def perform_ctfestimate(self):
            import os
            iname=self.shortname+'/'+self.downname+'.raw'
            pname=self.shortname+'/'+self.shortname+'_input.param'
            print '*********************************************************************'
            print '*  Estimate CTF for micrograph: '+iname
            paramlist = []
            paramlist.append('image= '+iname+'\n')
            paramlist.append('micrograph_averaging= yes \n')
            paramlist.append('voltage= '+str(self.KV)+'\n')
            paramlist.append('spherical_aberration= '+str(self.SA)+'\n')
            paramlist.append('sampling_rate= '+str(self.SAM)+'\n')
            paramlist.append('particle_horizontal= 128 \n')
	    paramlist.append('Q0= '+str(self.AMPCONTR)+'\n')
            paramlist.append('N_horizontal= 512 \n')
	    paramlist.append('selfile= '+self.downname+'.raw.sel \n')
            paramlist.append('periodogram= yes \n')
    
            fh=open(pname,"w")
            fh.writelines(paramlist)
            fh.close()
            command='xmipp_assign_CTF -i '+pname
            print '* ',command
            self.log.add("command\t:",command)
            os.system(command )
            oname=self.shortname+'/'+self.downname+'_Periodogramavg.ctfmodel'
            self.ctfselfile.append(oname+' 1 \n')

        def visualize_ctfs(self):
            import os
            print '*********************************************************************'
            print '*  Visualizing all CTFs: '
            fh=open("all_ctfs.sel","w")
            fh.writelines(self.ctfselfile)
	    fh.close()
            command='xmipp_show -ctfsel all_ctfs.sel &'
            print '* ',command
            self.log.add("command\t:",command)
            os.system(command )

        def perform_extract_and_phaseflip(self):
            import os
            iname=self.shortname+'/'+self.downname+'.raw' 
            oname=self.downname+'.raw.sel' 
            posname=self.shortname+'/'+self.downname+'.'+FAMILY 
            imgsubdir=self.PROJECTDIR+'/'+self.IMAGESDIR+'/'+self.shortname
            rootname=imgsubdir+'/'+self.shortname+'_'
            logname=self.shortname+'/scissor.log'
            print '*********************************************************************'
            if(os.path.exists(posname)==False):
                print '*  No coordinates-file! Skipping extract particles for micrograph: '+iname
                have_extracted=False
                return have_extracted
            else:
                print '*  Extract partciles from micrograph: '+iname
                have_extracted=True
        
            if (self.DO_FLIP):
                size=2*self.SIZE
            else:
                size=self.SIZE

            if(os.path.exists(imgsubdir)==False):
                self.log.add("command\t:",command)
                os.system('mkdir '+imgsubdir)

            command='xmipp_scissor -i '+iname+' -pos '+posname+' -root '+rootname+' -Xdim '+str(size)+'|grep "corresponding image is set to blank"> '+logname
            print '* ',command
            self.log.add("command\t:",command)
            os.system(command)

        
            fh=open(logname,"r")
            text=fh.read()
	    fh.close()
            num_lines=text.count("\n")
            print '*  Removed '+str(num_lines)+' particles from selfile because they were too near to the border'
            fh=open(oname,"r")
            text = fh.readlines()
	    fh.close()
            fh=open(oname, 'w')
	    fh.writelines([l for l in text if ' -1' not in l])
	    fh.close()

            if (self.DO_FLIP):
                ctfparam=self.shortname+'/'+self.downname+'_Periodogramavg.ctfparam'
                print '*********************************************************************'
                print '*  Flipping phases for images in: '+oname
                command='xmipp_correctphase -i '+oname+' -ctf '+ctfparam
                print '* ',command
                self.log.add("command\t:",command)
                os.system(command)
                command='xmipp_window -i '+oname+' -size '+str(self.SIZE)
                print '* ',command
                self.log.add("command\t:",command)
                os.system(command)
            
            command='cat '+oname+' >> '+self.OUT_SELFILE
            print '* ',command
            self.log.add("command\t:",command)
            os.system(command)


        def perform_normalize(self):
            import os
            iname=self.downname+'.raw.sel'
            print '*********************************************************************'
            print '*  Normalize partciles in: '+iname
            param=' -i ' +iname+' -background circle '+str(self.BG_RAD)
            if (self.DO_USE_RAMP):
                param=param+' -method Ramp'
            if (self.DO_REMOVE_BLACK_DUST):
                param=param+'  -remove_black_dust'
            if (self.DO_REMOVE_WHITE_DUST):
                param=param+'  -remove_white_dust'
            command='xmipp_normalize '+param
            print '* ',command
            self.log.add("command\t:",command)
            os.system(command)


        def perform_sort_junk(self):
            import os
            print '*********************************************************************'
            print '*  Sorting images to identify junk particles in: '+self.OUT_SELFILE
            command='xmipp_sort_junk -i '+self.OUT_SELFILE
            print '* ',command
            self.log.add("command\t:",command)
            os.system(command)
	    if (self.DO_VISUALIZESORT):
		    command='xmipp_show -sel sort_junk.sel & '
		    print '* ',command
		    self.log.add("command\t:",command)
		    os.system(command)


        def process_all_micrographs(self):
            import os
            import glob
            print '*********************************************************************'
            print '*  Pre-processing the following micrographs: '
            for self.filename in glob.glob(self.ALLMICROGRAPHS):
                (self.filepath, self.name) = os.path.split(self.filename)
                print '*  '+self.name

            if(os.path.exists(self.PROJECTDIR+'/'+self.IMAGESDIR)==False):
		    command='mkdir '+self.PROJECTDIR+'/'+self.IMAGESDIR
		    self.log.add("command\t:",command)
		    os.system(command)

	    self.ctfselfile = []
            for self.filename in glob.glob(self.ALLMICROGRAPHS):
		    (self.filepath, self.name) = os.path.split(self.filename)
		    self.shortname=self.name.replace ( '.tif', '' )
		    self.downname='down'+str(self.DOWN)+'_'+self.shortname
		    if(os.path.exists(self.shortname+'/'+self.downname+'.raw')==False):
			    self.downname=self.shortname

		    if (self.DO_CTFESTIMATE):
			    self.perform_ctfestimate()

	    if (self.DO_SHOWCTF):
		    self.visualize_ctfs()

            for self.filename in glob.glob(self.ALLMICROGRAPHS):
                (self.filepath, self.name) = os.path.split(self.filename)

                self.shortname=self.name.replace ( '.tif', '' )
                self.downname='down'+str(self.DOWN)+'_'+self.shortname
                if(os.path.exists(self.shortname+'/'+self.downname+'.raw')==False):
                    self.downname=self.shortname

                if (self.DO_EXTRACT):
                    self.perform_extract_and_phaseflip()
                if (self.DO_NORM):
                    self.perform_normalize()

            if (self.DO_SORT):
                self.perform_sort_junk()

            print '*  Done pre-processing all'
            print '*********************************************************************'

	def close(self):
		self.log.close()

#		
# Main
#     
if __name__ == '__main__':

   	# create preprocess_B_class object

	preprocessB=preprocess_B_class(ALLMICROGRAPHS,
                                       PROJECTDIR,
                                       DOWN,
                                       DO_CTFESTIMATE,
                                       KV,
                                       SA,
                                       SAM,
                                       AMPCONTR,
                                       DO_SHOWCTF,
                                       DO_EXTRACT,
                                       OUT_SELFILE,
                                       IMAGESDIR,
                                       FAMILY,
                                       SIZE,
                                       DO_FLIP,
                                       DO_NORM,
                                       BG_RAD,
                                       DO_USE_RAMP,
                                       DO_REMOVE_BLACK_DUST,
                                       DO_REMOVE_WHITE_DUST,
                                       DO_SORT,
				       LogPath,\
				       LogFileName,\
				       'off',\
				       'mail.cnb.uam.es'\
                                       )

	# close 
	preprocessB.close()

