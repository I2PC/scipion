#!/usr/bin/env python
#
# Copyright (C) 2014 Diamond Light Source, Karl Levik
#
# 2015-06-22
#
# Usage examples:
# python em_put_movie.py --visit=em11069-2 --sampleid=181182 --microscopeid=1 --mfile=/dls/m01/data/2014/em11069-2/raw/May11_10.01.12.mrc
#
# Assuming the previous returns --movieid=2:
# python em_put_movie.py --movieid=2 --micrograph=/dls/m01/data/2014/em11069-2/.ispyb/pngs/May10_13.08.38_SumCorr.png --powerspectrum1=/dls/m01/data/2014/em11069-2/.ispyb/pngs/May10_13.08.38_RawFFT.png --powerspectrum2=/dls/m01/data/2014/em11069-2/.ispyb/pngs/May10_13.08.38_CorrFFT.png --drift=/dls/m01/data/2014/em11069-2/.ispyb/May10_13.08.38_Log.dat --nimages=22 --framelen=0.2 --totexp=201.7 --magnification=900010 --samplepixsize=20.4 --doseperframe=301.2 --totdose=52788.6 --runstatus='Successful' --comments='Lucky one' --rundir=/dls/m01/data/2014/em11069-2/raw/ --binning=1 --particlediameter=23.9 --pixelsize=5.5 --boxsize=43.7 --minresol=2.0 --maxresol=3.2 --mindefocus=7.0 --maxdefocus=11.0 --defocusstepsize=0.5 --astigmatism=0.042 --extractsize=202.1 --bgradius=320.8 --stime='1998-06-23 14:00:01' --etime='1998-06-23 14:00:59'

#

import string
import logging
from logging.handlers import RotatingFileHandler
import time
import os
import sys
import datetime

from ispyb_api.dbconnection import dbconnection
from ispyb_api.core import core
from ispyb_api.mxacquisition import mxacquisition

from datetime import datetime

if __name__ == '__main__' :
    
    def exit(code, message=None):
        dbconnection.disconnect()
        if not message is None:
            print(message)
        sys.exit(code)
    
    logging.info("test")
    
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--movieid", dest="movieid", help="Id for movie", metavar="INTEGER")
    parser.add_option("--mfile", dest="filename", help="Path to movie file", metavar="FILE")
    parser.add_option("--visit", dest="visit", help="Visit name", metavar="STRING")
    parser.add_option("--groupid", dest="groupid", help="Id for group", metavar="INTEGER")
    parser.add_option("--sampleid", dest="sampleid", help="Id for sample", metavar="INTEGER")
    parser.add_option("--microscopeid", dest="microscopeid", help="Id for microscope", metavar="INTEGER")
    # parser.add_option("--micrograph", dest="micrograph", help="Path to micrograph image", metavar="FILE")
    # parser.add_option("--powerspectrum1", dest="powerspectrum1", help="Path to 1st power spectrum image", metavar="FILE")
    # parser.add_option("--powerspectrum2", dest="powerspectrum2", help="Path to 2nd power spectrum image", metavar="FILE")
    # parser.add_option("--drift", dest="drift", help="Path to drift .dat file", metavar="FILE")
    # parser.add_option("--nimages", dest="n_images", help="Number of images in .mrc file", metavar="INTEGER")
    # #parser.add_option("--framelen", dest="frame_len", help="Frame length (s)", metavar="FLOAT")
    # parser.add_option("--totexp", dest="tot_exp", help="Total exposure (s)", metavar="FLOAT")
    # parser.add_option("--magnification", dest="magnification", help="(X)", metavar="FLOAT")
    # #parser.add_option("--samplepixsize", dest="sample_pix_size", help="Sample pix size (A/pix)", metavar="FLOAT")
    # #parser.add_option("--doseperframe", dest="dose_per_frame", help="Dose per frame (e-/A^2)", metavar="FLOAT")
    # parser.add_option("--totdose", dest="tot_dose", help="Total dose (e-/A^2)", metavar="FLOAT")
    # parser.add_option("--runstatus", dest="run_status", help="Status of movie", metavar="STRING")
    # parser.add_option("--comments", dest="comments", help="User comments", metavar="STRING")
    # #parser.add_option("--rundir", dest="run_dir", help="Directory for raw data collected", metavar="STRING")
    # parser.add_option("--binning", dest="binning", help="Binning: 1 or 2 (X)", metavar="INTEGER")
    # parser.add_option("--particlediameter", dest="particle_diameter", help="Particle diameter (nm)", metavar="FLOAT")
    # #parser.add_option("--pixelsize", dest="pixel_size", help="Pixel size (nm)", metavar="FLOAT")
    # parser.add_option("--boxsize", dest="box_size", help="Box size (nm)", metavar="FLOAT")
    # parser.add_option("--minresol", dest="min_resol", help="Min resolution (A)", metavar="FLOAT")
    # parser.add_option("--maxresol", dest="max_resol", help="Max resolution (A)", metavar="FLOAT")
    # parser.add_option("--mindefocus", dest="min_defocus", help="Min defocus (A)", metavar="FLOAT")
    # parser.add_option("--maxdefocus", dest="max_defocus", help="Max defocus (A)", metavar="FLOAT")
    # parser.add_option("--defocusstepsize", dest="defocus_step_size", help="Defocus step size (A)", metavar="FLOAT")
    # parser.add_option("--astigmatism", dest="astigmatism", help="Amount of astigmatism (A)", metavar="FLOAT")
    # parser.add_option("--extractsize", dest="extract_size", help="Extract size (nm)", metavar="FLOAT")
    # parser.add_option("--bgradius", dest="bg_radius", help="Background radius (nm)", metavar="FLOAT")
    # parser.add_option("--stime", dest="stime", help="Start time (yyyy-mm-dd hh24:mi:ss)", metavar="TIME")
    # parser.add_option("--etime", dest="etime", help="End time (yyyy-mm-dd hh24:mi:ss)", metavar="TIME")
    parser.add_option("--db", dest="db", help="Database to use: dev, test or prod (default)", metavar="STRING")

    (opts, args) = parser.parse_args()

    cursor = None
    if opts.db is None or opts.db == "prod": 
        cursor = dbconnection.connect_to_prod()
    elif opts.db == "dev":
        cursor = dbconnection.connect_to_dev()
    elif opts.db == "test":
        cursor = dbconnection.connect_to_test()
    else:
        exit(1, "ERROR: Invalid database")

    # Find the id for the visit
    visitid = core.retrieve_visit_id(cursor, opts.visit)
    
    result = None
    if visitid is not None or opts.movieid is not None:
	# EMMovieId, blsessionId,blsampleId,movieFile,p_NOImages,frameLength,totalExposure,magnification,samplePixSize,dosePerFrame,totalDose,runStatus,comments,runDirectory,binning,particleDiameter,pixelSize,boxSize,minResolution,maxResolution, minDefocus, maxDefocus, defocusStepSize, amountAstigmatism, extractSize, starttime, endtime
        
        # Store a movie data collection ...
        params = mxacquisition.get_data_collection_params()
        params['parentid'] = opts.groupid
        params['visitid'] = visitid
        params['sampleid'] = opts.sampleid
        params['detectorid'] = opts.microscopeid
        
        if not opts.filename is None:
            imgdir = opts.filename[:opts.filename.rfind('/')]
        
            params['imgdir'] = imgdir
            params['imgprefix'] = opts.filename[opts.filename.rfind('/') + 1 : opts.filename.rfind('_')]
            params['imgsuffix'] = opts.filename[opts.filename.rfind('.')+1:]
            params['file_template'] = opts.filename

        # TODO: finish the mapping here ...
        # params['xtal_snapshot1'] = opts.powerspectrum2
        # params['xtal_snapshot2'] = opts.powerspectrum1
        # params['dat_file'] = opts.drift # not sure where drift should go yet
        # params['xtal_snapshot4'] = opts.micrograph
        #
        # start_time = None
        # end_time = None
        # if not opts.stime is None:
        #    params['starttime'] = datetime.datetime.strptime(opts.stime, '%Y-%m-%d %H:%M:%S')
        # if not opts.etime is None:
        #    params['endtime'] = datetime.datetime.strptime(opts.etime, '%Y-%m-%d %H:%M:%S')
        #
        # params['n_images'] = opts.n_images
        # params['exp_time'] = opts.tot_exp
        # params['magnification'] = opts.magnification
        # params['total_absorbed_dose'] = opts.tot_dose
        # params['run_status'] = opts.run_status
        # params['comments'] = opts.comments
        # params['binning'] = opts.binning
        # params['particle_diameter'] = opts.particle_diameter
        # params['box_size_ctf'] = opts.box_size
        # params['min_resolution'] = opts.min_resol
        # params['resolution'] = opts.max_resol
        # params['min_defocus'] = opts.min_defocus
        # params['max_defocus'] = opts.max_defocus
        # params['defocus_step_size'] = opts.defocus_step_size
        # params['amount_astigmatism'] = opts.astigmatism
        # params['extract_size'] = opts.extract_size
        # params['bg_radius'] = opts.bg_radius

        dc_id = None        
        if opts.movieid is None:
            dc_id = mxacquisition.insert_data_collection(cursor, params.values())
        else:
            dc_id = mxacquisition.update_data_collection(cursor, params.values())

        exit(0, "--movieid=%d" % dc_id)
    else:
        exit(1, "ERROR: Neither visit nor movieid found/provided.")
 
    


    
