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
import argparse

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
    
    parser = argparse.ArgumentParser()
    add = parser.add_argument

    add("--movieid", help="Id for movie", type=int)
    add("--mfile", help="Path to movie file")
    add("--visit", help="Visit name")
    add("--parentid", help="Id for group", type=int)
    add("--sampleid", help="Id for sample", type=int)
    add("--detectorid", help="Id for microscope", type=int)
    # add("--micrograph", help="Path to micrograph image")
    # add("--powerspectrum1", help="Path to 1st power spectrum image")
    # add("--powerspectrum2", help="Path to 2nd power spectrum image")
    # add("--drift", help="Path to drift .dat file")
    # add("--nimages", help="Number of images in .mrc file", type=int)
    # #add("--framelen", help="Frame length (s)", type=float)
    # add("--totexp", help="Total exposure (s)", type=float)
    # add("--magnification", help="(X)", type=float)
    # #add("--samplepixsize", help="Sample pix size (A/pix)", type=float)
    # #add("--doseperframe", help="Dose per frame (e-/A^2)", type=float)
    # add("--totdose", help="Total dose (e-/A^2)", type=float)
    # add("--runstatus", help="Status of movie")
    # add("--comments", help="User comments")
    # #add("--rundir", help="Directory for raw data collected")
    # add("--binning", help="Binning: 1 or 2 (X)", type=int)
    # add("--particlediameter", help="Particle diameter (nm)", type=float)
    # #add("--pixelsize", help="Pixel size (nm)", type=float)
    # add("--boxsize", help="Box size (nm)", type=float)
    # add("--minresol", help="Min resolution (A)", type=float)
    # add("--maxresol", help="Max resolution (A)", type=float)
    # add("--mindefocus", help="Min defocus (A)", type=float)
    # add("--maxdefocus", help="Max defocus (A)", type=float)
    # add("--defocusstepsize", help="Defocus step size (A)", type=float)
    # add("--astigmatism", help="Amount of astigmatism (A)", type=float)
    # add("--extractsize", help="Extract size (nm)", type=float)
    # add("--bgradius", help="Background radius (nm)", type=float)
    # add("--stime", help="Start time (yyyy-mm-dd hh24:mi:ss)")
    # add("--etime", help="End time (yyyy-mm-dd hh24:mi:ss)")
    add("--db", help="Database to use: dev, test or prod (default)")

    args = parser.parse_args()

    try:
        connect_func = getattr(dbconnection, 'connect_to_' + args.db)
        cursor = connect_func()
    except Exception as ex:
        exit(1, "ERROR: Could not connect to database")

    # Find the id for the visit
    visitid = core.retrieve_visit_id(cursor, args.visit)
    
    result = None
    if visitid is not None or args.movieid is not None:

        # Store a movie data collection ...
        params = mxacquisition.get_data_collection_params()

        # Update the data collection dictionary with all arguments
        # passed from the command line
        for k, v in vars(args).iteritems():
            if k in params:
                params[k] = v

        params['visitid'] = visitid

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
        #    params['starttime'] = opts.stime
        # if not opts.etime is None:
        #    params['endtime'] = opts.etime
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
 
    


    
