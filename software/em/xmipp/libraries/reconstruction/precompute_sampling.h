/***************************************************************************
 *
 * Authors:    Roberto Marabini 
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef _PRECOMPUTE_SAMPLING_H
#define _PRECOMPUTE_SAMPLING_H

#include <data/sampling.h>
#include <data/symmetries.h>
#include <data/args.h>
#include <fstream>
#include <iostream>
#include <string>

/**@defgroup precomputeSampling precompute_sampling (Precompute samplin points in the projection sphere)
   @ingroup ReconsLibrary */
//@{
/* PDB Phantom Program Parameters ------------------------------------------ */
/** Parameter class for the PDB Phantom program */
class Prog_Sampling_Parameters
{
public:
    /** sampling object */
    Sampling mysampling;

    /** Sampling rate. Distance between sampling points in degrees*/
    double sampling;

    /** A sampling point is neighbor if closer than this value in degrees */
    double neighborhood;


    /** root for output files */
    std::string sampling_file_root;

    /** Symmetry. One of the 17 possible symmetries in
       single particle electron microscopy.
       See details at url
       Possible values are: c1, ci, cs, cn, cnv, cnh, sn,
       dn, dnv, dnh, t, td, th, o, oh, i, ih */
    FileName        fn_sym;

    /** 1-17 symmetry grouips */ 
    int symmetry;

    /** For infinite groups symmetry order*/
    int sym_order;

    /** fil vector with symmetry axis */
    // std::vector <Matrix1D<double> > symmetry_vectors;

    /** maximun tilt angle */
    double max_tilt_angle;
    
    /** minimum tilt angle */
    double min_tilt_angle;


public:
    /** Empty constructor */
    Prog_Sampling_Parameters();

    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, const char **argv);

    /** Usage message.
        This function shows the way of introdustd::cing this parameters. */
    void usage();

    /** Show parameters. */
    void show();

    /** Run. */
    void run();

    /** get all directions related by symmetry to (1,0,0)  */
    void get_sym_vectors(std::vector< Matrix1D<double > > &sym_points);

};
//@}
#endif
