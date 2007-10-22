/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#ifndef _PROG_PDBPHANTOM_HH
#define _PROG_PDBPHANTOM_HH

#include "sampling.h"
#include <reconstruction/symmetries.h>
#include <data/args.h>
#include <fstream>
#include <iostream>
#include <string>
#include <data/projection.h>
#include <reconstruction/projection.h>
#include <data/volume.h>
#include <data/funcs.h>
/**@name PDB Phantom program */
//@{
/* PDB Phantom Program Parameters ------------------------------------------ */
/** Parameter class for the PDB Phantom program */
class Prog_Angular_Projection_Matching_Project_Parameters
{
public:
    /** sampling object */
    XmippSampling mysampling;

    /** Sampling rate. Distance between sampling points in degrees*/
    double sampling;

    /** produce projections in psi? */
    double psi_sampling;

    /** maximun tilt angle */
    double max_tilt_angle;
    
    /** minimum tilt angle */
    double min_tilt_angle;
    
    /** root for output files */
    FileName output_file_root;

    /** Symmetry. One of the 17 possible symmetries in
       single particle electron microscopy.
       See details at url
       Possible values are: c1, ci, cs, cn, cnv, cnh, sn,
       dn, dnv, dnh, t, td, th, o, oh, i, ih */
    string symmetry;

    /** For infinite groups symmetry order*/
    int sym_order;
    
    /** Input volume name */
    FileName input_volume;
    
    /**volume to be projecte */
    VolumeXmipp inputVol;
    
    /** projection x dim */
    int Xdim;
    
    /** projection y dim */
    int Ydim;

    /** fil vector with symmetry axis */
    // vector <Matrix1D<double> > symmetry_vectors;
public:
    /** Empty constructor */
    Prog_Angular_Projection_Matching_Project_Parameters();

    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introducing this parameters. */
    void usage();

    /** Show parameters. */
    void show();

    /** Run. */
    void run();

    /** get all directions related by symmetry to (1,0,0)  */
    void get_sym_vectors(vector< Matrix1D<double > > &sym_points);

    /** Project in all the directions between indexes init and end*/
    void project_angle_vector(int my_init, int my_end);
};
//@}
#endif
