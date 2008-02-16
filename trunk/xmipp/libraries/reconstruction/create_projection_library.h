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
class Prog_create_projection_library_Parameters
{
public:
    /** sampling object */
    XmippSampling mysampling;

    /** Sampling rate. Distance between sampling points in degrees*/
    double sampling;

    /**Perturb angles of reference projections*/
    double perturb_projection_vector;

    /** produce projections in psi? */
    double psi_sampling;

    /** maximun tilt angle */
    double max_tilt_angle;
    
    /** minimum tilt angle */
    double min_tilt_angle;
    
    /** root for output files */
    FileName output_file_root;

    /** Quiet */
    bool quiet;

    /** docfile with experimental images angles. 
        This information is used to generate only angles close to experimental
        images */
    DocFile DFexperimentalImages;

    /** vector with valid proyection directions after looking for 
        directions close to experimental data. 
         */
    std::vector <Matrix1D<double> > close_points_angles;

    /** filename with experimental images angles. 
        This information is used to generate only angles close to experimental
        images */
    FileName FnexperimentalImages;
    
    /** enabled angular_distance */
    bool angular_distance_bool;
   
    /** enabled angular_distance */
    double angular_distance;

   /** Symmetry. One of the 17 possible symmetries in
       single particle electron microscopy.
       See details at url
       Possible values are: c1, ci, cs, cn, cnv, cnh, sn,
       dn, dnv, dnh, t, td, th, o, oh, i, ih */
    std::string symmetry;

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
    // std::vector <Matrix1D<double> > symmetry_vectors;
public:
    /** Empty constructor */
    Prog_create_projection_library_Parameters();

    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introdustd::cing this parameters. */
    void usage();

    /** Show parameters. */
    void show();

    /** Run. */
    void run();

    /** get all directions related by symmetry to (1,0,0)  */
    void get_sym_vectors(std::vector< Matrix1D<double > > &sym_points);

    /** Project in all the directions between indexes init and end*/
    void project_angle_vector(int my_init, int my_end, bool verbose = true);

    /** Remove projection points no closer to experimetal data than*/
    void remove_points_not_close_to_experimental_points(void);
};
//@}
#endif
