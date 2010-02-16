/***************************************************************************
 *
 * Authors:     Slavica JONIC (slavica.jonic@impmc.jussieu.fr, slavica.jonic@a3.epfl.ch)
 *              Jean-Noël PIOCHE (jnp95@hotmail.com) 
 *        
 * Biomedical Imaging Group, EPFL (Lausanne, Suisse).
 * Structures des Assemblages Macromoléculaires, IMPMC UMR 7590 (Paris, France).
 * IUT de Reims-Châlons-Charleville (Reims, France).
 *
 * Last modifications by JNP the 27/05/2009 15:52:45  
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

#ifndef __PROJECTION_REAL_SHEARS_H__
#define __PROJECTION_REAL_SHEARS_H__

/***************************************************************************** 
 *    System includes 
 ****************************************************************************/ 
#include <stdlib.h> 
#include <stdio.h>
#include <stddef.h>
#include <math.h> 
#include <float.h> 
#include <string.h>
#include <iostream>
 
/***************************************************************************** 
 *    Toolbox defines 
 ****************************************************************************/ 
#include <external/bilib/configs.h> 
#include <external/bilib/headers/messagedisplay.h> 
#include <data/selfile.h> 
 
/***************************************************************************** 
 *    Toolbox includes 
 ****************************************************************************/ 
#include <external/bilib/headers/linearalgebra.h>
#include <external/bilib/headers/getput.h>
#include <external/bilib/headers/getputd.h>
#include <external/bilib/headers/changebasis.h>
 
/***************************************************************************** 
 *    New toolbox includes 
 ****************************************************************************/
#include <data/docfile.h>
#include <data/selfile.h>
#include <data/projection.h>
#include <data/volume.h>


#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503130808472633361816000000000e-16
#endif

/**@name ProjRealShears Projection library program using Real-Shears */
//@{
/// Structure for holding a volume
typedef struct 
{ 
    double *Volume; 
    long   nx_Volume, ny_Volume, nz_Volume; 
    short  *Proj_dims; 
    double *Identity_orientN; 
    double *Identity_orientV; 
    double *Identity_orientW; 
    double *IdentityOrigin; 
    double *K123; 
    double *Lambda123; 
    double *Gama123; 
    double *InitPsiThetaPhi; // Angles 
    double *InitDelta123; //Shifts 
    double *Output; 
    double PeriodOfSamplingInVDirection; 
    double PeriodOfSamplingInWDirection; 
} VolumeStruct;

/// Prepare a volume to be projected
void prepareStructVolume(const Matrix3D<double> &V, VolumeStruct &Data);

/// Allocates and fixes some VolumeStruct fields 
void allocAndInit_VolumeStruct(VolumeStruct &Data2);

///Main class of this program
class Projection_real_shears
{
    //--------------- Fields ---------------
    public :
        /// Filename of the projection parameters file (input).
        FileName fn_proj_param;

        /// Selection file with all projections (output).
        FileName fn_sel_file;

        ///Tell if the displaying is active
        bool display;

        ///X Translation parameter
        double shiftX;
        ///Y Translation parameter
        double shiftY;

        ///Projection Xdim
        int proj_Xdim;
        ///Projection Ydim
        int proj_Ydim;
        ///Projection Zdim
        int proj_Zdim;

        ///Volume file
        FileName fnPhantom;

        ///Seed name of the projections files
        std::string fnProjectionSeed;

        ///Starting number for the projections files
        int starting;
        ///Current number of the projection file to save
        int num_file;
        ///Extension name of the projections files
        std::string fn_projection_extension;
        ///Angular file
        FileName fn_angle;
        
        ///Projection to save
        Projection proj;
        ///Selection file to fill and save
        SelFile SF;
        ///Content of the angle file
        DocFile DF;
        ///Basics Data for projections
        VolumeStruct Data;

    //------------------ Functions -------------------
    public :
        ///Main function
        int ROUT_project_real_shears();

        ///Does starting instructions
        int start_to_process();

        ///Reads a DocLine and fills Data fields
        int read_a_DocLine();

        ///Executes instructions for one projection
        int do_oneProjection(VolumeStruct &Data2);

        ///Writes the projection file obtained
        int write_projection_file(int numFile);

        ///Does finish instructions
        int finish_to_process();

        ///Read input and output file parameters only
        void read(int argc, char **argv);

        /// Usage message. This function shows the way of introducing this parameters.
        void usage();

        ///"Overloaded" function in order to use translation parameters
        void read(const FileName &fn_proj_param); //read_withShift [...]

        ///Desallocates VolumeStruct fields
        void del_VolumeStruct(VolumeStruct &Data2);

        ///Transforms angles from (Ry, Rz, Ry) to (Rx, Ry, Rz) system. Returns possible error.
        int angles_transcription(double *angles, double *Lambda123);

        ///Returns a pointer of the multiplication of the 5 matrices.
        double *MatrixBem(double phi, double theta, double psi, double x0, double y0, double scale_x, double scale_y, double scale_z);

        ///Main compute function.
        int ROUT_project_execute(VolumeStruct &Data2);

        ///Computes projection
        int Compute_projection(
            double *Parameters, 
            double *Coef_x, 
            double *Coef_y, 
            double *Coef_z, 
            long    Nx, 
            long    Ny, 
            long    Nz, 
            short  *Proj_dims, 
            double *Identity_orientN, 
            double *Identity_orientV, 
            double *Identity_orientW, 
            double *IdentityOrigin, 
            double *PeriodOfSamplingInVDirection, 
            double *PeriodOfSamplingInWDirection, 
            double *RightOperHlp, 
            double *Ac, 
            double *Projection, 
            double *B);
        ///Computes one iteration of the projection
        int do_compute_projection(
            double *VWOrigin, 
            long    Ndv, 
            long    Ndw, 
            double *Identity_orientV, 
            double *Identity_orientW, 
            double  dv, 
            double  dw, 
            double *CoefVolume, 
            double  absscale, 
            double *Binv, 
            double *BinvCscaled, 
            long    ksimax, 
            int    *arr, 
            long    CoefVolumeNx, 
            long    CoefVolumeNy, 
            long    lmax, 
            long    mmax, 
            double *Projection);
};
//@}
#endif
