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
#ifndef _PROJECTION_HH
#  define _PROJECTION_HH

#include <XmippData/xmippImages.hh>
#include <XmippData/xmippProjection.hh>
#include "grids.hh"

/*---------------------------------------------------------------------------*/
/* PROJECTION                                                                */
/*---------------------------------------------------------------------------*/
/**@name Projections */
//@{
// Projecting functions ====================================================
#define FORWARD  1
#define BACKWARD 0
/**@name Single particle projections */
//@{
#define ARTK     1
#define CAVK     2
#define COUNT_EQ 3
#define CAV      4

/** From blob volumes (Forward & Backward(ART)).
    Project a grid volume with a blob basis.
    The Grid volume is projected onto a projection plane defined by
    (rot, tilt, psi) (1st, 2nd and 3rd Euler angles). The projection
    is previously is resized to Ydim x Xdim and initialized to 0.
    The projection itself, from now on, will keep the Euler angles.
    
    FORWARD process:
       Each volume of the grid is projected on to the projection plane.
       The output is the projection itself and a normalising image, the
       normalising image is the projection of the same grid supposing
       that all blobs are of value 1. This normalising image is used by
       the ART process
    
    BACKWARD process:
       During the backward process the normalising projection contains
       the correction image to apply to the volume (in the ART sense).
       The output is the volume itself, the projection image is useless
       in this case, and the normalising projection is not modified at
       all.

    As for the mode, valid modes are ARTK, CAVK, COUNT_EQ.
    
    M is the matrix corresponding to the projection process.
    */
template <class T>
void project_Volume(GridVolumeT<T> &vol,
   const ImageOver &footprint,const ImageOver &footprint2,
   Projection &proj, Projection &norm_proj, int Ydim, int Xdim,
   double rot, double tilt, double psi, int FORW, int eq_mode=ARTK,
   GridVolumeT<int> *GVNeq=NULL, matrix2D<double> *M=NULL,
   GridVolumeT<T> *vol_var=NULL, double ray_length=-1);

/** From voxel volumes.
    The voxel volume is projected onto a projection plane defined by
    (rot, tilt, psi) (1st, 2nd and 3rd Euler angles) . The projection
    is previously is resized to Ydim x Xdim and initialized to 0.
    The projection itself, from now on, will keep the Euler angles.
    
    An exception is thrown if the SPIDER environment variable is not set */
void project_Volume(matrix3D<double> &V, Projection &P, int Ydim, int Xdim,
   double rot, double tilt, double psi) _THROW;

/** Count equations in volume.
   For Component AVeraing (CAV), the number of equations in which
   each blob is involved is needed. */
void count_eqs_in_projection(GridVolumeT<int> &GVNeq,
   const ImageOver &footprint, const ImageOver &footprint2,
   Projection &read_proj);
//@}

/**@name Crystal projections */
//@{
/** Project a crystal blob volume.
    This function projects a crystal deformed blob volume, ie, in the
    documentation volume g. However the angles given must be those for
    volume f, the undeformed one. You must supply the deformed lattice
    vectors, and the matrix to pass from the deformed to the undeformed
    vectors (D and Dinv). a=D*ai;
    
    Valid eq_modes are ARTK and CAV.
*/
void project_Crystal_Volume(GridVolume &vol,
   const ImageOver &footprint, const ImageOver &footprint2,
   Projection &proj, Projection &norm_proj,
   int Ydim, int Xdim,
   double rot, double tilt, double psi, const matrix1D<double> &shift,
   const matrix1D<double> &aint, const matrix1D<double> &bint,
   const matrix2D<double> &D, const matrix2D<double> &Dinv,
   const matrix2D<int> &mask, int FORW, int eq_mode=ARTK);
//@}

//@}
#endif
