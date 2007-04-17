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

#ifndef _BLOBS_HH
   #define _BLOBS_HH
   #include "grids.hh"
   #include <XmippData/xmippImages.hh>
   #include <XmippData/xmippVolumes.hh>
   
/* ========================================================================= */
/* BLOBS                                                                     */
/* ========================================================================= */
/**@name Blobs */
//@{
// Blob structure ----------------------------------------------------------
/** Blob definition.
    The blob is a space limited function (click here for a theoretical
    explanation) which is used as basis function for the ART reconstructions.
    There are several parameters which define the shape of the blob.
    The following structure holds all needed information for a blob, a
    variable can be of this type and it is passed to the different functions
    containing all we need to know about the blob. As a type definition,
    we can work with several kind of blobs in the same program at the same
    time.
    
    The common way of defining a blob is as follows:
    \begin{verbatim}
    struct blobtype blob;                  // Definition of the blob
    blob.radius = 2;                       // Blob radius in voxels
    blob.order  = 2;                       // Order of the Bessel function
    blob.alpha  = 3.6;                     // Smoothness parameter
    \end{verbatim}
    
    Sometimes it is useful to plot any quantity related to the blobs. In the
    following example you have how to plot their Fourier transform in the
    continuous frequency space.
    
    \begin{verbatim}
      #include <Reconstruction/blobs.hh>
      #include <XmippData/xmippArgs.hh>

      int main(int argc, char **argv) {
	  struct blobtype blob;                  // Definition of the blob
	  blob.radius = 2;                       // Blob radius in voxels
	  blob.order  = 2;                       // Order of the Bessel function
	  blob.alpha  = AtoF(argv[1]);           // Smoothness parameter


	  double M=blob_Fourier_val (0, blob);
	  for (double w=0; w<=2; w += 0.05)
	     cout << w << " " <<  blob_Fourier_val (w, blob)/M << endl;

	  return 0;
      }
    \end{verbatim}
*/
struct blobtype {
   /// Spatial radius in \Ref{Universal System} units
   double radius;

   /// Derivation order and Bessel function order
   int   order;
   
   /// Smoothness parameter
   double alpha;
};

// Blob value --------------------------------------------------------------
/** Blob value.
    This function returns the value of a blob at a given distance from its
    center (in \Ref{Universal System} units). The distance must be
    always positive. Remember that a blob is spherically symmetrycal so
    the only parameter to know the blob value at a point is its distance
    to the center of the blob. It doesn't matter if this distance is
    larger than the real blob spatial extension, in this case the function
    returns 0 as blob value.
    \\ Ex:
    \begin{verbatim}
    struct blobtype blob; blob.radius = 2; blob.order = 2; blob.alpha = 3.6;
    matrix1D<double> v=vector_R3(1,1,1);
    cout << "Blob value at (1,1,1) = " << blob_val(v.mod(),blob) << endl;
    \end{verbatim} */
#define blob_val(r, blob) kaiser_value(r, blob.radius, blob.alpha, blob.order)
double kaiser_value(double r, double a, double alpha, int m);

// Blob projection ---------------------------------------------------------
/** Blob projection.
    This function returns the value of the blob line integral through a
    straight line which passes at a distance 'r' (in \Ref{Universal System}
    units) from the center of the
    blob. Remember that a blob is spherically symmetrycal so
    the only parameter to know this blob line integral is its distance
    to the center of the blob. It doesn't matter if this distance is
    larger than the real blob spatial extension, in this case the function
    returns 0.
    \\ Ex:
    \begin{verbatim}
    struct blobtype blob; blob.radius = 2; blob.order = 2; blob.alpha = 3.6;
    matrix1D<double> v=vector_R3(1,1,1);
    cout << "Blob line integral through (1,1,1) = " << blob_proj(v.mod(),blob)
         << endl;
    \end{verbatim} */
#define blob_proj(r, blob) kaiser_proj(r, blob.radius, blob.alpha, blob.order)
double kaiser_proj(double r, double a, double alpha, int m);

/** Fourier transform of a blob.
    This function returns the value of the Fourier transform of the blob
    at a given frequency (w). This frequency must be normalized by the
    sampling rate. For instance, for computing the Fourier Transform of 
    a blob at 1/Ts (Ts in Amstrongs) you must provide the frequency Tm/Ts,
    where Tm is the sampling rate.
    
    The Fourier Transform can be computed only for blobs with m=2. */
#define blob_Fourier_val(w, blob) \
   kaiser_Fourier_value(w, blob.radius, blob.alpha, blob.order)
double kaiser_Fourier_value(double w, double a, double alpha, int m);
/** Formula for a volume integral of a blob (n is the blob dimension)
    */
#define blob_mass(blob) \
    basvolume(blob.radius, blob.alpha, blob.order,3)
double  basvolume (double a, double alpha, int m, int n);

/** Limit (z->0) of (1/z)^n I_n(z) (needed by basvolume)*/    
double in_zeroarg(int n);

/** Limit (z->0) of (1/z)^(n+1/2) I_(n+1/2) (z) (needed by basvolume)*/    
double inph_zeroarg(int n);

/** Bessel function I_(n+1/2) (x),  n = 0, 1, 2, ... */
double i_nph(int n, double x);

/** Bessel function I_n (x),  n = 0, 1, 2, ...
	Use ONLY for small values of n	*/
double	i_n(int n, double x);

/** Blob pole.
    This is the normalized frequency at which the blob goes to 0. */
double blob_freq_zero(struct blobtype b);

/** Attenuation of a blob.
    The Fourier transform of the blob at w is the Fourier transform at w=0
    multiplied by the attenuation. This is the value returned. Remind that
    the frequency must be normalized by the sampling rate. Ie, Tm*w(cont) */
double blob_att(double w, struct blobtype b);

/** Number of operations for a blob.
    This is a number proportional to the number of operations that ART
    would need to make a reconstruction with this blob. */
double blob_ops(double w, struct blobtype b);

/** Optimal CC grid spacing.
    This function returns the optimal grid relative size for the blob
    selected. */
double optimal_CC_grid_relative_size(struct blobtype b);

/** Optimal BCC grid spacing.
    This function returns the optimal grid relative size for the blob
    selected. */
double optimal_BCC_grid_relative_size(struct blobtype b);

/** Optimal FCC grid spacing.
    This function returns the optimal grid relative size for the blob
    selected. */
double optimal_FCC_grid_relative_size(struct blobtype b);

/* ========================================================================= */
/* TOOLS                                                                     */
/* ========================================================================= */
/** Choose the best blob within a region.
    You must select a resolution limit, the maximum attenuation allowed at
    that frequency and then the blob with less
    computational cost is returned. m=2 always. The resolution criterion
    is given by w (remind that this w=Tm*w(cont)).
    
    If there is no blob meeting the specification then alpha=a=-1.*/
blobtype best_blob(double alpha_0, double alpha_F, double inc_alpha,
   double a_0, double a_F, double inc_a, double w, double *target_att,
   int target_length=1);

/** Footprint of a blob.
    In the actual implementation of the 3D reconstruction one of the most
    important things to save time is to have a precomputed blob projection.
    As it is spherically symmetrical all projections from different
    directions of the same blob happen to be the same. This function makes
    this precalculation storing it in an oversampled image. This is done
    so because the blob footprint might not be centered with respect to
    the universal grid needing so, the footprint at non-integer positions.
    As parameters to this function you may give the sampling rate in each
    direction, and an optional normalisation (usually disabled) at the end
    of the computation which divides the whole image by the sum of the whole
    image.
    \\ Ex:
    \begin{verbatim}
    ImageOver             blobprint;
    struct blobtype       blob;
    blob.radius = 1.7;
    blob.order  = 2;
    blob.alpha  = 3.6;
    footprint_blob(blobprint, blob);
    \end{verbatim}
    As the blob radius is
    1.7, the function will create a normalised footprint of logical corners
    (-2,-2) to (2,2) (in the universal coord. system) (this size is the minimum
    integer number of pixels which contain entirely the blob), with 50
    samples each pixel. The real size of the footprint is (5·50)x(5·50).
*/
void footprint_blob(ImageOver &blobprint, const struct blobtype &blob,
   int istep=50, int normalise=0);

/** Sum of a single blob over a grid.
    As a normalisation factor, the sum of the blob values over a given grid
    is needed. What this function does is put a blob at coordinate (0,0,0)
    and sums all the blob values at points of the grid which are inside the
    blob. It doesn't matter if the grid is compound or not.
    \\ Ex:
    \begin{verbatim}
    // Blob definition
    struct blobtype       blob;
    blob.radius = 2;
    blob.order  = 2;
    blob.alpha  = 3.6;
    
    // Grid definition
    Grid BCCgrid;
    BCCgrid=BCC_grid(1.41,vector_R3(-5,-5,-5),vector_R3( 5, 5, 5));
    
    cout << "The sum of a single blob over the grid is " <<
         << sum_blob_grid(blob, BCCgrid);
    \end{verbatim}    

   D is a 3x3 matrix specifying a volume deformation. It means that the
   sample at logical position (i,j,k) is really placed at space position
   D*(i,j,k)'.
*/
double sum_blob_Grid(const struct blobtype &blob, const Grid &grid,
   const matrix2D<double> *D=NULL);


/** Voxel shape for a blob volume.
    Given a blob volume this function returns the logical origin and
    size for the minimum voxel volume which holds it. See \Ref{blobs2voxels}
    for an explanation of the limit and V parameters.*/
void voxel_volume_shape(const GridVolume &vol_blobs,
    const struct blobtype &blob, const matrix2D<double> *D,
    matrix1D<int> &corner1, matrix1D<int> &size);

/** Blobs ---> Voxels.
    The voxel size is defined between two coordinates (Gcorner1 and Gcorner2)
    which accomplish
    \begin{verbatim}
       XX(Gcorner1)=MIN(XX(SGcorner1(i)));
       YY(Gcorner1)=MIN(YY(SGcorner1(i)));
       ZZ(Gcorner1)=MIN(ZZ(SGcorner1(i)));

       XX(Gcorner2)=MAX(XX(SGcorner2(i)));
       YY(Gcorner2)=MAX(YY(SGcorner2(i)));
       ZZ(Gcorner2)=MAX(ZZ(SGcorner2(i)));
    \end{verbatim}
    where SGcorner1(i) and SGcorner2(i) are the lowest and highest corners
    of each subgrid. These corners are expressed in Universal coordinates.
   
   This quantity is then enlarged to be an integer voxel position and to take
   into account that the former computed Gcorner1 and 2 are the furthest
   centers of blobs, ie, there will be voxels even further affected by these
   blobs.
   \begin{verbatim}
   Gcorner1 = CEILnD (Gcorner1 - blob.radius);
   Gcorner2 = FLOORnD(Gcorner2 + blob.radius);
   \end{verbatim}
   
   However, you might give a size, usually set to 0, i.e., no external size.
   If no size is provided a size is produced such that all blob centers
   fit into the output volume.
   
   D is a 3x3 matrix specifying a volume deformation. It means that the
   sample at logical position (i,j,k) is really placed at space position
   D*(i,j,k)'.
*/
void blobs2voxels(const GridVolume &vol_blobs,
   const struct blobtype &blob, matrix3D<double> *vol_voxels,
   const matrix2D<double> *D=NULL, int Zdim=0, int Ydim=0, int Xdim=0);

/** Blob coefficients as a voxel volume.
    This function returns a volume with the blob coefficients in their right position
    in space. The voxel size is g/2 */
void blobs2space_coefficients(const GridVolume &vol_blobs, const struct blobtype &blob,
   matrix3D<double> *vol_coefs);

#define SHOW_CONVERSION 1
/** Voxels ---> Blobs.
    The voxels to blobs procedure is an iterative process resolved by
    ART as an iterative technique to solve an equation system. The blobs to
    voxels process is the BASIC way to solve this problem.
    This applies the whole ART process.
    You can set a debugging level by setting the flag SHOW_CONVERSION in the
    tell argument.
    
    If the mask is NULL then it is assumed that all voxels in the input
    voxel volume must be adjusted by the blobs. If the input voxel volume
    is NULL the it is assumed that it is equal to 0 and only those
    corresponding to the 1 positions within the mask must be adjusted.
    
    R is the interest radius. If it is -1 two superimposed grids are created,
    otherwise a single grid with tilted axis is.
    
    Valid grid types are defined in \Ref{Some useful grids}*/
void voxels2blobs(const matrix3D<double> *vol_voxels,
   const struct blobtype &blob, GridVolume &vol_blobs,
   int grid_type, double grid_relative_size,
   double lambda, const matrix3D<double> *vol_mask=NULL,
   const matrix2D<double> *D=NULL,
   double final_error_change=0.01, int tell=0, double R=-1);

#define VARTK 1
#define VMAXARTK 2

/** Voxels ---> Blobs, single step.
    This is the basic step of the voxels to blobs conversion. It applies
    an ART step to the blob volume and produces its output in a different
    (or the same) volume. You must provide the blob structure, lambda and
    deformation matrix (typical for crystals). The theoretical volume and
    correction volume are resized inside and are output volumes meaning the
    conversion from blobs to voxels before this step and the correction
    applied to the blob volume. The read volume is the one to reproduce
    with blobs and the mask volume is used to select only one region to
    adjust.
    
    If the mask is NULL then it is assumed that all voxels in the input
    voxel volume must be adjusted by the blobs. If the input voxel volume
    is NULL the it is assumed that it is equal to 0 and only those
    corresponding to the 1 positions within the mask must be adjusted.
    An exception is thrown if the mask and voxels volume are both NULL.
    
    Then mean and maximum error commited for each blob are returned.
    
    Valid eq_modes are
     \\VARTK:     ART by blocks for volumes.
     \\VMAXARTK:  update with the maximum update.*/
void ART_voxels2blobs_single_step(
   GridVolume &vol_in, GridVolume *vol_out, const struct blobtype &blob,
   const matrix2D<double> *D, double lambda, matrix3D<double> *theo_vol,
   const matrix3D<double> *read_vol, matrix3D<double> *corr_vol,
   const matrix3D<double> *mask_vol,
   double &mean_error, double &max_error, int eq_mode=VARTK);
//@}
#endif
