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
#ifndef _XMIPP_FILTERS_HH
#  define _XMIPP_FILTERS_HH

#define LOG2 0.693147181

#include "xmippImages.hh"
#include "xmippVolumes.hh"
#include "xmippMasks.hh"

/**@name Filters */
//@{
/** Substract background.
    The background is computed as the plane which best fits all density
    values, then this plane is substracted from the image. */
    void substract_background_plane(Image *I);

/** Constrast enhancement.
    The minimum density value is brought to 0 and the maximum to 255. */
    void contrast_enhancement(Image *I);

/** Region growing for images.
    Given a position inside an image, this function grows
    a region with (filling_colour) until it finds a border of value
    (stop_colour). If the point is outside the image
    then nothing is done.
    
    If less is TRUE the region is grown in sucha a way that all
    voxels in its border are greater than the region voxels.
    If less is FALSE the region is grown so that all voxels on its 
    border are smaller than the region voxels.

    Valid neighbourhoods are 4 or 8.
    */
    void region_growing(const matrix2D<double> &I_in, matrix2D<double> &I_out,
       int i, int j,
       float stop_colour=1, float filling_colour=1, bool less=TRUE,
       int neighbourhood=8);
    
/** Region growing for volumes.
    Given a position inside a volume this function grows
    a region with (filling_colour) until it finds a border of value.
    (stop_colour). If the point is outside the volume
    then nothing is done.
    
    If less is TRUE the region is grown in sucha a way that all
    voxels in its border are greater than the region voxels.
    If less is FALSE the region is grown so that all voxels on its 
    border are smaller than the region voxels.
    */
    void region_growing(const matrix3D<double> &V_in, matrix3D<double> &V_out,
       int k, int i, int j,
       float stop_colour=1, float filling_colour=1, bool less=TRUE);

/** Label a binary image.
    This function receives a binary volume and labels all its connected
    components. The background is labeled as 0, and the components as 
    1, 2, 3, ....*/
    int label_image(const matrix2D<double> &I, matrix2D<double> &label,
       int neighbourhood=8);

/** Label a binary volume.
    This function receives a binary image and labels all its connected
    components. The background is labeled as 0, and the components as 
    1, 2, 3, ....*/
    int label_volume(const matrix3D<double> &V, matrix3D<double> &label);

/** Remove connected components smaller than a given size. They are set to 0.*/
    void remove_small_components(matrix2D<double> &I, int size,
       int neighbourhood=8);

/** Keep the biggest connected component.
    If the biggest component does not cover the percentage required
    (by default, 0), more big components are taken until this is
    accomplished.*/
    void keep_biggest_component(matrix2D<double> &I, double percentage=0,
       int neighbourhood=8);

/** Fill object.
    Everything that is not background is assumed to be object. */
    void fill_binary_object(matrix2D<double> &I,int neighbourhood=8);

/** Correlation 1D.
    This function returns the product of both signals in the common
    positions. Notice that it is not the correlation what is usually
    needed but the covariance that is the product of the two signals
    minus their means.
    
    This function returns the number of objects (different from background)*/
    template <class T>
       double correlation(const matrix1D<T> &x, const matrix1D<T> &y, 
	    				  const matrix1D<int> *mask=NULL,int l=0);

/** Correlation 2D. */
    template <class T>
       double correlation(const matrix2D<T> &x, const matrix2D<T> &y,
                          const matrix2D<int> *mask=NULL,int l=0, int m=0);

/** Correlation 3D. */
    template <class T>
       double correlation(const matrix3D<T> &x, const matrix3D<T> &y,
                          const matrix3D<int> *mask=NULL,int l=0, int m=0,
						  int q=0);

/** correlation_index 1D.
    Return the sum{(x-mean_x)*(y-mean_y)}/(stddev_x*stddev_y*n)
	in the common positions. */
    template <class T>
       double correlation_index(const matrix1D<T> &x, const matrix1D<T> &y);

/** correlation_index 2D. */
    template <class T>
       double correlation_index(const matrix2D<T> &x, const matrix2D<T> &y,
           const matrix2D<int> *mask=NULL,
		   matrix2D<double> *Contributions=NULL);

/** correlation_index 3D. */
    template <class T>
       double correlation_index(const matrix3D<T> &x, const matrix3D<T> &y,
           const matrix3D<int> *mask=NULL,
		   matrix3D<double> *Contributions=NULL);

/** Translational search.
    This function returns the best interpolated shift for the alignment
    of two images. You can restrict the shift to a region defined by a mask
    (the maximum will be sought where the mask is 1). */
    void best_shift(const matrix2D<double> &I1, const matrix2D<double> &I2,
       double &shiftX, double &shiftY, const matrix2D<int> *mask=NULL);

/** euclidian_distance 1D.
    Return the SQRT[sum{(x-y)*(x-y)}]
	in the common positions. */
    template <class T>
       double euclidian_distance(const matrix1D<T> &x, const matrix1D<T> &y);

/** euclidian distance 2D. */
    template <class T>
       double euclidian_distance(const matrix2D<T> &x, const matrix2D<T> &y,
           const matrix2D<int> *mask=NULL);

/** euclidian distance 3D. */
    template <class T>
       double euclidian_distance(const matrix3D<T> &x, const matrix3D<T> &y,
           const matrix3D<int> *mask=NULL);

/** mutual information 1D.
    Return the mutual information:
    MI = sum [ P(x,y)*log2{P(x,y)/(P(x)*P(y))} ] 
	in the common positions. 
    P(x), P(y) are 1D-histograms of the values of matrix x and y.
    P(x,y)     is the 2D-histogram, i.e. the count of times that a certain 
               combination of values in matrices x and y has ocurred.
    The sum runs over all histogram bins. 

    The histograms are calculated using the number of bins nx and ny.
    If no values (or zeros) are given, a Gaussian distribution of the values 
    in the matrices is assumed, and the number of bins is calculated as: log2(n)+1.
    (according to: Tourassi et al. (2001) Med. Phys. 28 pp. 2394-2402.) */
    template <class T>
       double mutual_information(const matrix1D<T> &x, const matrix1D<T> &y,
           int nx=0, int ny=0);

/** mutual information 2D. */
    template <class T>
       double mutual_information(const matrix2D<T> &x, const matrix2D<T> &y,
           int nx=0, int ny=0, const matrix2D<int> *mask=NULL);

/** mutual information 3D. */
    template <class T>
       double mutual_information(const matrix3D<T> &x, const matrix3D<T> &y,
           int nx=0, int ny=0, const matrix3D<int> *mask=NULL);

/** RMS 1D.
    Return the sqrt(sum{(x-y)*(x-y)}/n) in the common positions. */
    template <class T>
       double rms(const matrix1D<T> &x, const matrix1D<T> &y,
          const matrix1D<int> *mask=NULL, matrix1D<double> *Contributions=NULL);

/** RMS 2D. */
    template <class T>
       double rms(const matrix2D<T> &x, const matrix2D<T> &y,
          const matrix2D<int> *mask=NULL, matrix2D<double> *Contributions=NULL);

/** RMS 3D. */
    template <class T>
       double rms(const matrix3D<T> &x, const matrix3D<T> &y,
          const matrix3D<int> *mask=NULL,  matrix3D<double> *Contributions=NULL);

/** Fourier-Bessel decomposition.
    The Fourier-Bessel decomposition of those pixels in img_in whose radius is
    between r1 and r2 is computed. r1 and r2 are supposed to fit in the image
    shape, and the image logical origin is used for the decomposition.
    k1 and k2 determines the harmonic coefficients to be computed.
    */
    void Fourier_Bessel_decomposition(const matrix2D<double> &img_in,
       matrix2D<double> &m_out, double r1, double r2, int k1, int k2);

/** Harmonic decomposition.
    */
    void harmonic_decomposition(const matrix2D<double> &img_in,
       matrix1D<double> &v_out);

/** Median_filter with a 3x3 window */
    template <class T> 
       void median_filter3x3(matrix2D<T> &m,matrix2D<T> &out);

/** Mumford-Shah smoothing.
    Purpose:  This function simultaneously smooths and segments an image	 
	      using non-linear diffusion.  Mumford-&-Shah's functional	 
	      minimization algorithm is used to detect region boundaries	  
	      and relax image smoothness constraints near these		 
	      discontinuities.						 
	      The functional minimized is: 				 

	      E = W0*(f-d)*(f-d)	              (data matching)	 
		+ W1*(fx*fx + fy*fy)*(1-s)*(1-s)      (1st deriv smooth) 
		+ W2*(s*s) 			      (edge strengths)   
		+ W3*(sx*sx + sy*sy)		      (edge smoothness)  


    Comments:								 
	      The program diffusion from KUIM (developed by J. Gaush,	 
	      U. Kansas) was used as the "seed".
	      Paper: Teboul, et al. IEEE-Trans. on Image Proc. Vol. 7, 387-397.
*/			 
void Smoothing_Shah(matrix2D<double> &img, 
   matrix2D<double> &surface_strength, matrix2D<double> &edge_strength,
   const matrix1D<double> &W, int OuterLoops, int InnerLoops,
   int RefinementLoops, bool adjust_range=true);

/** Rotational invariant moments.
    The mask and the image are supposed to be of the same shape.
    If no mask is provided, the moments are computed on the whole image.
    The moments are measured with respect to the origin of the image.

    These moments have been taken from
    http://www.cs.cf.ac.uk/Dave/Vision_lecture/node36.html (moments 1 to 5).
*/
void rotational_invariant_moments(const matrix2D<double> &img,
   const matrix2D<int> *mask, matrix1D<double> &v_out);

/** Inertia moments.
    They are measured with respect to the center of the image, and not
    with respect to the center of mass. For an image there are only two
    inertia moments. */
void inertia_moments(const matrix2D<double> &img,
   const matrix2D<int> *mask, matrix1D<double> &v_out);

/** Fill a triangle defined by three points.
    The points are supplied as a pointer to three integer positions. They can
    be negative. */
void fill_triangle(matrix2D<double> &img, int *tx, int *ty, double color);
//@}
#endif
