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
    */
    void region_growing(const Image *I_in, Image *I_out, int i, int j,
       float stop_colour=1, float filling_colour=1, bool less=TRUE);
    
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
    void region_growing(const Volume *V_in, Volume *V_out, int k, int i, int j,
       float stop_colour=1, float filling_colour=1, bool less=TRUE);

/** Label a binary image.
    This function receives a binary volume and labels all its connected
    components. The background is labeled as 0, and the components as 
    1, 2, 3, ....*/
    int label_image(const Image *I, Image *label);

/** Label a binary volume.
    This function receives a binary image and labels all its connected
    components. The background is labeled as 0, and the components as 
    1, 2, 3, ....*/
    int label_volume(const Volume *V, Volume *label);

/** Correlation 1D.
    This function returns the product of both signals in the common
    positions. Notice that it is not the correlation what is usually
    needed but the covariance that is the product of the two signals
    minus their means.
    
    This function returns the number of objects (different from background)*/
    template <class T>
       double correlation(matrix1D<T> &x, matrix1D<T> &y, 
	    				  const matrix1D<int> *mask=NULL,int l=0);

/** Correlation 2D. */
    template <class T>
       double correlation(matrix2D<T> &x, matrix2D<T> &y,
                          const matrix2D<int> *mask=NULL,int l=0, int m=0);

/** Correlation 3D. */
    template <class T>
       double correlation(matrix3D<T> &x, matrix3D<T> &y,
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

//@}
#endif
