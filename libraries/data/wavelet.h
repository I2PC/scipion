/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Antonio Jose Rodriguez Sanchez (ajr@cnb.csic.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
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

#ifndef WAVELET_H
#define WAVELET_H

#include "multidim_array.h"
#include "numerical_recipes.h"

/// @defgroup Wavelets Wavelets
/// @ingroup DataLibrary
//@{
#define HAAR 2
#define DAUB4 4
#define DAUB12 12
#define DAUB20 20

/// @name Wavelet transform
//@{
/** B-spline Wavelet transform of a vector.
 *
 * The B-spline wavelet transform of the input array is computed. The size of
 * the array must be so as to allow the downsampling by 2 as many times as the
 * number of iterations. For instance, if iterations is 1, then it must be a
 * multiple of 2. If iterations is 2, then it must be a multiple of 4. If
 * iterations if 3, then it must be a multiple of 8. And so on.
 *
 * If the isign=-1 then the inverse wavelet transform is performed.
 */
void Bilib_DWT(const MultidimArray< double >& input,
               MultidimArray< double >& result,
               int iterations,
               int isign = 1);

/** Set DWT type.
 *
 * The DWT type should be set before starting making transforms. Valid types
 * are: HAAR, DAUB4, DAUB12, DAUB20
 */
void set_DWT_type(int DWT_type);

/** Get maximum scale.
 *
 * This function returns the maximum scale achievable by the DWT transform of
 * a given size.
 */
inline int Get_Max_Scale(int size)
{
    return ROUND(log10(static_cast< double >(size)) / log10(2.0));
}

/** DWT of a MultidimArray
 *
 * The output vector can be the same as the input one. Previously the type of
 * DWT must be set with set_DWT_type. If isign=1 the direct DWT is performed,
 * if isign=-1 the inverse DWT is done.
 */
template<typename T>
void DWT(const MultidimArray< T >& v, MultidimArray< double >& result, int isign = 1)
{
    unsigned long int nn[3];
    unsigned long int *ptr_nn = nn - 1;
    int dim;
    nn[2] = ZSIZE(v);
    nn[1] = YSIZE(v);
    nn[0] = XSIZE(v);

    if (YSIZE(v) == 1 && ZSIZE(v) == 1)
        dim = 1;
    else if (ZSIZE(v) == 1)
        dim = 2;
    else
        dim = 3;

    typeCast(v, result);
    double* ptr_result = MULTIDIM_ARRAY(result) - 1;
    wtn(ptr_result, ptr_nn, dim, isign, pwt);
}

/** IDWT of a MultidimArray.
 *
 * The output volume can be the same as the input one. Previously the type of
 * DWT must be set with set_DWT_type.
 */
void IDWT(const MultidimArray< double >& v, MultidimArray< double >& result);
//@}

/// @name Wavelet related functions
//@{
/** DWT Low pass versions.
 *
 * This function returns the low pass versions at different scales. The low
 * pass version of the image at scale s is stored in the 01 quadrant of that
 * scale.
 */
void DWT_lowpass2D(const MultidimArray< double >& v, MultidimArray< double >& result);

#define DWT_Imin(s, smax, l) static_cast< int >(((l == '0') ? 0 : \
        pow(2.0, smax - s - 1)))
#define DWT_Imax(s, smax, l) static_cast< int>(((l == '0') ? \
        pow(2.0, smax - s - 1) - 1 : pow(2.0, smax - s) - 1))

/** Select Block 1D.
 *
 * Given the scale (s=0 is the finest) and the quadrant "0" (Lower frequencies)
 * or "1"(Higher frequencies) this routine returns the indices that should be
 * explored for this block (x1 and x2 should be included in the for).
 */
template<typename T>
void SelectDWTBlock(int scale,
                    const MultidimArray< T >& I,
                    const std::string& quadrant,
                    int& x1,
                    int& x2)
{
    double Nx = Get_Max_Scale(XSIZE(I));
    I.toLogical(DWT_Imin(scale, Nx, quadrant[0]), x1);
    I.toLogical(DWT_Imax(scale, Nx, quadrant[0]), x2);
}

/** Select Block 2D.
 *
 * Given the scale (s=0 is the finest) and the quadrant "xy"="00" (Upper left),
 * "01" (Upper right), "10" (Lower left), "11" (Lower right). This routine
 * returns the indices that should be explored for this block (the extremes
 * should be included in the for).
 */
template<typename T>
void SelectDWTBlock(int scale,
                    const MultidimArray< T >& I,
                    const std::string& quadrant,
                    int& x1,
                    int& x2,
                    int& y1,
                    int& y2)
{
    double Nx = Get_Max_Scale(XSIZE(I));
    double Ny = Get_Max_Scale(YSIZE(I));

    x1 = DWT_Imin(scale, Nx, quadrant[0]);
    y1 = DWT_Imin(scale, Ny, quadrant[1]);
    x2 = DWT_Imax(scale, Nx, quadrant[0]);
    y2 = DWT_Imax(scale, Ny, quadrant[1]);

    I.toLogical(y1, x1, y1, x1);
    I.toLogical(y2, x2, y2, x2);
}


/** Select Block 3D.
 *
 * Given the scale (s=0 is the finest) and the quadrant "xyz"="000", "001",
 * "010", "011", "100", "101", "110", "111". This routine returns the indices
 * that should be explored for this block (the extremes should be included in
 * the for).
 */
template<typename T>
void SelectDWTBlock(int scale,
                    const MultidimArray< T >& I,
                    const std::string& quadrant,
                    int& x1,
                    int& x2,
                    int& y1,
                    int& y2,
                    int& z1,
                    int& z2)
{
    double Nx = Get_Max_Scale(XSIZE(I));
    double Ny = Get_Max_Scale(YSIZE(I));
    double Nz = Get_Max_Scale(ZSIZE(I));

    x1 = DWT_Imin(scale, Nx, quadrant[0]);
    y1 = DWT_Imin(scale, Ny, quadrant[1]);
    z1 = DWT_Imin(scale, Nz, quadrant[2]);
    x2 = DWT_Imax(scale, Nx, quadrant[0]);
    y2 = DWT_Imax(scale, Ny, quadrant[1]);
    z2 = DWT_Imax(scale, Nz, quadrant[2]);

    I.toLogical(z1, y1, x1, z1, y1, x1);
    I.toLogical(z2, y2, x2, z2, y2, x2);
}


/** Given a quadrant number it returns the string associated to it.
 *
 * That is nothing more than its corresponding binary representation.
 */
std::string Quadrant2D(int q);

/** Given a quadrant number it returns the string associated to it.
 *
 * That is nothing more than its corresponding binary representation.
 */
std::string Quadrant3D(int q);

/** Get scale and quadrant 1D.
 *
 * Given a point and the maximum size of the image, this routine returns the
 * scale and quadrant it belongs.
 */
void Get_Scale_Quadrant(int size_x, int x, int& scale, std::string& quadrant);

/** Get scale and quadrant 2D.
 *
 * Given a point and the maximum size of the image, this routine returns the
 * scale and quadrant it belongs.
 */
void Get_Scale_Quadrant(int size_x,
                        int size_y,
                        int x,
                        int y,
                        int& scale,
                        std::string& quadrant);

/** Get scale and quadrant 3D.
 *
 * Given a point and the maximum size of the image, this routine returns the
 * scale and quadrant it belongs.
 */
void Get_Scale_Quadrant(int size_x,
                        int size_y,
                        int size_z,
                        int x,
                        int y,
                        int z,
                        int& scale,
                        std::string& quadrant);
//@}

/// @name Wavelet Denoising
//@{

/** Remove all information within a quadrant and scale.
 */
void clean_quadrant2D(MultidimArray< double >& I,
                      int scale,
                      const std::string& quadrant);

/** Remove all information within a quadrant and scale.
 */
void clean_quadrant3D(MultidimArray< double >& I,
                      int scale,
                      const std::string& quadrant);

/** Soft thresholding .
 *
 * Substract a value from all coefficients, if the the value is greater than
 * the absolute value of the coefficient, that coefficient is set to 0.
 */
void soft_thresholding(MultidimArray< double >& I, double th);

/** Adaptive soft thresholding 2D.
 *
 * Chang, Yu, Betterli. IEEE Int. Conf. Image Processing
 */
void adaptive_soft_thresholding2D(MultidimArray< double >& I, int scale);

/** Keep central part 2D.
 *
 * Keep those coefficients in a certain radius.
 */
void DWT_keep_central_part(MultidimArray< double >& I, double R);

/** Bayesian, Wiener filtering.
 *
 * Bijaoui, Signal Processing 2002, 82: 709-712. The Denoising procedure is
 * applied up to the scale given. SNR0 is the smallest SNR and SNRF is the
 * largest SNR.
 *
 * This function returns the estimated coefficients for S and N at each scale.
 * If denoise is set to false, then S and N coefficients are estimated but they
 * are not applied to the image.
 */
Matrix1D< double > bayesian_wiener_filtering2D(MultidimArray< double >& WI,
        int allowed_scale,
        double SNR0 = 0.1,
        double SNRF = 0.2,
        bool white_noise = false,
        int tell = 0,
        bool denoise = true);

/** Bayesian, Wiener filtering.
 *
 * This is the function that really denoise.
 */
void bayesian_wiener_filtering2D(MultidimArray< double >& WI,
                                 int allowed_scale,
                                 Matrix1D< double >& estimatedS);

/** Bayesian, Wiener filtering.
 *
 * Bijaoui, Signal Processing 2002, 82: 709-712. The denoising procedure is
 * applied up to the scale given. SNR0 is the smallest SNR and SNRF is the
 * largest SNR.
 *
 * This function returns the estimated coefficients for S and N at each scale.
 * If denoise is set to false, then S and N coefficients are estimated but they
 * are not applied to the image.
 */
Matrix1D< double > bayesian_wiener_filtering3D(MultidimArray< double >& WI,
        int allowed_scale,
        double SNR0 = 0.1,
        double SNRF = 0.2,
        bool white_noise = false,
        int tell = 0,
        bool denoise = true);

/** Bayesian, Wiener filtering.
 *
 * This is the function that really denoise.
 */
void bayesian_wiener_filtering3D(MultidimArray< double >& WI,
                                 int allowed_scale,
                                 Matrix1D< double >& estimatedS);

/** Phase congruency filtering.
 *
 * Phase congruency of an image using monogenic filters.
 *
 * 	   I                     - Image to be processed
 * 	   PC                    - Phase congruency indicating edge significance
 * 	   FT                    - Local weighted mean phase angle at every point in the
 * 	   						   image.  A value of pi/2 corresponds to a bright line, 0
 *                              corresponds to a step and -pi/2 is a dark line.
 *     nscale            5    - Number of wavelet scales, try values 3-6
 *      minWaveLength    3    - Wavelength of smallest scale filter.
 *      mult             2.1  - Scaling factor between successive filters.
 *      sigmaOnf         0.55 - Ratio of the standard deviation of the Gaussian
 *                              describing the log Gabor filter's transfer function
 *                              in the frequency domain to the filter center frequency.
 *
 *
 * References:
 *
 *    Peter Kovesi, "Image Features From Phase Congruency". Videre: A
 *    Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
 *    Summer 1999 http://mitpress.mit.edu/e-journals/Videre/001/v13.html
 *
 *    Michael Felsberg and Gerald Sommer, "A New Extension of Linear Signal
 *    Processing for Estimating Local Properties and Detecting Features". DAGM
 *    Symposium 2000, Kiel
 *
 *    Michael Felsberg and Gerald Sommer. "The Monogenic Signal" IEEE
 *    Transactions on Signal Processing, 49(12):3136-3144, December 2001
 *
 *    Peter Kovesi, "Phase Congruency Detects Corners and Edges". Proceedings
 *          DICTA 2003, Sydney Dec 10-12
 *
 */
void phaseCongMono(MultidimArray< double >& I,
				   MultidimArray< double >& PC,
				   MultidimArray< double >& FT,
				   MultidimArray< double >& Energy,
				   MultidimArray< double >& lowPass,
				   MultidimArray< double >& Radius,
				   MultidimArray< std::complex <double> >& H,
                                 int nScale,
                                 double minWaveLength,
                                 double mult,
                                 double sigmaOnf
                                 );

//@}
//@}
#endif
