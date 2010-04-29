/***************************************************************************
*
* Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef FILTERS_H
#define FILTERS_H

#define LOG2 0.693147181

#include "image.h"
#include "volume.h"
#include "mask.h"

/// @defgroup Filters Filters
/// @ingroup DataLibrary

/** Substract background
 * @ingroup Filters
 *
 * The background is computed as the plane which best fits all density values,
 * then this plane is substracted from the image.
 */
void substract_background_plane(Matrix2D<double> &I);

/** Substract background
 * @ingroup Filters
 *
 * The background is computed as a rolling ball operation with a ball
 * with this radius. The radius is typically Xdim/10.
 *
 * This code has been implemented after the one of "Subtract background" in
 * ImageJ.
 */
void substract_background_rolling_ball(Matrix2D<double> &I, int radius);

/** Constrast enhancement
 * @ingroup Filters
 *
 * The minimum density value is brought to 0 and the maximum to 255.
 */
void contrast_enhancement(Image* I);

/** Region growing for images
 * @ingroup Filters
 *
 * Given a position inside an image, this function grows a region with
 * (filling_colour) until it finds a border of value (stop_colour). If the point
 * is outside the image then nothing is done.
 *
 * If less is true the region is grown in sucha a way that all voxels in its
 * border are greater than the region voxels. If less is false the region is
 * grown so that all voxels on its border are smaller than the region voxels.
 *
 * Valid neighbourhoods are 4 or 8.
 */
void region_growing(const Matrix2D< double >& I_in,
                    Matrix2D< double >& I_out,
                    int i,
                    int j,
                    float stop_colour = 1,
                    float filling_colour = 1,
                    bool less = true,
                    int neighbourhood = 8);

/** Region growing for volumes
 * @ingroup Filter
 *
 * Given a position inside a volume this function grows a region with
 * (filling_colour) until it finds a border of value (stop_colour). If the point
 * is outside the volume then nothing is done.
 *
 * If less is true the region is grown in sucha a way that all voxels in its
 * border are greater than the region voxels. If less is false the region is
 * grown so that all voxels on its border are smaller than the region voxels.
 */
void region_growing(const Matrix3D< double >& V_in,
                    Matrix3D< double >& V_out,
                    int k,
                    int i,
                    int j,
                    float stop_colour = 1,
                    float filling_colour = 1,
                    bool less = true);

/** L1 distance transform
  * @ingroup Filters
  *
  * If wrap is set, the image borders are wrapped around.
  * This is useful if the image coordinates represent angles
  */
void distance_transform(const Matrix2D<int> &in,
                        Matrix2D<int> &out, bool wrap=false);

/** Label a binary image
 * @ingroup Filters
 *
 * This function receives a binary volume and labels all its connected
 * components. The background is labeled as 0, and the components as 1, 2, 3
 * ...
 */
int label_image(const Matrix2D< double >& I,
                Matrix2D< double >& label,
                int neighbourhood = 8);

/** Label a binary volume
 * @ingroup Filters
 *
 * This function receives a binary image and labels all its connected
 * components. The background is labeled as 0, and the components as 1, 2, 3
 * ...
 */
int label_volume(const Matrix3D< double >& V, Matrix3D< double >& label);

/** Remove connected components
 * @ingroup Filters
 *
 * Remove connected components smaller than a given size. They are set to 0.
 */
void remove_small_components(Matrix2D< double >& I,
                             int size,
                             int neighbourhood = 8);

/** Keep the biggest connected component
 * @ingroup Filters
 *
 * If the biggest component does not cover the percentage required (by default,
 * 0), more big components are taken until this is accomplished.
 */
void keep_biggest_component(Matrix2D< double >& I,
                            double percentage = 0,
                            int neighbourhood = 8);

/** Fill object
 * @ingroup Filters
 *
 * Everything that is not background is assumed to be object.
 */
void fill_binary_object(Matrix2D< double >&I, int neighbourhood = 8);

/** Segment an object using Otsu's method
 * @ingroup Filters
 *
 * Otsu's method determines a threshold such that the variance of the
 * two classes is minimized
 *
 * http://www.biomecardio.com/matlab/otsu.html
 */
void OtsuSegmentation(Matrix3D<double> &V);

/** Segment an object using Entropy method
 * @ingroup Filters
 *
 * Entropy method determines a threshold such that the entropy of the
 * two classes is maximized
 *
 * http://rsbweb.nih.gov/ij/plugins/download/Entropy_Threshold.java
 */
void EntropySegmentation(Matrix3D<double> &V);

/** Segment an object using a combination of Otsu and Entropy method
 * @ingroup Filters
 *
 * The combination aims at minimizing Z(t)=-log10(sigma2B(t))/H(t)
 * Minimizing the intraclass variance in Otsu is the same as
 * maximizing sigma2B. H is the entropy of the two classes in the entropy
 * method.
 *
 * Then, the lowest percentil of Z is computed. The threshold applied to the
 * volume is the first value in the curve Z(t) falling below this
 * percentil.
 */
void EntropyOtsuSegmentation(Matrix3D<double> &V, double percentil=0.05);

/** Correlation 1D
 * @ingroup Filters
 *
 * This function returns the product of both signals in the common positions.
 * Notice that it is not the correlation what is usually needed but the
 * covariance that is the product of the two signals minus their means.
 *
 * This function returns the number of objects (different from background)
 */
template <typename T>
double correlation(const Matrix1D< T >& x,
                   const Matrix1D< T >& y,
                   const Matrix1D< int >* mask = NULL,
                   int l = 0)
{
    SPEED_UP_temps;

    double retval = 0; // returned value
    int i, ip; // indexes
    int Rows; // of the matrices

    Rows = XSIZE(x);
    for (i = 0; i < Rows; i++)
    {
        ip = i - l;
        if (ip >= 0 && ip < Rows)
        {
            if (mask != NULL)
                if (!DIRECT_VEC_ELEM((*mask), i))
                    continue;

            retval += DIRECT_VEC_ELEM(x, i) * DIRECT_VEC_ELEM(y, ip);
        }
    }

    return retval / Rows;
}

/** Correlation 2D
 * @ingroup Filters
 */
template <typename T>
double correlation(const Matrix2D< T >& x,
                   const Matrix2D< T >& y,
                   const Matrix2D< int >* mask = NULL,
                   int l = 0, int m = 0)
{
    /* Note: l index is for rows and m index for columns */

    SPEED_UP_temps;

    double retval = 0; // returned value
    int i, j, ip, jp; // indexes
    int Rows, Cols; // of the matrices

    // do the computation
    Cols = x.colNumber();
    Rows = x.rowNumber();

    for (i = 0; i < Rows; i++)
        for (j = 0; j < Cols; j++)
        {
            ip = i - l;
            jp = j - m;

            if (ip >= 0 && ip < Rows && jp >= 0 && jp < Cols)
                if (mask != NULL)
                    if (!DIRECT_MAT_ELEM((*mask), i, j))
                        continue;

            retval += DIRECT_MAT_ELEM(x, i, j) * DIRECT_MAT_ELEM(y, ip, jp);
        }

    return retval / (Cols * Rows);
}

/** Correlation 3D
 * @ingroup Filters
 */
template <typename T>
double correlation(const Matrix3D< T >& x,
                   const Matrix3D< T >& y,
                   const Matrix3D< int >* mask = NULL,
                   int l = 0,
                   int m = 0,
                   int q = 0)
{
    SPEED_UP_temps;

    double retval = 0; // returned value
    int i, j, k, ip, jp, kp; // indexes
    int Rows, Cols, Slices; // of the volumes

    // do the computation
    Cols = x.colNumber();
    Rows = x.rowNumber();
    Slices = x.sliceNumber();

    long N = 0;
    for (k = 0; k < Slices; k++)
        for (i = 0; i < Rows; i++)
            for (j = 0; j < Cols; j++)
            {
                ip = i - l;
                jp = j - m;
                kp = k - q;

                if (ip >= 0 && ip < Rows && jp >= 0 && jp < Cols && kp >= 0 &&
                    kp < Slices)
                {
                    if (mask != NULL)
                        if (!DIRECT_VOL_ELEM((*mask), k, i, j))
                            continue;

                    retval += DIRECT_VOL_ELEM(x, k, i, j) *
                              DIRECT_VOL_ELEM(y, kp, ip, jp);
                }
            }

    return retval / (Slices * Rows * Cols);
}

/** correlation_index 1D
 * @ingroup Filters
 *
 * Return the sum{(x-mean_x)*(y-mean_y)}/(stddev_x*stddev_y*n) in the common
 * positions.
 */
template <typename T>
double correlation_index(const Matrix1D< T >& x, const Matrix1D< T >& y)
{
    SPEED_UP_temps;

    double retval = 0;
    double mean_x, mean_y;
    double stddev_x, stddev_y;
    double aux;
    T dummy;
    long n = 0;

    x.computeStats(mean_x, stddev_x, dummy, dummy);
    y.computeStats(mean_y, stddev_y, dummy, dummy);
    if (ABS(stddev_x)<XMIPP_EQUAL_ACCURACY ||
        ABS(stddev_y)<XMIPP_EQUAL_ACCURACY) return 0;

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x, y)
    {
        retval += (VEC_ELEM(x, i) - mean_x) * (VEC_ELEM(y, i) - mean_y);
        n++;
    }

    if (n != 0)
        return retval / ((stddev_x * stddev_y) * n);
    else
        return 0;
}

/** correlation_index 2D
 * @ingroup Filters
 */
template <typename T>
double correlation_index(const Matrix2D< T >& x,
                         const Matrix2D< T >&y,
                         const Matrix2D< int >* mask = NULL,
                         Matrix2D< double >* Contributions = NULL)
{
    SPEED_UP_temps;

    double retval = 0, aux;
    double mean_x, mean_y;
    double stddev_x, stddev_y;
    T dummy;
    long n = 0;

    if (mask == NULL)
    {
        x.computeStats(mean_x, stddev_x, dummy, dummy);
        y.computeStats(mean_y, stddev_y, dummy, dummy);
    }
    else
    {
        computeStats_within_binary_mask(*mask, x, dummy, dummy, mean_x,
                                        stddev_x);
        computeStats_within_binary_mask(*mask, y, dummy, dummy, mean_y,
                                        stddev_y);
    }
    if (ABS(stddev_x)<XMIPP_EQUAL_ACCURACY ||
        ABS(stddev_y)<XMIPP_EQUAL_ACCURACY) return 0;

    // If contributions are desired. Please, be careful interpreting individual
    // contributions to the covariance! One pixel value afect others.
    if (Contributions != NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(i, j))
                    continue;

            aux = (MAT_ELEM(x, i, j) - mean_x) * (MAT_ELEM(y, i, j) - mean_y);
            MAT_ELEM(*Contributions, i, j) = aux;
            retval += aux;
            n++;
        }

        FOR_ALL_ELEMENTS_IN_MATRIX2D(*Contributions)
        MAT_ELEM(*Contributions, i, j) /= ((stddev_x * stddev_y) * n);

    }
    // In other case, normal process (less computationaly expensive)
    else
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(i, j))
                    continue;

            retval += (MAT_ELEM(x, i, j) - mean_x) * (MAT_ELEM(y, i, j) -
                      mean_y);
            n++;
        }
    }

    if (n != 0)
        return retval / ((stddev_x * stddev_y) * n);
    else
        return 0;
}

/** correlation_index 3D
 * @ingroup Filters
 */
template <typename T>
double correlation_index(const Matrix3D< T >& x,
                         const Matrix3D< T >& y,
                         const Matrix3D< int >* mask = NULL,
                         Matrix3D< double >* Contributions = NULL)
{
    SPEED_UP_temps;

    double retval = 0, aux;
    T dummy;
    double mean_x, mean_y;
    double stddev_x, stddev_y;

    long n = 0;

    if (mask == NULL)
    {
        x.computeStats(mean_x, stddev_x, dummy, dummy);
        y.computeStats(mean_y, stddev_y, dummy, dummy);
    }
    else
    {
        computeStats_within_binary_mask(*mask, x, dummy, dummy, mean_x,
                                        stddev_x);
        computeStats_within_binary_mask(*mask, y, dummy, dummy, mean_y,
                                        stddev_y);
    }
    if (ABS(stddev_x)<XMIPP_EQUAL_ACCURACY ||
        ABS(stddev_y)<XMIPP_EQUAL_ACCURACY) return 0;

    // If contributions are desired. Please, be careful interpreting individual
    // contributions to the covariance! One pixel value afect others.
    if (Contributions != NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(k, i, j))
                    continue;

            aux = (VOL_ELEM(x, k, i, j) - mean_x) * (VOL_ELEM(y, k, i, j) -
                    mean_y);
            VOL_ELEM(*Contributions, k, i, j) = aux;
            retval += aux;
            n++;
        }

        FOR_ALL_ELEMENTS_IN_MATRIX3D(*Contributions)
        VOL_ELEM(*Contributions, k, i, j) /= ((stddev_x * stddev_y) * n);
    }
    else
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(k, i, j))
                    continue;

            retval += (VOL_ELEM(x, k, i, j) - mean_x) * (VOL_ELEM(y, k, i, j) -
                      mean_y);
            n++;
        }
    }

    if (n != 0)
        return retval / ((stddev_x * stddev_y) * n);
    else
        return 0;
}

/** Correntropy 1D
 * @ingroup Filters
 */
template <typename T>
double correntropy(const Matrix1D<T> &x, const Matrix1D<T> &y,
                   double sigma)
{
    double retval=0;
    double K=-0.5/(sigma*sigma);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(x)
    {
        double diff=DIRECT_VEC_ELEM(x,i)-DIRECT_VEC_ELEM(y,i);
        retval+=exp(K*diff*diff);
    }
    retval/=XSIZE(x);
    return retval;
}

/** Fast Correntropy 1D
 * @ingroup Filters
 */
template <typename T>
double fastCorrentropy(const Matrix1D<T> &x, const Matrix1D<T> &y,
                       double sigma, const GaussianInterpolator &G)
{
    double retval=0;
    double isigma=1.0/sigma;
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(x)
    retval+=G.getValue(isigma*(DIRECT_VEC_ELEM(x,i)-DIRECT_VEC_ELEM(y,i)));
    retval/=XSIZE(x);
    return retval;
}

/** Correntropy 2D
 * @ingroup Filters
 */
template <typename T>
double correntropy(const Matrix2D<T> &x, const Matrix2D<T> &y,
                   double sigma)
{
    double retval=0;
    double K=-0.5/(sigma*sigma);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(x)
    {
        double diff=DIRECT_MAT_ELEM(x,i,j)-DIRECT_MAT_ELEM(y,i,j);
        retval+=exp(K*diff*diff);
    }
    retval/=XSIZE(x)*YSIZE(x);
    return retval;
}

/** Correntropy 3D
 * @ingroup Filters
 */
template <typename T>
double correntropy(const Matrix3D<T> &x, const Matrix3D<T> &y,
                   double sigma)
{
    double retval=0;
    double K=-0.5/(sigma*sigma);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(x)
    {
        double diff=DIRECT_VOL_ELEM(x,k,i,j)-DIRECT_VOL_ELEM(y,k,i,j);
        retval+=exp(K*diff*diff);
    }
    retval/=XSIZE(x)*YSIZE(x)*ZSIZE(x);
    return retval;
}

/** Translational search
 * @ingroup Filters
 *
 * This function returns the best interpolated shift for the alignment of two
 * images. You can restrict the shift to a region defined by a mask (the maximum
 * will be sought where the mask is 1).
 *
 * To apply these results you must shift I1 by (-shiftX,-shiftY) or
 * I2 by (shiftX, shiftY)
 */
void best_shift(const Matrix2D< double >& I1,
                const Matrix2D< double >& I2,
                double& shiftX,
                double& shiftY,
                const Matrix2D< int >* mask = NULL);

/** Translational search (non-wrapping)
 * @ingroup Filters
 *
 * This function returns the best interpolated shift for the alignment of two
 * images. You can restrict the shift to a region defined by a mask (the maximum
 * will be sought where the mask is 1).
 */
void best_nonwrapping_shift(const Matrix2D< double >& I1,
                            const Matrix2D< double >& I2,
                            double& shiftX,
                            double& shiftY);

/** Align two images
* @ingroup Filters
 *
 * This function modifies I2 to be aligned with I1. Translational and
 * rotational alignments are both considered. The matrix transforming I2
 * into I1 is returned.
 */
void alignImages(const Matrix2D< double >& Iref,
                 Matrix2D< double >& I,
                 Matrix2D< double >&M);

/** Unnormalized 2D gaussian value using covariance
 * @ingroup NumericalFunctions
 *
 * This function returns the value of a multivariate (2D) gaussian function at
 * the point r (column vector of dimension 2).
 *
 * G(r,mu,sigma)=exp(-0.5 * (r-mu)^t sigma^-1 (r-mu))
 */
double unnormalizedGaussian2D(const Matrix1D<double> &r,
                              const Matrix1D<double> &mu,
                              const Matrix2D<double> &sigmainv);

/** Fit Gaussian spot to an image.
 * @ingroup Filters
 *
 * The fitted Gaussian is a*G(r,mu,sigma)+b where
 * G(r,mu,sigma)=exp(-0.5 * (r-mu)^t sigma^-1 (r-mu))
 *
 * You can choose if the center is estimated or it is assumed to be 0.
 * You can choose the number of iterations for the estiamtion.
 */
void estimateGaussian2D(const Matrix2D<double> &I,
                        double &a, double &b, Matrix1D<double> &mu, Matrix2D<double> &sigma,
                        bool estimateMu=true, int iterations=10);

/** euclidian_distance 1D
 * @ingroup Filters
 *
 * Return the SQRT[sum{(x-y)*(x-y)}] in the common positions.
 */
template <typename T>
double euclidian_distance(const Matrix1D< T >& x, const Matrix1D< T >& y)
{
    SPEED_UP_temps;

    double retval = 0;
    long n = 0;

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x, y)
    {
        retval += (VEC_ELEM(x, i) - VEC_ELEM(y, i)) * (VEC_ELEM(x, i) -
                  VEC_ELEM(y, i));
        n++;
    }

    if (n != 0)
        return sqrt(retval);
    else
        return 0;
}

/** euclidian distance 2D
 * @ingroup Filters
 *
 */
template <typename T>
double euclidian_distance(const Matrix2D< T >& x,
                          const Matrix2D< T >& y,
                          const Matrix2D< int >* mask = NULL)
{
    SPEED_UP_temps;

    double retval = 0;
    long n = 0;

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x, y)
    {
        if (mask != NULL)
            if (!(*mask)(i, j))
                continue;

        retval += (MAT_ELEM(x, i, j) - MAT_ELEM(y, i, j)) * (MAT_ELEM(x, i, j)
                  - MAT_ELEM(y, i, j));
        n++;
    }

    if (n != 0)
        return sqrt(retval);
    else
        return 0;
}

/** euclidian distance 3D
 * @ingroup Filters
 */
template <typename T>
double euclidian_distance(const Matrix3D< T >& x,
                          const Matrix3D< T >& y,
                          const Matrix3D< int >* mask = NULL)
{
    SPEED_UP_temps;

    double retval = 0;
    long n = 0;

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x, y)
    {
        if (mask != NULL)
            if (!(*mask)(k, i, j))
                continue;

        retval += (VOL_ELEM(x, k, i, j) - VOL_ELEM(y, k, i, j)) *
                  (VOL_ELEM(x, k, i, j) - VOL_ELEM(y, k, i, j));
        n++;
    }

    if (n != 0)
        return sqrt(retval);
    else
        return 0;
}

/** mutual information 1D
 * @ingroup Filters
 *
 * Return the mutual information:
 * MI = sum [ P(x,y)*log2{P(x,y)/(P(x)*P(y))} ]
 * in the common positions.
 * P(x), P(y) are 1D-histograms of the values of matrix x and y.
 * P(x,y)     is the 2D-histogram, i.e. the count of times that a certain
 *            combination of values in matrices x and y has ocurred.
 * The sum runs over all histogram bins.
 *
 * The histograms are calculated using the number of bins nx and ny. If no
 * values (or zeros) are given, a Gaussian distribution of the values in the
 * matrices is assumed, and the number of bins is calculated as: log2(n)+1.
 * (according to: Tourassi et al. (2001) Med. Phys. 28 pp. 2394-2402.)
 */
template <typename T>
double mutual_information(const Matrix1D< T >& x,
                          const Matrix1D< T >& y,
                          int nx = 0,
                          int ny = 0)
{
    SPEED_UP_temps;

    long n = 0;
    histogram1D histx, histy;
    histogram2D histxy;
    Matrix1D< T > aux_x, aux_y;
    Matrix1D< double > mx, my;
    Matrix2D< double > mxy;
    int xdim, ydim;
    double MI = 0.0;
    double HAB = 0.0;
    double retval = 0.0;

    aux_x.resize(x);
    aux_y.resize(y);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x, y)
    {
        aux_x(n) = VEC_ELEM(x, i);
        aux_y(n) = VEC_ELEM(y, i);
        n++;
    }

    aux_x.resize(n);
    aux_y.resize(n);

    if (n != 0)
    {
        if (nx == 0)
            //Assume Gaussian distribution
            nx = (int)((log((double) n) / LOG2) + 1);

        if (ny == 0)
            //Assume Gaussian distribution
            ny = (int)((log((double) n) / LOG2) + 1);

        compute_hist(aux_x, histx, nx);
        compute_hist(aux_y, histy, ny);
        compute_hist(aux_x, aux_y, histxy, nx, ny);

        mx = histx;
        my = histy;
        mxy = histxy;

        for (int i = 0; i < nx; i++)
        {
            double histxi = (histx(i)) / n;
            for (int j = 0; j < ny; j++)
            {
                double histyj = (histy(j)) / n;
                double histxyij = (histxy(i, j)) / n;

                if (histxyij > 0)
                    retval += histxyij * log(histxyij / (histxi * histyj)) /
                              LOG2;
            }
        }

        return retval;
    }
    else
        return 0;
}

/** mutual information 2D
 * @ingroup Filters
 */
template <typename T>
double mutual_information(const Matrix2D< T >& x,
                          const Matrix2D< T >& y,
                          int nx = 0,
                          int ny = 0,
                          const Matrix2D< int >* mask = NULL)
{
    SPEED_UP_temps;

    long n = 0;
    histogram1D histx, histy;
    histogram2D histxy;
    Matrix1D< T > aux_x, aux_y;
    Matrix1D< double > mx, my;
    Matrix2D< double > mxy;
    int xdim, ydim;
    double retval = 0.0;

    xdim=XSIZE(x);
    ydim=YSIZE(x);
    aux_x.resize(xdim * ydim);
    xdim=XSIZE(y);
    ydim=YSIZE(y);
    aux_y.resize(xdim * ydim);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x, y)
    {
        if (mask != NULL)
            if (!(*mask)(i, j))
                continue;

        aux_x(n) = MAT_ELEM(x, i, j);
        aux_y(n) = MAT_ELEM(y, i, j);
        n++;
    }

    aux_x.resize(n);
    aux_y.resize(n);

    if (n != 0)
    {
        if (nx == 0)
            //Assume Gaussian distribution
            nx = (int)((log((double) n) / LOG2) + 1);

        if (ny == 0)
            //Assume Gaussian distribution
            ny = (int)((log((double) n) / LOG2) + 1);

        compute_hist(aux_x, histx, nx);
        compute_hist(aux_y, histy, ny);
        compute_hist(aux_x, aux_y, histxy, nx, ny);

        mx = histx;
        my = histy;
        mxy = histxy;

        for (int i = 0; i < nx; i++)
        {
            double histxi = (histx(i)) / n;
            for (int j = 0; j < ny; j++)
            {
                double histyj = (histy(j)) / n;
                double histxyij = (histxy(i, j)) / n;
                if (histxyij > 0)
                    retval += histxyij * log(histxyij / (histxi * histyj)) /
                              LOG2;
            }
        }

        return retval;
    }
    else
        return 0;
}

/** mutual information 3D
 * @ingroup Filters
 */
template <typename T>
double mutual_information(const Matrix3D< T >& x,
                          const Matrix3D< T >& y,
                          int nx = 0,
                          int ny = 0,
                          const Matrix3D< int >* mask = NULL)
{
    SPEED_UP_temps;

    long n = 0;
    histogram1D histx, histy;
    histogram2D histxy;
    Matrix1D< T > aux_x, aux_y;
    Matrix1D< double > mx, my;
    Matrix2D< double > mxy;
    int xdim, ydim, zdim;
    double retval = 0.0;

    xdim=XSIZE(x);
    ydim=YSIZE(x);
    zdim=ZSIZE(x);
    aux_x.resize(xdim * ydim * zdim);
    xdim=XSIZE(y);
    ydim=YSIZE(y);
    zdim=ZSIZE(y);
    aux_y.resize(xdim * ydim * zdim);

    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x, y)
    {
        if (mask != NULL)
            if (!(*mask)(k, i, j))
                continue;

        aux_x(n) = VOL_ELEM(x, k, i, j);
        aux_y(n) = VOL_ELEM(y, k, i, j);
        n++;
    }

    aux_x.resize(n);
    aux_y.resize(n);

    if (n != 0)
    {
        if (nx == 0)
            //Assume Gaussian distribution
            nx = (int)((log((double) n) / LOG2) + 1);

        if (ny == 0)
            //Assume Gaussian distribution
            ny = (int)((log((double) n) / LOG2) + 1);

        compute_hist(aux_x, histx, nx);
        compute_hist(aux_y, histy, ny);
        compute_hist(aux_x, aux_y, histxy, nx, ny);

        mx = histx;
        my = histy;
        mxy = histxy;
        for (int i = 0; i < nx; i++)
        {
            double histxi = (histx(i)) / n;
            for (int j = 0; j < ny; j++)
            {
                double histyj = (histy(j)) / n;
                double histxyij = (histxy(i, j)) / n;
                if (histxyij > 0)
                    retval += histxyij * log(histxyij / (histxi * histyj)) /
                              LOG2;
            }
        }

        return retval;
    }
    else
        return 0;
}

/** RMS 1D
 * @ingroup Filters
 *
 * Return the sqrt(sum{(x-y)*(x-y)}/n) in the common positions.
 */
template <typename T>
double rms(const Matrix1D< T >& x,
           const Matrix1D< T >& y,
           const Matrix1D< int >* mask = NULL,
           Matrix1D< double >* Contributions = NULL)
{
    SPEED_UP_temps;

    double retval = 0, aux;
    int n = 0;

    // If contributions are desired
    if (Contributions != NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(i))
                    continue;

            aux = (VEC_ELEM(x, i) - VEC_ELEM(y, i)) * (VEC_ELEM(x, i) -
                    VEC_ELEM(y, i));
            VEC_ELEM(*Contributions, i) = aux;
            retval += aux;
            n++;
        }

        FOR_ALL_ELEMENTS_IN_MATRIX1D(*Contributions)
        VEC_ELEM(*Contributions, i) = sqrt(VEC_ELEM(*Contributions, i) / n);
    }
    // In other case, normal process (less computationaly expensive)
    else
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(i))
                    continue;

            retval += (VEC_ELEM(x, i) - VEC_ELEM(y, i)) * (VEC_ELEM(x, i) -
                      VEC_ELEM(y, i));
            n++;
        }
    }

    if (n != 0)
        return sqrt(retval / n);
    else
        return 0;
}

/** RMS 2D
 * @ingroup Filters
 */
template <typename T>
double rms(const Matrix2D< T >& x,
           const Matrix2D< T >& y,
           const Matrix2D< int >* mask = NULL,
           Matrix2D< double >* Contributions = NULL)
{
    SPEED_UP_temps;

    double retval = 0, aux;
    int n = 0;

    // If contributions are desired
    if (Contributions != NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(i, j))
                    continue;

            aux = (MAT_ELEM(x, i, j) - MAT_ELEM(y, i, j)) * (MAT_ELEM(x, i, j)
                    - MAT_ELEM(y, i, j));
            MAT_ELEM(*Contributions, i, j) = aux;
            retval += aux;
            n++;
        }

        FOR_ALL_ELEMENTS_IN_MATRIX2D(*Contributions)
        MAT_ELEM(*Contributions, i, j) = sqrt(MAT_ELEM(*Contributions, i, j) /
                                              n);
    }
    // In other case, normal process (less computationaly expensive)
    else
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(i, j))
                    continue;

            retval += (MAT_ELEM(x, i, j) - MAT_ELEM(y, i, j)) *
                      (MAT_ELEM(x, i, j) - MAT_ELEM(y, i, j));
            n++;
        }
    }

    if (n != 0)
        return sqrt(retval / n);
    else
        return 0;
}

/** RMS 3D
 * @ingroup Filters
 */
template <typename T>
double rms(const Matrix3D< T >& x,
           const Matrix3D< T >& y,
           const Matrix3D< int >* mask = NULL,
           Matrix3D< double >* Contributions = NULL)
{
    SPEED_UP_temps;

    double retval = 0;
    double aux;
    int n = 0;

    // If contributions are desired
    if (Contributions != NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(k, i, j))
                    continue;

            aux = (VOL_ELEM(x, k, i, j) - VOL_ELEM(y, k, i, j)) *
                  (VOL_ELEM(x, k, i, j) - VOL_ELEM(y, k, i, j));
            VOL_ELEM(*Contributions, k, i, j) = aux;
            retval += aux;
            n++;
        }

        FOR_ALL_ELEMENTS_IN_MATRIX3D(*Contributions)
        VOL_ELEM(*Contributions, k, i, j) = sqrt(VOL_ELEM(*Contributions,
                                            k, i, j) / n);
    }
    else
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(k, i, j))
                    continue;

            retval += (VOL_ELEM(x, k, i, j) - VOL_ELEM(y, k, i, j)) *
                      (VOL_ELEM(x, k, i, j) - VOL_ELEM(y, k, i, j));
            n++;
        }
    }

    if (n != 0)
        return sqrt(retval / n);
    else
        return 0;
}

/** Fourier-Bessel decomposition
 * @ingroup Filters
 *
 * The Fourier-Bessel decomposition of those pixels in img_in whose radius is
 * between r1 and r2 is computed. r1 and r2 are supposed to fit in the image
 * shape, and the image logical origin is used for the decomposition.
 * k1 and k2 determines the harmonic coefficients to be computed.
 */
void Fourier_Bessel_decomposition(const Matrix2D< double >& img_in,
                                  Matrix2D< double >& m_out,
                                  double r1,
                                  double r2,
                                  int k1,
                                  int k2);

/** Harmonic decomposition
 * @ingroup Filters
 */
void harmonic_decomposition(const Matrix2D< double >& img_in,
                            Matrix1D< double >& v_out);

// TODO Document, check indentation
template <typename T>
void sort(T a, T b, T c, Matrix1D< T >& v)
{
    if (a < b)
        if (b < c)
        {
            v(0) = a;
            v(1) = b;
            v(2) = c;
        }
        else if (a < c)
        {
            v(0) = a;
            v(1) = c;
            v(2) = b;
        }
        else
        {
            v(0) = c;
            v(1) = a;
            v(2) = b;
        }
    else if (a < c)
    {
        v(0) = b;
        v(1) = a;
        v(2) = c;
    }
    else if (b < c)
    {
        v(0) = b;
        v(1) = c;
        v(2) = a;
    }
    else
    {
        v(0) = c;
        v(1) = b;
        v(2) = a;
    }
}

template <typename T>
void merge_sort(Matrix1D< T >& v1, Matrix1D< T >& v2, Matrix1D< T >& v)
{
    int i1 = 0, i2 = 0, i = 0;

    while ((i1 < 3) && (i2 < 3))
    {
        if (v1(i1) < v2(i2))
            v(i++) = v1(i1++);
        else
            v(i++) = v2(i2++);
    }

    while (i1 < 3)
        v(i++) = v1(i1++);

    while (i2 < 3)
        v(i++) = v2(i2++);
}

// This UGLY function performs a fast merge sort for the case of vectors of 3
// elements. This way is guaranteed a minimum number of comparisons (maximum
// number of comparisons to perform the sort, 5)
template <typename T>
void fast_merge_sort(Matrix1D< T >& x, Matrix1D< T >& y, Matrix1D< T >& v)
{
    if (x(0) < y(0))
    {
        v(0) = x(0);
        if (x(1) < y(0))
        {
            v(1) = x(1);
            if (x(2) < y(0))
            {
                v(2) = x(2);
                v(3) = y(0);
                v(4) = y(1);
                v(5) = y(2);
            }
            else
            {
                v(2) = y(0);
                if (x(2) < y(1))
                {
                    v(3) = x(2);
                    v(4) = y(1);
                    v(5) = y(2);
                }
                else
                {
                    v(3) = y(1);
                    if (x(2) < y(2))
                    {
                        v(4) = x(2);
                        v(5) = y(2);
                    }
                    else
                    {
                        v(4) = y(2);
                        v(5) = x(2);
                    }
                }
            }
        }
        else
        {
            v(1) = y(0);
            if (x(1) < y(1))
            {
                v(2) = x(1);
                if (x(2) < y(1))
                {
                    v(3) = x(2);
                    v(4) = y(1);
                    v(5) = y(2);
                }
                else
                {
                    v(3) = y(1);
                    if (x(2) < y(2))
                    {
                        v(4) = x(2);
                        v(5) = y(2);
                    }
                    else
                    {
                        v(4) = y(2);
                        v(5) = x(2);
                    }
                }
            }
            else
            {
                v(2) = y(1);
                if (x(1) < y(2))
                {
                    v(3) = x(1);
                    if (x(2) < y(2))
                    {
                        v(4) = x(2);
                        v(5) = y(2);
                    }
                    else
                    {
                        v(4) = y(2);
                        v(5) = x(2);
                    }
                }
                else
                {
                    v(3) = y(2);
                    v(4) = x(1);
                    v(5) = x(2);
                }
            }
        }
    }
    else
    {
        v(0) = y(0);
        if (x(0) < y(1))
        {
            v(1) = x(0);
            if (x(1) < y(1))
            {
                v(2) = x(1);
                if (x(2) < y(1))
                {
                    v(3) = x(2);
                    v(4) = y(1);
                    v(5) = y(2);
                }
                else
                {
                    v(3) = y(1);
                    if (x(2) < y(2))
                    {
                        v(4) = x(2);
                        v(5) = y(2);
                    }
                    else
                    {
                        v(4) = y(2);
                        v(5) = x(2);
                    }
                }
            }
            else
            {
                v(2) = y(1);
                if (x(1) < y(2))
                {
                    v(3) = x(1);
                    if (x(2) < y(2))
                    {
                        v(4) = x(2);
                        v(5) = y(2);
                    }
                    else
                    {
                        v(4) = y(2);
                        v(5) = x(2);
                    }
                }
                else
                {
                    v(3) = y(2);
                    v(4) = x(1);
                    v(5) = x(2);
                }
            }
        }
        else
        {
            v(1) = y(1);
            if (x(0) < y(2))
            {
                v(2) = x(0);
                if (x(1) < y(2))
                {
                    v(3) = x(1);
                    if (x(2) < y(2))
                    {
                        v(4) = x(2);
                        v(5) = y(2);
                    }
                    else
                    {
                        v(4) = y(2);
                        v(5) = x(2);
                    }
                }
                else
                {
                    v(3) = y(2);
                    v(4) = x(1);
                    v(5) = x(2);
                }
            }
            else
            {
                v(2) = y(2);
                v(3) = x(0);
                v(4) = x(1);
                v(5) = x(2);
            }
        }
    }
}

// TODO Document
template <typename T>
void median(Matrix1D< T >& x, Matrix1D< T >& y, T& m)
{
    if (x(0) < y(1))
        if (x(1) < y(1))
            if (x(2) < y(1))
                m = y(1);
            else
                m = x(2);
        else if (x(1) < y(2))
            m = y(2);
        else
            m = x(1);
    else if (x(0) < y(2))
        if (x(1) < y(2))
            m = y(2);
        else
            m = x(1);
    else if (x(0) < y(3))
        m = y(3);
    else if (x(0) < y(4))
        m = x(0);
    else
        m = y(4);
}

/** Median_filter with a 3x3 window
 * @ingroup Filters
 */
template <typename T>
void median_filter3x3(Matrix2D< T >&m, Matrix2D< T >& out)
{
    int backup_startingx = STARTINGX(m);
    int backup_startingy = STARTINGY(m);

    STARTINGX(m) = STARTINGY(m) = 0;
    Matrix1D< T > v1(3), v2(3), v3(3), v4(3);
    Matrix1D< T > v(6);

    // Set the output matrix size
    out.resize(m);

    // Set the initial and final matrix indices to explore
    int initialY = 1, initialX = 1;
    int finalY = m.rowNumber() - 2;
    int finalX = m.colNumber() - 2;

    // For every row
    for (int i = initialY; i <= finalY; i++)
    {
        // For every pair of pixels (mean is computed obtaining
        // two means at the same time using an efficient method)
        for (int j = initialX; j <= finalX; j += 2)
        {
            // If we are in the border case
            if (j == 1)
            {
                // Order the first and second vectors of 3 elements
                sort(DIRECT_MAT_ELEM(m, i - 1, j - 1), DIRECT_MAT_ELEM(m, i, j - 1),
                     DIRECT_MAT_ELEM(m, i + 1, j - 1), v1);
                sort(DIRECT_MAT_ELEM(m, i - 1, j), DIRECT_MAT_ELEM(m, i, j),
                     DIRECT_MAT_ELEM(m, i + 1, j), v2);
            }
            else
            {
                // Simply take ordered vectors from previous
                v1 = v3;
                v2 = v4;
            }

            // As we are computing 2 medians at the same time, if the matrix has
            // an odd number of columns, the last column isn't calculated. It is
            // done here
            if (j == finalX)
            {
                v1 = v3;
                v2 = v4;
                sort(DIRECT_MAT_ELEM(m, i - 1, j + 1), DIRECT_MAT_ELEM(m, i, j + 1),
                     DIRECT_MAT_ELEM(m, i + 1, j + 1), v3);
                fast_merge_sort(v2, v3, v);
                median(v1, v, out(i, j));
            }
            else
            {
                // Order the third and fourth vectors of 3 elements
                sort(DIRECT_MAT_ELEM(m, i - 1, j + 1), DIRECT_MAT_ELEM(m, i, j + 1),
                     DIRECT_MAT_ELEM(m, i + 1, j + 1), v3);
                sort(DIRECT_MAT_ELEM(m, i - 1, j + 2), DIRECT_MAT_ELEM(m, i, j + 2),
                     DIRECT_MAT_ELEM(m, i + 1, j + 2), v4);

                // Merge sort the second and third vectors
                fast_merge_sort(v2, v3, v);

                // Find the first median and assign it to the output
                median(v1, v, out(i, j));

                // Find the second median and assign it to the output
                median(v4, v, out(i, j + 1));
            }
        }
    }

    STARTINGX(m) = STARTINGX(out) = backup_startingx;
    STARTINGY(m) = STARTINGY(out) = backup_startingy;
}

/** Mumford-Shah smoothing
 * @ingroup Filters
 *
 * This function simultaneously smooths and segments an image using non-linear
 * diffusion. Mumford-&-Shah's functional minimization algorithm is used to
 * detect region boundaries and relax image smoothness constraints near these
 * discontinuities. The functional minimized is:
 *
 *   E = W0*(f-d)*(f-d)               (data matching)
 * + W1*(fx*fx + fy*fy)*(1-s)*(1-s)      (1st deriv smooth)
 * + W2*(s*s)          (edge strengths)
 * + W3*(sx*sx + sy*sy)        (edge smoothness)
 *
 * The program diffusion from KUIM (developed by J. Gaush, U. Kansas) was used
 * as the "seed".
 *
 * Paper: Teboul, et al. IEEE-Trans. on Image Proc. Vol. 7, 387-397.
 */
void Smoothing_Shah(Matrix2D< double >& img,
                    Matrix2D< double >& surface_strength,
                    Matrix2D< double >& edge_strength,
                    const Matrix1D< double >& W,
                    int OuterLoops,
                    int InnerLoops,
                    int RefinementLoops,
                    bool adjust_range = true);

/** Tomographic diffusion
 * @ingroup Filters
 *
 * The direction of the tilt axis must be taken into account in the
 * definition of the diffusion constants alpha.
 *
 * The function returns the value of the regularization term.
 */
double tomographicDiffusion(Matrix3D< double >& V,
                            const Matrix1D< double >& alpha, double lambda);

/** Rotational invariant moments
 * @ingroup Filters
 *
 * The mask and the image are supposed to be of the same shape. If no mask is
 * provided, the moments are computed on the whole image. The moments are
 * measured with respect to the origin of the image.
 *
 * These moments have been taken from
 * http://www.cs.cf.ac.uk/Dave/Vision_lecture/node36.html (moments 1 to 5).
 */
void rotational_invariant_moments(const Matrix2D< double >& img,
                                  const Matrix2D< int >* mask,
                                  Matrix1D< double >& v_out);

/** Inertia moments
 * @ingroup Filters
 *
 * They are measured with respect to the center of the image, and not with
 * respect to the center of mass. For an image there are only two inertia
 * moments. v_out contains the inertia moments while the columns of u contain
 * the directions of the principal axes.
 */
void inertia_moments(const Matrix2D< double >& img,
                     const Matrix2D< int >* mask,
                     Matrix1D< double >& v_out,
                     Matrix2D< double >& u);

/** Fill a triangle defined by three points
 * @ingroup Filters
 *
 * The points are supplied as a pointer to three integer positions. They can be
 * negative
 */
void fill_triangle(Matrix2D< double >&img, int* tx, int* ty, double color);

/** Local thresholding
 * @ingroup Filters
 *
 * This function implements the local thresholding as described in
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/adpthrsh.htm
 * A mask can be supplied to limit the image area.
 *
 * The procedure employed is the following:
 * - Convolving the image with the statistical operator,
 * i.e. the mean or median (the size of the convolution kernel is dimLocal)
 * - Subtracting the original from the convolved image.
 * - Thresholding the difference image with C.
 * - Inverting the thresholded image.
 */
void local_thresholding(Matrix2D< double >& img,
                        double C,
                        double dimLocal,
                        Matrix2D< int >& result,
                        Matrix2D< int >* mask = NULL);

/** Center an image translationally
 * @ingroup Filters
 *
 * Given an image, this function returns the image that has been centered
 * translationally. For doing so, it compares this image with its mirrored
 * (X, Y, XY) versions.
 */
void centerImageTranslationally(Matrix2D<double> &I);

/** Center an image rotationally
 * @ingroup Filters
 *
 * Given an image, this function returns the image that has been centered
 * rotationally. For doing so, it compares this image with its mirrored
 * (X) version.
 */
void centerImageRotationally(Matrix2D<double> &I);

/** Center an image both translationally and rotationally
 * @ingroup Filters
 *
 * Given an image, this function returns the image that has been centered
 * rotationally and translationally. For doing so, it compares this image
 * with its mirrored (X, Y, XY) versions. The image is aligned translationally
 * and then rotationally Niter times.
 */
void centerImage(Matrix2D<double> &I, int Niter=10, bool limitShift=true);

#endif
