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

#include <queue>
#include "xmipp_image.h"
#include "numerical_tools.h"
#include "histogram.h"
#include "xmipp_program.h"
#include "mask.h"
#include "polar.h"

/// @defgroup Filters Filters
/// @ingroup DataLibrary

/** Substract background
 * @ingroup Filters
 *
 * The background is computed as the plane which best fits all density values,
 * then this plane is substracted from the image.
 */
void substractBackgroundPlane(MultidimArray<double> &I);

/** Substract background
 * @ingroup Filters
 *
 * The background is computed as a rolling ball operation with a ball
 * with this radius. The radius is typically Xdim/10.
 *
 * This code has been implemented after the one of "Subtract background" in
 * ImageJ.
 */
void substractBackgroundRollingBall(MultidimArray<double> &I, int radius);

/** Detect background
 * @ingroup Filters
 * This function receives a Matrix3D vol, and try to find the background
 * assuming that all the outside planes contain background, and apply
 * interval confidence, were alpha is the probabity to fail.
 * Mask is of the same size of vol, and is the solution, mask have
 * value 1 if background else value 0
*/
void detectBackground(const MultidimArray<double> &vol, MultidimArray<double> &mask, double alpha,
                      double &final_mean);

/** Constrast enhancement
 * @ingroup Filters
 *
 * The minimum density value is brought to 0 and the maximum to 255.
 */
void contrastEnhancement(Image<double>* I);

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
void regionGrowing2D(const MultidimArray< double >& I_in,
                     MultidimArray< double >& I_out,
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
void regionGrowing3D(const MultidimArray< double >& V_in,
                     MultidimArray< double >& V_out,
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
void distanceTransform(const MultidimArray<int> &in,
                       MultidimArray<int> &out, bool wrap=false);

/** Label a binary image
 * @ingroup Filters
 *
 * This function receives a binary volume and labels all its connected
 * components. The background is labeled as 0, and the components as 1, 2, 3
 * ...
 */
int labelImage2D(const MultidimArray< double >& I,
                 MultidimArray< double >& label,
                 int neighbourhood = 8);

/** Label a binary volume
 * @ingroup Filters
 *
 * This function receives a binary image and labels all its connected
 * components. The background is labeled as 0, and the components as 1, 2, 3
 * ...
 */
int labelImage3D(const MultidimArray< double >& V, MultidimArray< double >& label);

/** Remove connected components
 * @ingroup Filters
 *
 * Remove connected components smaller than a given size. They are set to 0.
 */
void removeSmallComponents(MultidimArray< double >& I,
                           int size,
                           int neighbourhood = 8);

/** Keep the biggest connected component
 * @ingroup Filters
 *
 * If the biggest component does not cover the percentage required (by default,
 * 0), more big components are taken until this is accomplished.
 */
void keepBiggestComponent(MultidimArray< double >& I,
                          double percentage = 0,
                          int neighbourhood = 8);

/** Fill object
 * @ingroup Filters
 *
 * Everything that is not background is assumed to be object.
 */
void fillBinaryObject(MultidimArray< double >&I, int neighbourhood = 8);

/** Segment an object using Otsu's method
 * @ingroup Filters
 *
 * Otsu's method determines a threshold such that the variance of the
 * two classes is minimized
 *
 * http://www.biomecardio.com/matlab/otsu.html
 */
void OtsuSegmentation(MultidimArray<double> &V);

/** Segment an object using Entropy method
 * @ingroup Filters
 *
 * Entropy method determines a threshold such that the entropy of the
 * two classes is maximized
 *
 * http://rsbweb.nih.gov/ij/plugins/download/Entropy_Threshold.java
 */
void EntropySegmentation(MultidimArray<double> &V);

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
 *
 * The binarization threshold is returned
 */
double EntropyOtsuSegmentation(MultidimArray<double> &V, double percentil=0.05,
                               bool binarizeVolume=true);

/** Correlation nD
 * @ingroup Filters
 */
template <typename T>
double correlation(const MultidimArray< T >& x,
                   const MultidimArray< T >& y,
                   const MultidimArray< int >* mask = NULL,
                   int l = 0,
                   int m = 0,
                   int q = 0)
{
    SPEED_UP_temps;

    double retval = 0; // returned value
    int i, j, k, ip, jp, kp; // indexes
    int Rows, Cols, Slices; // of the volumes

    // do the computation
    Cols = XSIZE(x);
    Rows = YSIZE(x);
    Slices = ZSIZE(x);

    long N = 0;
    for (k = 0; k < Slices; k++)
    {
        kp = k - q;
        if (kp < 0 || kp >= Slices)
            continue;
        for (i = 0; i < Rows; i++)
        {
            ip = i - l;
            if (ip < 0 || ip >= Rows)
                continue;
            for (j = 0; j < Cols; j++)
            {
                jp = j - m;

                if (jp >= 0 && jp < Cols)
                {
                    if (mask != NULL)
                        if (!DIRECT_A3D_ELEM((*mask), k, i, j))
                            continue;

                    retval += DIRECT_A3D_ELEM(x, k, i, j) *
                              DIRECT_A3D_ELEM(y, kp, ip, jp);
                    ++N;
                }
            }
        }
    }

    return retval / N;
}

/** Correlation nD
 * @ingroup Filters
 */
template <typename T>
double fastMaskedCorrelation(const MultidimArray< T >& x,
                             const MultidimArray< T >& y,
                             const MultidimArray< int >& mask)
{
    double retval = 0; // returned value
    long N = 0;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(x)
    {
        if (DIRECT_MULTIDIM_ELEM(mask,n))
        {
            retval += DIRECT_MULTIDIM_ELEM(x, n) * DIRECT_MULTIDIM_ELEM(y, n);
            ++N;
        }
    }
    return retval / N;
}

template <typename T>
double fastCorrelation(const MultidimArray< T >& x,
                       const MultidimArray< T >& y)
{
    double retval = 0; // returned value
    size_t nmax=4*(MULTIDIM_SIZE(x)/4);
    // loop unrolling
    for (size_t n=0; n<nmax; n+=4)
    {
        size_t n_1=n+1;
        size_t n_2=n+2;
        size_t n_3=n+3;
        retval += DIRECT_MULTIDIM_ELEM(x, n)   * DIRECT_MULTIDIM_ELEM(y, n)+
                  DIRECT_MULTIDIM_ELEM(x, n_1) * DIRECT_MULTIDIM_ELEM(y, n_1)+
                  DIRECT_MULTIDIM_ELEM(x, n_2) * DIRECT_MULTIDIM_ELEM(y, n_2)+
                  DIRECT_MULTIDIM_ELEM(x, n_3) * DIRECT_MULTIDIM_ELEM(y, n_3);
    }
    for (size_t n=nmax; n<MULTIDIM_SIZE(x); ++n)
        retval += DIRECT_MULTIDIM_ELEM(x, n)   * DIRECT_MULTIDIM_ELEM(y, n);
    return retval / MULTIDIM_SIZE(x);
}

/** correlationIndex nD
 * @ingroup Filters
 */
template <typename T>
double correlationIndex(const MultidimArray< T >& x,
                        const MultidimArray< T >& y,
                        const MultidimArray< int >* mask = NULL,
                        MultidimArray< double >* Contributions = NULL)
{
	SPEED_UP_tempsInt;

    double retval = 0, aux;
    double mean_x, mean_y;
    double stddev_x, stddev_y;

    long N = 0;

    if (mask == NULL)
    {
        x.computeAvgStdev(mean_x, stddev_x);
        y.computeAvgStdev(mean_y, stddev_y);
    }
    else
    {
        x.computeAvgStdev_within_binary_mask(*mask, mean_x,stddev_x);
        y.computeAvgStdev_within_binary_mask(*mask, mean_y,stddev_y);
    }
    if (ABS(stddev_x)<XMIPP_EQUAL_ACCURACY ||
        ABS(stddev_y)<XMIPP_EQUAL_ACCURACY)
        return 0;

    // If contributions are desired. Please, be careful interpreting individual
    // contributions to the covariance! One pixel value afect others.
    if (Contributions != NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(x, y)
        {
            if (mask != NULL)
                if (!A3D_ELEM(*mask,k, i, j))
                    continue;

            aux = (A3D_ELEM(x, k, i, j) - mean_x) * (A3D_ELEM(y, k, i, j) -
                    mean_y);
            A3D_ELEM(*Contributions, k, i, j) = aux;
            retval += aux;
            ++N;
        }

        FOR_ALL_ELEMENTS_IN_ARRAY3D(*Contributions)
        A3D_ELEM(*Contributions, k, i, j) /= ((stddev_x * stddev_y) * N);
    }
    else
    {
        if (mask==NULL && x.sameShape(y))
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(x)
                retval += DIRECT_MULTIDIM_ELEM(x, n)*DIRECT_MULTIDIM_ELEM(y, n);
            N=MULTIDIM_SIZE(x);
            retval-=N*mean_x*mean_y;
        }
        else
        {
            FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(x, y)
            {
                if (mask != NULL)
                    if (!A3D_ELEM(*mask,k, i, j))
                        continue;

                retval += (A3D_ELEM(x, k, i, j) - mean_x) *
                          (A3D_ELEM(y, k, i, j) - mean_y);
                ++N;
            }
        }
    }

    if (N != 0)
        return retval / ((stddev_x * stddev_y) * N);
    else
        return 0;
}

/** Fast Correntropy 1D
 * @ingroup Filters
 */
template <typename T>
double fastCorrentropy(const MultidimArray<T> &x, const MultidimArray<T> &y,
                       double sigma, const GaussianInterpolator &G)
{
    double retval=0;
    double isigma=1.0/sigma;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(x)
    retval+=G.getValue(isigma*(DIRECT_MULTIDIM_ELEM(x,n)-
                               DIRECT_MULTIDIM_ELEM(y,n)));
    retval/=XSIZE(x);
    return retval;
}

/** Correntropy with mask. */
double fastCorrentropy(const MultidimArray<double> &x, const MultidimArray<double> &y,
                       double sigma, const GaussianInterpolator &G, const MultidimArray<int> &mask);

/** Correntropy nD
 * @ingroup Filters
 */
template <typename T>
double correntropy(const MultidimArray<T> &x, const MultidimArray<T> &y,
                   double sigma)
{
    double retval=0;
    double K=-0.5/(sigma*sigma);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(x)
    {
        double diff=DIRECT_MULTIDIM_ELEM(x,n)-DIRECT_MULTIDIM_ELEM(y,n);
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
 * I2 by (shiftX, shiftY).
 *
 * You can limit the maximum achievable shift by using maxShift. If it is set to -1,
 * any shift is valid.
 *
 * The function returns the maximum correlation found.
 */
double bestShift(const MultidimArray< double >& I1,
               const MultidimArray< double >& I2,
               double& shiftX,
               double& shiftY,
               CorrelationAux &aux,
               const MultidimArray< int >* mask = NULL,
               int maxShift=-1);

/** Translational search (3D)
 * @ingroup Filters
 *
 * This function returns the best interpolated shift for the alignment of two
 * volumes. You can restrict the shift to a region defined by a mask (the maximum
 * will be sought where the mask is 1).
 *
 * To apply these results you must shift I1 by (-shiftX,-shiftY,-shiftZ) or
 * I2 by (shiftX, shiftY,shiftZ).
 */
void bestShift(const MultidimArray<double> &I1, const MultidimArray<double> &I2,
               double &shiftX, double &shiftY, double &shiftZ, CorrelationAux &aux,
               const MultidimArray<int> *mask=NULL);

/** Translational search (non-wrapping)
 * @ingroup Filters
 *
 * This function returns the best interpolated shift for the alignment of two
 * images. You can restrict the shift to a region defined by a mask (the maximum
 * will be sought where the mask is 1).
 */
void bestNonwrappingShift(const MultidimArray< double >& I1,
                          const MultidimArray< double >& I2,
                          double& shiftX,
                          double& shiftY,
                          CorrelationAux &aux);

/** Translational search (non-wrapping).
 * @ingroup Filters
 *
 * Search is performed in real-space
 */
double bestShiftRealSpace(const MultidimArray<double> &I1, MultidimArray<double> &I2,
               double &shiftX, double &shiftY,
               const MultidimArray<int> *mask=NULL, int maxShift=5, double shiftStep=1.0);

/** Auxiliary class for fast image alignment */
class AlignmentAux
{
public:
    Matrix2D<double> ARS, ASR, R;
    MultidimArray<double> IauxSR, IauxRS, rotationalCorr;
    Polar_fftw_plans *plans;
    Polar< std::complex<double> > polarFourierIref, polarFourierI;
    AlignmentAux();
    ~AlignmentAux();
};

/** Align two images
 * @ingroup Filters
 *
 * This function modifies I to be aligned with Iref. Translational and
 * rotational alignments are both considered. The matrix transforming I
 * into Iref is returned.
 *
 * The function returns the correlation between the two aligned images.
 */
double alignImages(const MultidimArray< double >& Iref,
                   MultidimArray< double >& I,
                   Matrix2D< double >&M,
                   bool wrap=WRAP);

/** Align two images considering mirrors */
double alignImagesConsideringMirrors(const MultidimArray<double>& Iref, MultidimArray<double>& I,
                   Matrix2D<double>&M, bool wrap);

/** Fast version of align two images
 * @ingroup Filters
 */
double alignImages(const MultidimArray< double >& Iref,
                   MultidimArray< double >& I,
                   Matrix2D< double >&M,
                   bool wrap,
                   AlignmentAux &aux,
                   CorrelationAux &aux2,
                   RotationalCorrelationAux &aux3);

/** Auxiliary class for fast volume alignment */
class VolumeAlignmentAux
{
public:
    MultidimArray<double> IrefCyl, Icyl, corr, I1, I12, I123;
    Matrix2D<double> R1, R2, R3;
};

/** Align two volumes by applying a rotation around Z.
 * The rotation must be applied on I to get Iref.
 * @ingroup Filters
 */
double bestRotationAroundZ(const MultidimArray< double >& Iref,
                           const MultidimArray< double >& I,
                           CorrelationAux &aux,
                           VolumeAlignmentAux &aux2);

/** Find a ZYZ rotation that transforms I into Iref.
 * The rotation matrix is returned in R. I is modified to be aligned.
 */
void fastBestRotation(const MultidimArray<double>& IrefCylZ,
                      const MultidimArray<double>& IrefCylY,
                      const MultidimArray<double>& IrefCylX,
                      MultidimArray<double>& I,
                      const String &eulerAngles,
                      Matrix2D<double> &R, CorrelationAux &aux, VolumeAlignmentAux &aux2);

/** Fast best rotation around Z */
double fastBestRotationAroundZ(const MultidimArray<double>& IrefCylZ,
                               const MultidimArray<double>& I, CorrelationAux &aux,
                               VolumeAlignmentAux &aux2);

/** Fast best rotation around Y */
double fastBestRotationAroundY(const MultidimArray<double>& IrefCylY,
                               const MultidimArray<double>& I, CorrelationAux &aux,
                               VolumeAlignmentAux &aux2);

/** Fast best rotation around X */
double fastBestRotationAroundX(const MultidimArray<double>& IrefCylX,
                               const MultidimArray<double>& I, CorrelationAux &aux,
                               VolumeAlignmentAux &aux2);

/** Align two images considering also the mirrors
 * @ingroup Filters
 *
 * This function modifies I to be aligned with Iref. Translational and
 * rotational alignments are both considered. The correlation coefficient
 * between I transformed and Iref is returned. A mask can be supplied for
 * computing this correlation.
 */
double alignImagesConsideringMirrors(const MultidimArray< double >& Iref,
                                     MultidimArray< double >& I,
                                     Matrix2D<double> &M,
                                     AlignmentAux &aux,
                                     CorrelationAux &aux2,
                                     RotationalCorrelationAux &aux3,
                                     bool wrap=WRAP,
                                     const MultidimArray< int >* mask = NULL);

/** Align a set of images.
 * Align a set of images and produce a class average as well as the set of
 * alignment parameters. The output is in Iavg. The metadata is modified by adding
 * the alignment information (transformation and euclidean distance to the
 * average. The process is iterative, first an average is computed. All
 * images are aligned to the average, and this is updated. The process
 * is run for a given number of iterations.
 */
void alignSetOfImages(MetaData &MD, MultidimArray< double >& Iavg,
                      int Niter=10, bool considerMirror=true);

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
void estimateGaussian2D(const MultidimArray<double> &I,
                        double &a, double &b, Matrix1D<double> &mu, Matrix2D<double> &sigma,
                        bool estimateMu=true, int iterations=10);

/** euclidian distance nD
 * @ingroup Filters
 */
template <typename T>
double euclidianDistance(const MultidimArray< T >& x,
                         const MultidimArray< T >& y,
                         const MultidimArray< int >* mask = NULL)
{
    SPEED_UP_temps;

    double retval = 0;
    long n = 0;

    FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(x, y)
    {
        if (mask != NULL)
            if (!(*mask)(k, i, j))
                continue;

        retval += (A3D_ELEM(x, k, i, j) - A3D_ELEM(y, k, i, j)) *
                  (A3D_ELEM(x, k, i, j) - A3D_ELEM(y, k, i, j));
        n++;
    }

    if (n != 0)
        return sqrt(retval);
    else
        return 0;
}

/** mutual information nD
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
double mutualInformation(const MultidimArray< T >& x,
                         const MultidimArray< T >& y,
                         int nx = 0,
                         int ny = 0,
                         const MultidimArray< int >* mask = NULL)
{
    SPEED_UP_temps;

    long n = 0;
    Histogram1D histx, histy;
    Histogram2D histxy;
    MultidimArray< T > aux_x, aux_y;
    MultidimArray< double > mx, my;
    MultidimArray< double > mxy;
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

    FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(x, y)
    {
        if (mask != NULL)
            if (!(*mask)(k, i, j))
                continue;

        aux_x(n) = A3D_ELEM(x, k, i, j);
        aux_y(n) = A3D_ELEM(y, k, i, j);
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

/** RMS nD
 * @ingroup Filters
 */
template <typename T>
double rms(const MultidimArray< T >& x,
           const MultidimArray< T >& y,
           const MultidimArray< int >* mask = NULL,
           MultidimArray< double >* Contributions = NULL)
{
    SPEED_UP_tempsInt;

    double retval = 0;
    double aux;
    int n = 0;

    // If contributions are desired
    if (Contributions != NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(k, i, j))
                    continue;

            aux = (A3D_ELEM(x, k, i, j) - A3D_ELEM(y, k, i, j)) *
                  (A3D_ELEM(x, k, i, j) - A3D_ELEM(y, k, i, j));
            A3D_ELEM(*Contributions, k, i, j) = aux;
            retval += aux;
            n++;
        }

        FOR_ALL_ELEMENTS_IN_ARRAY3D(*Contributions)
        A3D_ELEM(*Contributions, k, i, j) = sqrt(A3D_ELEM(*Contributions,
                                            k, i, j) / n);
    }
    else
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(x, y)
        {
            if (mask != NULL)
                if (!(*mask)(k, i, j))
                    continue;

            retval += (A3D_ELEM(x, k, i, j) - A3D_ELEM(y, k, i, j)) *
                      (A3D_ELEM(x, k, i, j) - A3D_ELEM(y, k, i, j));
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
void fourierBesselDecomposition(const MultidimArray< double >& img_in,
                                MultidimArray< double >& m_out,
                                double r1,
                                double r2,
                                int k1,
                                int k2);

/** Harmonic decomposition
 * @ingroup Filters
 */
void harmonicDecomposition(const MultidimArray< double >& img_in,
                           MultidimArray< double >& v_out);

// Function needed by median filtering
template <typename T>
void sort(T a, T b, T c, MultidimArray< T >& v)
{
    if (a < b)
        if (b < c)
        {
            DIRECT_MULTIDIM_ELEM(v,0) = a;
            DIRECT_MULTIDIM_ELEM(v,1) = b;
            DIRECT_MULTIDIM_ELEM(v,2) = c;
        }
        else if (a < c)
        {
            DIRECT_MULTIDIM_ELEM(v,0) = a;
            DIRECT_MULTIDIM_ELEM(v,1) = c;
            DIRECT_MULTIDIM_ELEM(v,2) = b;
        }
        else
        {
            DIRECT_MULTIDIM_ELEM(v,0) = c;
            DIRECT_MULTIDIM_ELEM(v,1) = a;
            DIRECT_MULTIDIM_ELEM(v,2) = b;
        }
    else if (a < c)
    {
        DIRECT_MULTIDIM_ELEM(v,0) = b;
        DIRECT_MULTIDIM_ELEM(v,1) = a;
        DIRECT_MULTIDIM_ELEM(v,2) = c;
    }
    else if (b < c)
    {
        DIRECT_MULTIDIM_ELEM(v,0) = b;
        DIRECT_MULTIDIM_ELEM(v,1) = c;
        DIRECT_MULTIDIM_ELEM(v,2) = a;
    }
    else
    {
        DIRECT_MULTIDIM_ELEM(v,0) = c;
        DIRECT_MULTIDIM_ELEM(v,1) = b;
        DIRECT_MULTIDIM_ELEM(v,2) = a;
    }
}

// Function needed by median filtering
template <typename T>
void mergeSort(MultidimArray< T >& v1, MultidimArray< T >& v2, MultidimArray< T >& v)
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

// Function needed by median filtering
// This UGLY function performs a fast merge sort for the case of vectors of 3
// elements. This way is guaranteed a minimum number of comparisons (maximum
// number of comparisons to perform the sort, 5)
template <typename T>
void fastMergeSort(MultidimArray< T >& x, MultidimArray< T >& y, MultidimArray< T >& v)
{
    if (DIRECT_MULTIDIM_ELEM(x,0) < DIRECT_MULTIDIM_ELEM(y,0))
    {
        DIRECT_MULTIDIM_ELEM(v,0) = DIRECT_MULTIDIM_ELEM(x,0);
        if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,0))
        {
            DIRECT_MULTIDIM_ELEM(v,1) = DIRECT_MULTIDIM_ELEM(x,1);
            if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,0))
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(x,2);
                DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(y,0);
                DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,1);
                DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
            }
            else
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(y,0);
                if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,1))
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(x,2);
                    DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,1);
                    DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(y,1);
                    if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,2))
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                    }
                }
            }
        }
        else
        {
            DIRECT_MULTIDIM_ELEM(v,1) = DIRECT_MULTIDIM_ELEM(y,0);
            if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,1))
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(x,1);
                if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,1))
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(x,2);
                    DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,1);
                    DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(y,1);
                    if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,2))
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                    }
                }
            }
            else
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(y,1);
                if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,2))
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(x,1);
                    if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,2))
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                    }
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(y,2);
                    DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,1);
                    DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                }
            }
        }
    }
    else
    {
        DIRECT_MULTIDIM_ELEM(v,0) = DIRECT_MULTIDIM_ELEM(y,0);
        if (DIRECT_MULTIDIM_ELEM(x,0) < DIRECT_MULTIDIM_ELEM(y,1))
        {
            DIRECT_MULTIDIM_ELEM(v,1) = DIRECT_MULTIDIM_ELEM(x,0);
            if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,1))
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(x,1);
                if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,1))
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(x,2);
                    DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,1);
                    DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(y,1);
                    if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,2))
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                    }
                }
            }
            else
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(y,1);
                if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,2))
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(x,1);
                    if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,2))
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                    }
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(y,2);
                    DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,1);
                    DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                }
            }
        }
        else
        {
            DIRECT_MULTIDIM_ELEM(v,1) = DIRECT_MULTIDIM_ELEM(y,1);
            if (DIRECT_MULTIDIM_ELEM(x,0) < DIRECT_MULTIDIM_ELEM(y,2))
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(x,0);
                if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,2))
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(x,1);
                    if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,2))
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(y,2);
                    }
                    else
                    {
                        DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(y,2);
                        DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                    }
                }
                else
                {
                    DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(y,2);
                    DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,1);
                    DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
                }
            }
            else
            {
                DIRECT_MULTIDIM_ELEM(v,2) = DIRECT_MULTIDIM_ELEM(y,2);
                DIRECT_MULTIDIM_ELEM(v,3) = DIRECT_MULTIDIM_ELEM(x,0);
                DIRECT_MULTIDIM_ELEM(v,4) = DIRECT_MULTIDIM_ELEM(x,1);
                DIRECT_MULTIDIM_ELEM(v,5) = DIRECT_MULTIDIM_ELEM(x,2);
            }
        }
    }
}

// Function needed by median filtering
template <typename T>
void median(MultidimArray< T >& x, MultidimArray< T >& y, T& m)
{
    if (DIRECT_MULTIDIM_ELEM(x,0) < DIRECT_MULTIDIM_ELEM(y,1))
        if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,1))
            if (DIRECT_MULTIDIM_ELEM(x,2) < DIRECT_MULTIDIM_ELEM(y,1))
                m = DIRECT_MULTIDIM_ELEM(y,1);
            else
                m = DIRECT_MULTIDIM_ELEM(x,2);
        else
            if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,2))
                m = DIRECT_MULTIDIM_ELEM(y,2);
            else
                m = DIRECT_MULTIDIM_ELEM(x,1);
    else
        if (DIRECT_MULTIDIM_ELEM(x,0) < DIRECT_MULTIDIM_ELEM(y,2))
            if (DIRECT_MULTIDIM_ELEM(x,1) < DIRECT_MULTIDIM_ELEM(y,2))
                m = DIRECT_MULTIDIM_ELEM(y,2);
            else
                m = DIRECT_MULTIDIM_ELEM(x,1);
        else
            if (DIRECT_MULTIDIM_ELEM(x,0) < DIRECT_MULTIDIM_ELEM(y,3))
                m = DIRECT_MULTIDIM_ELEM(y,3);
            else
                if (DIRECT_MULTIDIM_ELEM(x,0) < DIRECT_MULTIDIM_ELEM(y,4))
                    m = DIRECT_MULTIDIM_ELEM(x,0);
                else
                    m = DIRECT_MULTIDIM_ELEM(y,4);
}

/** Median_filter with a 3x3 selfWindow
 * @ingroup Filters
 */
template <typename T>
void medianFilter3x3(MultidimArray< T >&m, MultidimArray< T >& out)
{
    int backup_startingx = STARTINGX(m);
    int backup_startingy = STARTINGY(m);

    STARTINGX(m) = STARTINGY(m) = 0;
    MultidimArray< T > v1(3), v2(3), v3(3), v4(3);
    MultidimArray< T > v(6);

    // Set the output matrix size
    out.initZeros(m);

    // Set the initial and final matrix indices to explore
    int initialY = 1, initialX = 1;
    int finalY = YSIZE(m) - 2;
    int finalX = XSIZE(m) - 2;

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
                sort(DIRECT_A2D_ELEM(m, i - 1, j - 1), DIRECT_A2D_ELEM(m, i, j - 1),
                     DIRECT_A2D_ELEM(m, i + 1, j - 1), v1);
                sort(DIRECT_A2D_ELEM(m, i - 1, j), DIRECT_A2D_ELEM(m, i, j),
                     DIRECT_A2D_ELEM(m, i + 1, j), v2);
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
                sort(DIRECT_A2D_ELEM(m, i - 1, j + 1), DIRECT_A2D_ELEM(m, i, j + 1),
                     DIRECT_A2D_ELEM(m, i + 1, j + 1), v3);
                fastMergeSort(v2, v3, v);
                median(v1, v, DIRECT_A2D_ELEM(out,i, j));
            }
            else
            {
                // Order the third and fourth vectors of 3 elements
                sort(DIRECT_A2D_ELEM(m, i - 1, j + 1), DIRECT_A2D_ELEM(m, i, j + 1),
                     DIRECT_A2D_ELEM(m, i + 1, j + 1), v3);
                sort(DIRECT_A2D_ELEM(m, i - 1, j + 2), DIRECT_A2D_ELEM(m, i, j + 2),
                     DIRECT_A2D_ELEM(m, i + 1, j + 2), v4);

                // Merge sort the second and third vectors
                fastMergeSort(v2, v3, v);

                // Find the first median and assign it to the output
                median(v1, v, DIRECT_A2D_ELEM(out, i, j));

                // Find the second median and assign it to the output
                median(v4, v, DIRECT_A2D_ELEM(out, i, j + 1));
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
void smoothingShah(MultidimArray< double >& img,
                   MultidimArray< double >& surface_strength,
                   MultidimArray< double >& edge_strength,
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
double tomographicDiffusion(MultidimArray< double >& V,
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
void rotationalInvariantMoments(const MultidimArray< double >& img,
                                const MultidimArray< int >* mask,
                                MultidimArray< double >& v_out);

/** Inertia moments
 * @ingroup Filters
 *
 * They are measured with respect to the center of the image, and not with
 * respect to the center of mass. For an image there are only two inertia
 * moments. v_out contains the inertia moments while the columns of u contain
 * the directions of the principal axes.
 */
void inertiaMoments(const MultidimArray< double >& img,
                    const MultidimArray< int >* mask,
                    Matrix1D< double >& v_out,
                    Matrix2D< double >& u);

/** Fill a triangle defined by three points
 * @ingroup Filters
 *
 * The points are supplied as a pointer to three integer positions. They can be
 * negative
 */
void fillTriangle(MultidimArray< double >&img, int* tx, int* ty, double color);

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
void localThresholding(MultidimArray< double >& img,
                       double C,
                       double dimLocal,
                       MultidimArray< int >& result,
                       MultidimArray< int >* mask = NULL);

/** Center an image translationally
 * @ingroup Filters
 *
 * Given an image, this function returns the image that has been centered
 * translationally. For doing so, it compares this image with its mirrored
 * (X, Y, XY) versions.
 */
void centerImageTranslationally(MultidimArray<double> &I,
                                CorrelationAux &aux);

/** Center an image rotationally
 * @ingroup Filters
 *
 * Given an image, this function returns the image that has been centered
 * rotationally. For doing so, it compares this image with its mirrored
 * (X) version.
 */
void centerImageRotationally(MultidimArray<double> &I, RotationalCorrelationAux &aux);

/** Center an image both translationally and rotationally
 * @ingroup Filters
 *
 * Given an image, this function returns the image that has been centered
 * rotationally and translationally. For doing so, it compares this image
 * with its mirrored (X, Y, XY) versions. The image is aligned translationally
 * and then rotationally Niter times.
 */
void centerImage(MultidimArray<double> &I, CorrelationAux &aux,
                 RotationalCorrelationAux &aux2,
                 int Niter=10, bool limitShift=true);


/** Force positive.
 *  * @ingroup Filters
 *
 *  A median filter is applied at those negative values. Positive values are untouched.
 */
void forcePositive(MultidimArray<double> &V);


/** Remove bad pixels.
 *  * @ingroup Filters
 *
 *  A boundaries median filter is applied at those pixels given by the mask.
 */
template <typename T>
void boundMedianFilter(MultidimArray< T > &V, const MultidimArray<char> &mask, int n=0)
{
    bool badRemaining;
    T neighbours[125];
    T aux;
    int N = 0, index;

    do
    {
        badRemaining=false;

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V)
        if (DIRECT_A3D_ELEM(mask, k, i, j) != 0)
        {
            N = 0;
            for (int kk=-2; kk<=2; kk++)
            {
            	size_t kkk=k+kk;
                if (kkk<0 || kkk>=ZSIZE(V))
                    continue;
                for (int ii=-2; ii<=2; ii++)
                {
                	size_t iii=i+ii;
                    if (iii<0 || iii>=YSIZE(V))
                        continue;
                    for (int jj=-2; jj<=2; jj++)
                    {
                    	size_t jjj=j+jj;
                        if (jjj<0 || jjj>=XSIZE(V))
                            continue;

                        if (DIRECT_A3D_ELEM(mask, kkk, iii, jjj) == 0)
                        {
                            index = N++;
                            neighbours[index] = DIRECT_A3D_ELEM(V, kkk,iii,jjj);
                            //insertion sort
                            while (index > 0 && neighbours[index-1] > neighbours[index])
                            {
                                SWAP(neighbours[index-1], neighbours[index], aux);
                                --index;
                            }
                        }
                    }
                }
                if (N == 0)
                    badRemaining = true;
                else
                {
                    //std::sort(neighbours.begin(),neighbours.end());
                    if (N % 2 == 0)
                        DIRECT_A3D_ELEM(V, k, i, j) = (T)(0.5*(neighbours[N/2-1]+ neighbours[N/2]));
                    else
                        DIRECT_A3D_ELEM(V, k, i, j) = neighbours[N/2];
                    DIRECT_A3D_ELEM(mask, k, i, j) = false;
                }
            }
        }
    }
    while (badRemaining);
}

/** Remove bad pixels.
 *  * @ingroup Filters
 *
 *  A boundaries median filter is applied at those pixels whose value is out of range
 *  given by thresFactor * std.
  */
template <typename T>
void pixelDesvFilter(MultidimArray< T > &V, double thresFactor)
{
    if (thresFactor > 0 )
    {
        double avg, stddev, high, low;
        T dummy;
        MultidimArray<char> mask(ZSIZE(V), YSIZE(V), XSIZE(V));
        avg = stddev = low = high = 0;
        V.computeStats(avg, stddev, dummy, dummy);//min and max not used
        low  = (avg - thresFactor * stddev);
        high = (avg + thresFactor * stddev);

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
        {
            double x = DIRECT_MULTIDIM_ELEM(V, n);
            DIRECT_MULTIDIM_ELEM(mask, n) = (x < low || x > high) ? 1 : 0;
        }
        boundMedianFilter(V, mask);
    }
}

/** Compute logarithm.
 *  * @ingroup Filters
 *
 *  apply transformation a+b*ln(x+c). i.e:  4.431-0.4018*LN(P1+336.6)
 *
 */
//#include <math.h>

template <typename T>
void logFilter(MultidimArray< T > &V, double a, double b, double c)
{

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
    {
        double x = DIRECT_MULTIDIM_ELEM(V, n);
        //this casting is kind of risky
        DIRECT_MULTIDIM_ELEM(V, n) = (T)(a-b*log(x+c));
    }
}

/** Compute edges with Sobel */
void computeEdges (const MultidimArray <double>& vol, MultidimArray<double> &vol_edge);

/** Force sparsity.
 * Only eps*100  % of the DWT coefficients are kept. It is assumed that the input volume
 * is cubic.
 */
void forceDWTSparsity(MultidimArray<double> &V, double eps);

/** Abstract class that will be the base for all filters */
class XmippFilter
{
public:
    /** Read params from a program */
    virtual void readParams(XmippProgram *program)
    {}//do nothing by default
    /** Apply filter to an image */
    virtual void apply(MultidimArray<double> &img) = 0;

    /** Virtual destructor */
    virtual ~XmippFilter()
    {}
    ;

    /** Show some info before running */
    virtual void show()
    {}
    ;
};

/** Some concrete filters */


class BadPixelFilter: public XmippFilter
{
public:
    /** Apply filter on bad pixels */
    typedef enum { NEGATIVE, MASK, OUTLIER } BadPixelFilterType;

    BadPixelFilterType type; //type of filter
    double factor;    //for the case of outliers bad pixels
    Image<char> *mask; //for the case of mask bad pixels

    /** Define the parameters for use inside an Xmipp program */
    static void defineParams(XmippProgram * program);
    /** Read from program command line */
    void readParams(XmippProgram * program);
    /** Apply the filter to an image or volume*/
    void apply(MultidimArray<double> &img);
};

class LogFilter: public XmippFilter
{
public:
    /** Apply filter on bad pixels */
    //4.431-0.4018*LN(ABS(P1+336.6))
    //a-b*ln(x+c)
    double a,b,c;
    /** Define the parameters for use inside an Xmipp program */
    static void defineParams(XmippProgram * program);
    /** Read from program command line */
    void readParams(XmippProgram * program);
    /** Apply the filter to an image or volume*/
    void apply(MultidimArray<double> &img);
};

class BackgroundFilter: public XmippFilter
{
    /** Apply filter on bad pixels */
    typedef enum { PLANE, ROLLINGBALL } BackgroundType;
    BackgroundType type;
    int radius;

public:
    /** Define the parameters for use inside an Xmipp program */
    static void defineParams(XmippProgram * program);
    /** Read from program command line */
    void readParams(XmippProgram * program);
    /** Apply the filter to an image or volume*/
    void apply(MultidimArray<double> &img);
};

class MedianFilter: public XmippFilter
{
public:
    /** Define the parameters for use inside an Xmipp program */
    static void defineParams(XmippProgram * program);
    /** Read from program command line */
    void readParams(XmippProgram * program);
    /** Apply the filter to an image or volume*/
    void apply(MultidimArray<double> &img);
};

class DiffusionFilter: public XmippFilter
{
    /** Shah number of outer iterations
     */
    int Shah_outer;

    /** Shah number of inner iterations
     */
    int Shah_inner;

    /** Shah number of refinement iterations
     */
    int Shah_refinement;

    /** Shah weight.
     *
     * w0 = data matching (=0)
     * w1 = 1st derivative smooth (=50)
     * w2 = edge strength (=50)
     * w3 = edge smoothness (=0.02)
     */
    Matrix1D< double > Shah_weight;

    /** Produce Shah edge instead of Shah smooth.
     */
    bool Shah_edge;
public:
    /** Define the parameters for use inside an Xmipp program */
    static void defineParams(XmippProgram * program);
    /** Read from program command line */
    void readParams(XmippProgram * program);
    /** Show parameters */
    void show();
    /** Apply the filter to an image or volume*/
    void apply(MultidimArray<double> &img);
};

class BasisFilter: public XmippFilter
{
    /** Stack with the basis */
    FileName fnBasis;

    /// Number of basis to use.
    int Nbasis;

    // Stack with the basis
    Image<double> basis;
public:
    /** Define the parameters for use inside an Xmipp program */
    static void defineParams(XmippProgram * program);
    /** Read from program command line */
    void readParams(XmippProgram * program);
    /** Show parameters */
    void show();
    /** Apply the filter to an image or volume*/
    void apply(MultidimArray<double> &img);
};


#endif
