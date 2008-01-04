/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "matrix1d.h"
#include "matrix2d.h"
#include "matrix3d.h"

/// @defgroup Histograms Histograms
/// @ingroup DataLibrary

/** Histograms with 1 parameter
 * @ingroup Histograms
 *
 * This class of histograms are the usual ones where we want to count the number
 * of elements within a certain range of a variable, then we make the histogram
 * of that variable. The range is divided into small subranges within which the
 * values will be grouped. Any value outside the global range will not be
 * counted in the histogram.
 *
 * To see exactly which is the division between subranges let's have a look on
 * the following example where an histogram between 0 and 2 is computed with
 * 5 steps.
 *
 * @code
 * [......)                        [.......]
 * [      )                        [       ]
 * [      )[......)                [       ]
 * [      )[      )        [......)[       ]
 * [      )[      )[......)[      )[       ]
 * [      )[      )[      )[      )[       ]
 * [   0  )[   1  )[   2  )[   3  )[    4  ]
 * [      )[      )[      )[      )[       ]
 * |---------|---------|---------|---------|
 * 0.0       0.5       1.0       1.5       2.0
 * @endcode
 *
 * The border points are 0.0, 0.4, 0.8, 1.2, 1.6 and 2.0. The brackets and
 * parenthesis try to represent where the border point belongs to, and the
 * numbers within the bars are the index of each bar whithin the histogram. The
 * height of each bar is the number of times that a value within that subrange
 * has been inserted. Be careful that this is not a probability density function
 * (pdf), to be so it should be divided by the total number of values inserted.
 *
 * The following example shows how to work with the histograms. In it we will
 * compute which is the central range within which the 95% of the values of a
 * matrix are comprised. This example could be even simplified by using the
 * function compute_hist but it has been kept like this to show the idea behind
 * the histograms
 *
 * @code
 * // Variable definition
 * histogram1D hist;
 * Matrix2D A(50, 50);
 * double min_val, max_val;
 * double eff0, effF;
 *
 * // Matrix initialisation
 * A.init_random(0, 100);
 *
 * // Histogram initialisation
 * min_val = A.min();
 * max_val = A.max();
 * hist.init(min_val, max_val, 200);
 *
 * // Histogram computation
 * for (int i=STARTINGY(A); i<=FINISHINGY(A); i++)
 *     for (int j=STARTINGX(A); j<=FINISHINGX(A); j++)
 *         hist.insert_value(MAT_ELEM(A,i,j));
 *
 * // Effective range computation
 * eff0 = hist.percentil(2.5);
 * effF = hist.percentil(97.5);
 *
 * std::cout << "The effective range goes from " << eff0
 *     << " to " << effF << std::endl;
 * @endcode
 */
class histogram1D: public Matrix1D< double >
{
public:
    // Structure
    double hmin; // minimum value of the histogram
    double hmax; // maximum value of the histogram
    double step_size; // size of step
    int no_samples; // No. of points inside the histogram

    /** Empty constructor
     *
     * Creates an empty histogram. Before using it you must
     * initialise it with init.
     *
     * @code
     * histogram1D hist;
     * @endcode
     */
    histogram1D()
    {
        clear();
    }

    /** Copy constructor
     *
     * Makes an exact copy of the given histogram into another histogram.
     *
     * @code
     * histogram1D hist2(hist1);
     * @endcode
     */
    histogram1D(const histogram1D& H)
    {
        clear();
        *this = H;
    }

    /** Empties an histogram
     *
     * Remember to initialise it before using it again.
     *
     * @code
     * hist.clear();
     * @endcode
     */
    void clear();

    /** Assignment */
    histogram1D& operator=(const histogram1D& H);

    /** Another function for assignament.*/
    void assign(const histogram1D& H);

    /** Initialisation of the histogram
     *
     * This is the operation which allows the histogram to be used. This should
     * be performed before inserting any value in it. The information given to
     * this initialisation is the range within which the values will be counted,
     * and the number of steps (discrete bars) in this range. If the value is
     * outside this range it will not be taken into account although we have
     * asked for its insertion in the histogram.
     *
     * @code
     * hist.init(-3, 3, 100);
     * // 100 steps in the range -3...3
     * @endcode
     */
    void init(double min_val, double max_val, int n_steps);

    /** Insert a value within histogram
     *
     * The right interval is chosen according to the initialisation of the
     * histogram and the count of elements in that interval is incremented by 1.
     * If the value lies outside the global range of the histogram nothing is
     * done.
     *
     * @code
     * hist.insert_value(3);
     * @endcode
     */
    void insert_value(double val);

    /** Returns the percentil value
     *
     * This function returns the value within the range for which a given
     * percent mass of the total number of elements are below it. For instance,
     * if we have 120 values distributed between 0 and 45, and we ask for the
     * percentil of 60%, the returned value is that within 0 and 45 for which
     * 120 * 0.6 = 72 elements are smaller than it.
     *
     * @code
     * percentil60=hist.percentil(60);
     * @endcode
     */
    double percentil(double percent_mass);

    /** Mass below
     *
     * Returns the number of points which are below a certain value
     */
    double mass_below(double value);

    /** Mass above
     *
     * Returns the number of points which are above a certain value
     */
    double mass_above(double value)
    {
        return no_samples - mass_below(value);
    }

    /** Show an histogram
     *
     * The first column is the value associated to each histogram measure. The
     * second one is the histogram measure.
     */
    friend std::ostream& operator<<(std::ostream& o, const histogram1D& hist);

    /** Write an histogram to disk.
     */
    void write(const FileName& fn);

    /** Value --> Index
     *
     * Given a value it returns the code of the interval where it should be
     * counted. If it is outside the global range of the histogram it returns
     * -1
     *
     * @code
     * hist.val2index(12.3, interval_for_it);
     * @endcode
     */
    void val2index(double v, int& i) const
    {
        if (v == hmax)
            i = stepNo() - 1;
        else
            i = (int) FLOOR((v - hmin) / step_size);

        if (i < 0 || i >= stepNo())
            i = -1;
    }

    /** Index --> Value
     *
     * Given the code of one interval, this function returns the value of its
     * starting point (its left border point). If the intervals are defined as
     * [a,b), this function returns a.
     *
     * @code
     * hist.index2val(0, beginning_of_interval_0);
     * @endcode
     */
    void index2val(double i, double& v) const
    {
        v = hmin + i * step_size;
    }

    /** Minimum value where the histogram is defined
     *
     * @code
     * std::cout << "Minimum value for histogram " << hist.min() << std::endl;
     * @endcode
     */
    double hist_min() const
    {
        return hmin;
    }

    /** Maximum value where the histogram is defined
     *
     * @code
     * std::cout << "Maximum value for histogram " << hist.max() << std::endl;
     * @endcode
     */
    double hist_max() const
    {
        return hmax;
    }

    /** Step size for the histogram
     *
     * @code
     * std::cout << "Step size of the histogram " << hist.step() << std::endl;
     * @endcode
     */
    double step() const
    {
        return step_size;
    }

    /** Number of steps in the histogram
     *
     * @code
     * std::cout << "No. Steps in the histogram " << hist.stepNo() << std::endl;
     * @endcode
     */
    int stepNo() const
    {
        return XSIZE(*this);
    }

    /** Number of samples introduced in the histogram
     *
     * @code
     * std::cout << "No. Samples in the histogram " << hist.sampleNo() << std::endl;
     * @endcode
     */
    double sampleNo() const
    {
        return no_samples;
    }
};

/// @defgroup HistogramFunctions1D Functions related to histograms 1D
/// @ingroup Histograms

/** Compute histogram of a vector within its minimum and maximum value
 * @ingroup HistogramFunctions1D
 *
 * Given an array as input, this function returns its histogram within the
 * minimum and maximum of the array, in this way all the values in the array are
 * counted. The array can be of any numerical type (short int, int, double, ...)
 * and dimension. The number of steps must always be given.
 *
 * @code
 * histogram1D hist;
 * compute_hist(v, hist, 100);
 * @endcode
 */
template<typename T>
void compute_hist(T& array, histogram1D& hist, int no_steps)
{
    double min, max;
    array.computeDoubleMinMax(min, max);
    compute_hist(array, hist, min, max, no_steps);
}

/** Compute histogram of the array within two values
 * @ingroup HistogramFunctions1D
 *
 * Given a array as input, this function returns its histogram within two
 * values, the array values outside this range are not counted. This can be used
 * to avoid the effect of outliers which causes a "compression" in the
 * histogram. The array can be of any numerical type (short int, int, double,
 * ...). The number of steps must always be given.
 *
 * @code
 * histogram1D hist;
 * compute_hist(v, hist, 0, 1, 100);
 * @endcode
 */
template<typename T>
void compute_hist(const T& v, histogram1D& hist,
                  double min, double max, int no_steps)
{
    hist.init(min, max, no_steps);
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
    hist.insert_value(MULTIDIM_ELEM(v, i));
}

/** Compute histogram within a region
 * @ingroup HistogramFunctions1D
 *
 * The region is defined by its corners
 */
template<typename T>
void compute_hist(const Matrix2D< T >& v, histogram1D& hist,
                  const Matrix1D< int >& corner1,
                  const Matrix1D< int >& corner2,
                  int no_steps = 100)
{
    double min, max;
    v.computeDoubleMinMax(min, max, corner1, corner2);
    hist.init(min, max, no_steps);

    Matrix1D< int > r(2);
    FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1, corner2)
    hist.insert_value(v(r));
}

/** Compute histogram within a region 3D
 * @ingroup HistogramFunctions1D
 *
 * The region is defined by the corners (k0, i0, j0) and (kF, iF, jF).
 */
template<typename T>
void compute_hist(const Matrix3D< T >& v, histogram1D& hist,
                  const Matrix1D< int >& corner1,
                  const Matrix1D< int >& corner2,
                  int no_steps = 100)
{
    double min, max;
    v.computeDoubleMinMax(min, max, corner1, corner2);
    hist.init(min, max, no_steps);

    Matrix1D< int >r(3);
    FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1, corner2)
    hist.insert_value(v(r));
}

/** Compute the detectability error between two pdf's
 * @ingroup HistogramFunctions1D
 *
 * The input histograms are expressed as probability density functions
 * representing two different objects with the same parameter of variation, for
 * instance, the protein and background both defined by a grey-level. The first
 * histogram would be the histogram of the grey-levels associated to voxels
 * which we know to belong to the protein area, while the second histogram is
 * for the grey-level of voxels which we know to belong to the background. Then
 * the intersection area between the both associated pdf's represents the
 * probability of comitting an error when classifying a voxel, ie, the
 * probability of saying that a voxel is protein when it really is background
 * and viceversa. This function allows you to compute this probability error
 * when the two histograms are provided. Be careful that this probability
 * usually has got a very low value.
 *
 * @code
 * detect_error = detectability_error(hist1, hist2);
 * @endcode
 */
double detectability_error(const histogram1D& h1, const histogram1D& h2);

/** Returns the effective range of a multidimensional array
 * @ingroup HistogramFunctions1D
 *
 * The effective range is defined as the difference of those two values
 * comprising a given central percentage of the array histogram. This function
 * is used to compute the range removing outliers. The default central
 * percentage is 99.75%, although this value should be increased as the number
 * of values in the array decreases. For the default, for instance, the 0.125%
 * of the smaller values are left out as well as the 0.125% of the higher
 * values. The range is given always as a double number.
 *
 * @code
 * double range = v.effective_range();
 * // range for the 99.75% of the mass
 *
 * double range = v.effective_range(1);
 * // range for the 99% of the mass
 * @endcode
 */
template<typename T>
double effective_range(const T& v, double percentil_out = 0.25)
{
    histogram1D hist;
    compute_hist(v, hist, 200);
    double min_val = hist.percentil(percentil_out / 2);
    double max_val = hist.percentil(100 - percentil_out / 2);
    return max_val - min_val;
}

/** Clips the array values within the effective range
 * @ingroup HistogramFunctions1D
 *
 * Look at the documentation of effective_rage
 */
template<typename T>
void reject_outliers(T& v, double percentil_out = 0.25)
{
    histogram1D hist;
    compute_hist(v, hist, 200);
    double eff0 = hist.percentil(percentil_out / 2);
    double effF = hist.percentil(100 - percentil_out / 2);

    // FIXME No words for this...
#define vi MULTIDIM_ELEM(v, i)

    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
    if (vi < eff0)
        vi = eff0;
    else if (vi > effF)
        vi = effF;
#undef vi
}

/** Histogram equalization and re-quantization
 * @ingroup HistogramFunctions1D
 *
 * This function equalizes the histogram of the input multidimensional array,
 * and re-quantize the input array to a specified number of bins. The output
 * array is defined between 0 and bins-1.
 */
template<typename T>
void histogram_equalization(T& v, int bins = 8)
{
    const int hist_steps = 200;
    histogram1D hist;
    compute_hist(v, hist, hist_steps);

    // Compute the distribution function of the pdf
    Matrix1D< double > norm_sum(hist_steps);
    norm_sum(0) = hist(0);

    for (int i = 1; i < hist_steps; i++)
        norm_sum(i) = norm_sum(i - 1) + hist(i);

    norm_sum /= MULTIDIM_SIZE(v);

    // array to store the boundary pixels of bins
    Matrix1D< double > div(bins - 1);
    int index = 0;

    for (int current_bin = 1; current_bin < bins; current_bin++)
    {
        double current_value = (double) current_bin / bins;

        while (norm_sum(index) < current_value)
            index++;

        hist.index2val((double) index, div(current_bin - 1));
    }

    // requantize and equalize histogram
    // FIXME Oh, my, again...
#define vi MULTIDIM_ELEM(v, i)

    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
    if (vi < div(0))
        vi = 0;
    else if (vi > div(bins - 2))
        vi = bins - 1;
    else
    {
        index = 0;

        while (vi > div(index))
            index++;

        vi = index;
    }
#undef vi
}

/** Histograms with 2 parameters
 * @ingroup Histograms
 *
 * The histogram with 2 parameters can be regarded as an approximation to the
 * joint probability density function of two variables (just dividing the
 * histogram by its total mass). Ie, the 2D histogram is the count of times that
 * a certain combination of values in two variables has ocurred. For example,
 * this 2D histograms can be used to plot the projection distribution over the
 * topological sphere, in such situation only the first two Euler angles are
 * interesting and we could plot how many projections are there with first angle
 * equal to 45º and second equal to 0º, and so on covering the whole range for
 * each angle.
 *
 * The 2D histogram is defined, as in the 1D case, by the respective ranges for
 * the two variables and the number of intervals on each range. The distribution
 * and limits of intervals are just the 2D extension of the graph shown in
 * histograms 1D.
 *
 */
class histogram2D : public Matrix2D< double >
{
public:
    // Structure
    double imin; // minimum value of the i axis
    double imax; // maximum value of the i axis
    double istep_size; // size of step
    double jmin; // minimum value of the j axis
    double jmax; // maximum value of the j axis
    double jstep_size; // size of step
    int no_samples; // No. of points inside the histogram

    /** Empty constructor
     *
     * Creates an empty histogram. Before using it you must initialise it with
     * init
     *
     * @code
     * histogram2D hist;
     * @endcode
     */
    histogram2D()
    {
        clear();
    }

    /** Copy constructor
     *
     * Makes an exact copy of the given histogram into another histogram.
     *
     * @code
     * histogram2D hist2(hist1);
     * @endcode
     */
    histogram2D(const histogram2D& H)
    {
        *this = H;
    }

    /** Empties an histogram
     *
     * Remember to initialise it before using it again.
     *
     * @code
     * hist.clear();
     * @endcode
     */
    void clear();

    /** Assignment.
     */
    histogram2D& operator=(const histogram2D& H);

    /** Another function for assigment.
     */
    void assign(const histogram2D& H);

    /** Initialisation of the histogram
     *
     * This is the operation which allows the histogram to be used. This should
     * be performed before inserting any value in it. The information given to
     * this initialisation is the range for each variable within which the
     * values will be counted, and the number of steps (discrete bars) in these
     * ranges.If the value is outside the defined ranges it will not be taken
     * into account although we have asked for its insertion in the histogram.
     *
     * @code
     * hist.init(0, 90, 100, 0, 360, 200);
     * // 100 steps in the range V=0...90 and 200 steps for U=0...360
     * @endcode
     */
    void init(double imin_val, double imax_val, int in_steps,
              double jmin_val, double jmax_val, int jn_steps);

    /** Insert a value within histogram
     *
     * The right interval is chosen according to the initialisation of the
     * histogram and the count of elements in that interval is incremented by 1.
     * If the value lies outside the global range of the histogram nothing is
     * done. Notice that two values are needed, these are the two values of the
     * two variables in our example of the projection at with first Euler
     * angle=45 and second=0, the insertion of this projection in the 2D
     * histogram would be like in the following example.
     *
     * @code
     * hist.insert_value(45, 0);
     * @endcode
     */
    void insert_value(double v, double u);

    /** Show an histogram
     *
     * The first column and second column are the (X,Y) coordinates of each
     * histogram measure. The third one is the histogram measure.
     */
    friend std::ostream& operator<<(std::ostream& o, const histogram2D& hist);

    /** Write an histogram to disk
     */
    void write(const FileName& fn);

    /** Value --> Index
     *
     * Given two values for the two variables it returns the code of the
     * interval where it should be counted. If it is outside the global range of
     * the histogram it returns -1 for that variable (i or j). The following
     * example tells you that the interval corresponding to (45,0) is that with
     * code (i,j), i.e, you could access to hist()(i,j) to know how many
     * projections are there in the same interval as this projection.
     *
     * @code
     * hist.val2index(45, 0, i, j);
     * @endcode
     */
    void val2index(double v, double u, int& i, int& j) const
    {
        if (v == imax)
            i = IstepNo() - 1;
        else
            i = (int) FLOOR((v - imin) / istep_size);

        if (u == jmax)
            j = JstepNo() - 1;

        j = (int) FLOOR((u - jmin) / jstep_size);

        if (i < 0 || i >= IstepNo())
            i = -1;

        if (j < 0 || j >= JstepNo())
            j = -1;
    }

    /** Index --> Value
     *
     * Given the code of one interval, this function returns the value of its
     * starting point (its left-top border point, ie, its lowest corner). If the
     * intervals are defined as [v0,vF) and [u0,uF), this function returns the
     * point (v0,u0)
     *
     * @code
     * hist.index2val(5, 1, v, u);
     * @endcode
     */
    void index2val(double i, double j, double& v, double& u) const
    {
        v = imin + i * istep_size;
        u = jmin + j * jstep_size;
    }

    /** Minimum i value where the histogram is defined
     *
     * @code
     * std::cout << "Minimum value for histogram " << hist.Imin() << std::endl;
     * @endcode
     */
    double Ihist_min() const
    {
        return imin;
    }

    /** Maximum i value where the histogram is defined
     *
     * @code
     * std::cout << "Maximum value for histogram " << hist.Imax() << std::endl;
     * @endcode
     */
    double Ihist_max() const
    {
        return imax;
    }

    /** Step size in i for the histogram
     *
     * @code
     * std::cout << "Step size of the histogram " << hist.Istep() << std::endl;
     * @endcode
     */
    double Istep() const
    {
        return istep_size;
    }

    /** Number of steps in i in the histogram
     *
     * @code
     * std::cout << "No. Steps in the histogram " << hist.IstepNo() << std::endl;
     * @endcode
     */
    int IstepNo() const
    {
        return YSIZE(*this);
    }

    /** Minimum j value where the histogram is defined
     *
     * @code
     * std::cout << "Minimum value for histogram " << hist.Jmin() << std::endl;
     * @endcode
     */
    double Jhist_min() const
    {
        return jmin;
    }

    /** Maximum j value where the histogram is defined
     *
     * @code
     * std::cout << "Maximum value for histogram " << hist.Jmax() << std::endl;
     * @endcode
     */
    double Jhist_max() const
    {
        return jmax;
    }

    /** Step size in j for the histogram
     *
     * @code
     * std::cout << "Step size of the histogram " << hist.Jstep() << std::endl;
     * @endcode
     */
    double Jstep() const
    {
        return jstep_size;
    }

    /** Number of steps in j in the histogram
     *
     * @code
     * std::cout << "No. Steps in the histogram " << hist.JstepNo() << std::endl;
     * @endcode
     */
    int JstepNo() const
    {
        return XSIZE(*this);
    }

    /** Number of samples introduced in the histogram
     *
     * @code
     * std::cout << "No. Samples in the histogram " << hist.sampleNo() << std::endl;
     * @endcode
     */
    int sampleNo() const
    {
        return no_samples;
    }
};

/** @defgroup HistogramFunctions2D Functions related to histograms 2D
 * @ingroup Histograms
 *
 * The vectors can be of any numerical type (short int, int, double, ...), but
 * both of the same type. Vectors must be of the same shape, the first element
 * of v1 and the first of v2 define the position were the first point will be
 * inserted in the histogram, then the second of v1 and of v2, ... That is, the
 * elements of v1 and v2 serve as coordinates within the histogram. The number
 * of steps must always be given.
 */

/** Compute histogram of two arrays within their minimum and maximum values
 * @ingroup HistogramFunctions2D
 *
 * Given two arrays as input, this function returns their joint histogram within
 * the minimum and maximum of the arrays, in this way all the values in the
 * arrays are counted. Both arrays must have the same shape
 */
template<typename T>
void compute_hist(const T& v1, const T& v2,
                  histogram2D& hist, int no_steps1, int no_steps2)
{
    double min1, max1;
    v1.computeDoubleMinMax(min1, max1);

    double min2, max2;
    v2.computeDoubleMinMax(min2, max2);

    compute_hist(v1, v2, hist, min1, max1, min2, max2, no_steps1, no_steps2);
}

/** Compute histogram of two arrays within given values
 * @ingroup HistogramFunctions2D
 *
 * Given two arrays as input, this function returns their joint histogram
 * within the specified values, all the values lying outside are not counted
 */
template<typename T>
void compute_hist(const T& v1, const T& v2, histogram2D& hist,
                  double m1, double M1, double m2, double M2, int no_steps1,
                  int no_steps2)
{
    if (!v1.sameShape(v2))
        REPORT_ERROR(1, "compute_hist: v1 and v2 are of different shape");

    hist.init(m1, M1, no_steps1, m2, M2, no_steps2);

    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v1)
    hist.insert_value(MULTIDIM_ELEM(v1, i), MULTIDIM_ELEM(v2, i));
}

#endif
