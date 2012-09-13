/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "histogram.h"

/* ------------------------------------------------------------------------- */
/* HISTOGRAMS 1D                                                             */
/* ------------------------------------------------------------------------- */
/* Clear ------------------------------------------------------------------- */
void histogram1D::clear()
{
    hmin = 0;
    hmax = 0;
    step_size = 0;
    no_samples = 0;
    MultidimArray<double>::clear();
}

/* Assignment -------------------------------------------------------------- */
histogram1D & histogram1D::operator =(const histogram1D &H)
{
    if (this != &H)
    {
        this->MultidimArray<double>::operator = (H);
        hmin       = H.hmin;
        hmax       = H.hmax;
        step_size  = H.step_size;
        no_samples = H.no_samples;
    }
    return *this;
}

/* Another function for assignament ---------------------------------------- */
void histogram1D::assign(const histogram1D &H)
{
    *this = H;
}
/* Initialize -------------------------------------------------------------- */
void histogram1D::init(double min_val, double max_val, int n_steps)
{
    hmin = min_val;
    hmax = max_val;
    step_size = (double)(max_val - min_val) / (double) n_steps; // CO: n_steps-1->n_steps
    MultidimArray<double>::initZeros(n_steps);
    no_samples = 0;
}

/* Insert value ------------------------------------------------------------ */
//#define DEBUG
void histogram1D::insert_value(double val)
{
    int i;
    val2index(val, i);
    if (i == -1)
        return; // the value is outside our scope
    A1D_ELEM(*this, i)++;
    no_samples++;
#ifdef DEBUG

    std::cout << "   hmin " << hmin << " hmax " << hmax << " value " << val
    << " index " << i << " out of " << no_steps << std::endl;
#endif
}
#undef DEBUG

/* std::cout << hist ------------------------------------------------------------ */
std::ostream& operator << (std::ostream &o, const histogram1D &hist)
{
    MultidimArray<double> aux;
    aux.resize(hist.stepNo(), 2);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(hist)
    {
        hist.index2val(i, A2D_ELEM(aux, i, 0));
        A2D_ELEM(aux, i, 1) = A1D_ELEM(hist, i);
    }
    o << aux;
    return o;
}

/* Write to file ----------------------------------------------------------- */
void histogram1D::write(const FileName &fn)
{
    std::ofstream  fh;
    fh.open(fn.c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR(1, (std::string)"Histogram1D::write: File " + fn +
                     " cannot be openned for output");
    fh << *this;
    fh.close();
}

/* Percentil --------------------------------------------------------------- */
/* This function returns the value of the variable under which the mass% of
   the histogram is comprised */
double histogram1D::percentil(double percent_mass)
{
    int i = 0;
    double acc = 0;
    double required_mass;
    double percentil_i;
    double ret_val;

    // Check it is a correct mass
    if (percent_mass > 100)
        REPORT_ERROR(2001, "Asked for a percentil greater than 100");

    // Trivial cases
    if (percent_mass == 0)
        return(hmin);
    if (percent_mass == 100)
        return(hmax);

    // Any other case, find index of corresponding piece
    required_mass = (double)no_samples * percent_mass / 100.0;
    int N_diff_from_0 = 0;
    while (acc < required_mass)
    {
        acc += A1D_ELEM(*this, i);
        if (A1D_ELEM(*this, i) > 0)
            N_diff_from_0++;
        i++;
    }

    // If the sum is just the one we want OK
    if (acc == required_mass)
        percentil_i = i;
    // If there is only one sample different from 0
    // then there is no way of setting the threshold in the middle
    // Let's put it at the beginning of the bin
    else if (N_diff_from_0 == 1)
        percentil_i = i - 1;
    // If not, then go back a step and compute which fraction of the
    // bar is needed to finish the required mass
    else
    {
        /* CO: We cannot assure that there is at least what is supposed to be
               above this threshold. Let's move to the safe side
        i--;
        acc -= A1D_ELEM(*this,i);
        percentil_i=i+(required_mass-acc)/(double) A1D_ELEM(*this,i); */
        percentil_i = i - 1;
    }

    // Now translate from index to range
    index2val(percentil_i, ret_val);
    return ret_val;
}

/* Mass below -------------------------------------------------------------- */
double histogram1D::mass_below(double value)
{
    // Trivial cases
    if (value <= hmin)
        return 0;
    if (value >= hmax)
        return no_samples;

    // Any other case, find index of corresponding piece
    int i = 0;
    double acc = 0;
    double current_value;
    index2val(i, current_value);
    while (current_value <= value)
    {
        acc += A1D_ELEM(*this, i);
        i++;
        index2val(i, current_value);
    }
    return acc;
}

/* Entropy ----------------------------------------------------------------- */
double histogram1D::entropy() const
{
    MultidimArray<double> p;
    p.initZeros(XSIZE(*this));
    double pSum=0;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(p)
    {
        A1D_ELEM(p,i)=A1D_ELEM(*this,i)+1;
        pSum+=A1D_ELEM(p,i);
    }
    double entropy=0;
    double ipSum=1.0/pSum;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(p)
    {
        double pi=A1D_ELEM(p,i)*ipSum;
        entropy-=pi*log(pi);
    }
    return entropy;
}

/* Detectability error ----------------------------------------------------- */
// Given two histograms (probability density functions) this function returns
// the detection error for using both, ie, the probability that given a sample
// I would assign it to class 1 when it belongs to class 2 and viceversa.
// This is computed by the calculation of the overlapping areas.
double detectability_error(const histogram1D &h1, const histogram1D &h2)
{
    double hmin, hmax;
    double step;
    double v;
    double error = 0;
    int   ih1, ih2;               // Indexes within the histograms
    double p1, p2;                 // Probability associated

    // Find global range
    hmin = XMIPP_MAX(h1.hmin, h2.hmin);
    hmax = XMIPP_MIN(h1.hmax, h2.hmax);
    step = XMIPP_MIN(h1.step_size, h2.step_size) / 2;

    // Go over the range computing the errors
    v = hmin;
    int N = 0;
    while (v <= hmax)
    {
        h1.val2index(v, ih1);
        p1 = A1D_ELEM(h1, ih1) / h1.no_samples;
        h2.val2index(v, ih2);
        p2 = A1D_ELEM(h2, ih2) / h2.no_samples;
        //#define DEBUG
#ifdef DEBUG

        std::cout << "Comparing at " << v << " (" << ih1 << ") p1=" << p1 << " p2= " << p2 << std::endl;
        std::cout << "   hmin " << hmin << " hmax " << hmax << " stepsize " << h1.step_size << std::endl;
#endif//;

        if (p1 != 0 && p2 != 0)
            if (p1 > p2)
                error += p2;
            else
                error += p1;
        v += step;
        N++;
    }

    // Normalise such that the result is the area of a probability function
    error *= step / (hmax - hmin);
    error /= N;
#ifdef DEBUG

    std::cout << "Total error = " << error << std::endl;
#endif

    return error;
}

/* Kullback Leibler distance ----------------------------------------------- */
double KLDistance(const histogram1D& h1, const histogram1D& h2)
{
    if (XSIZE(h1)!=XSIZE(h2))
        REPORT_ERROR(1,"KLDistance: Histograms of different sizes");
    
    double retval=0;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(h1)
        if (h2(i)!=0.0 && h1(i)!=0.0) retval += h1(i)*log10(h1(i)/h2(i)); 
    return retval;
}

/* ------------------------------------------------------------------------- */
/* IRREGULAR HISTOGRAMS                                                      */
/* ------------------------------------------------------------------------- */
/* Initialization ---------------------------------------------------------- */
void IrregularHistogram1D::init(const histogram1D &hist,
    const MultidimArray<int> &bins)
{
    int steps_no = XSIZE(bins);
    __binsRightLimits.initZeros(steps_no);
    __hist.initZeros(steps_no);

    int k = 0;
    for (int i = 0; i < steps_no; i++)
    {
        hist.index2val(bins(i), __binsRightLimits(i));
        __hist(i) = 0;
        for (int j = k; j <= bins(i); j++)
            __hist(i) += hist(j);
        k = bins(i) + 1;
    }
}

/* val2index --------------------------------------------------------------- */
int IrregularHistogram1D::val2Index(double value)
{
    int binsNo = XSIZE(__binsRightLimits);
    for (int i = 0; i < binsNo; i++)
        if(value <= __binsRightLimits(i)) return i;
    
    //In case the value is greater, we return the last bin
    return binsNo - 1;
}

/* Normalization ----------------------------------------------------------- */
void IrregularHistogram1D::selfNormalize()
{
    __hist /= __hist.sum();
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out,
    const IrregularHistogram1D &_hist)
{
    for (int i = 0; i < XSIZE(_hist.__binsRightLimits); i++)
        _out << "\t" << _hist.__binsRightLimits(i) << "\t\t" << _hist.__hist(i) << std::endl;
    return _out;
}

/* Get value --------------------------------------------------------------- */
double IrregularHistogram1D::operator()(int i) const
{
    return __hist(i);
}

/* Get value --------------------------------------------------------------- */
const histogram1D& IrregularHistogram1D::getHistogram() const
{
    return __hist;
}

/* ------------------------------------------------------------------------- */
/* HISTOGRAMS 2D                                                             */
/* ------------------------------------------------------------------------- */
/* Clear ------------------------------------------------------------------- */
void histogram2D::clear()
{
    imin = 0;
    imax = 0;
    istep_size = 0;
    jmin = 0;
    jmax = 0;
    jstep_size = 0;
    no_samples = 0;
    MultidimArray<double>::clear();
}

/* Assignment -------------------------------------------------------------- */
histogram2D & histogram2D::operator = (const histogram2D &H)
{
    if (this != &H)
    {
        this->MultidimArray<double>::operator =(H);
        imin        = H.imin;
        imax        = H.imax;
        istep_size  = H.istep_size;
        jmin        = H.jmin;
        jmax        = H.jmax;
        jstep_size  = H.jstep_size;
        no_samples  = H.no_samples;
    }
    return *this;
}

/* Another function for assignment -------------------------------------------------------------- */
void histogram2D::assign(const histogram2D &H)
{
    *this = H;
}

/* Initialize -------------------------------------------------------------- */
void histogram2D::init(double imin_val, double imax_val, int in_steps,
                       double jmin_val, double jmax_val, int jn_steps)
{
    // V axis
    imin = imin_val;
    imax = imax_val;
    istep_size = (double)(imax_val - imin_val) / (double) in_steps;

    // U axis
    jmin = jmin_val;
    jmax = jmax_val;
    jstep_size = (double)(jmax_val - jmin_val) / (double) jn_steps;

    initZeros(in_steps, jn_steps);
    no_samples = 0;
}

/* Insert value ------------------------------------------------------------ */
void histogram2D::insert_value(double v, double u)
{
    int i, j;
    val2index(v, u, i, j);
    if (i == -1 || j == -1)
        return; // it is outside our scope
    i = CLIP(i, 0, YSIZE(*this));
    j = CLIP(j, 0, XSIZE(*this));
    A2D_ELEM(*this, i, j)++;
    no_samples++;
}

/* std::cout << hist ------------------------------------------------------------ */
std::ostream& operator << (std::ostream &o, const histogram2D &hist)
{
    MultidimArray<double> aux;
    aux.resize(hist.IstepNo()*hist.JstepNo(), 3);
    int n = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(hist)
    {
        hist.index2val(i, j, A2D_ELEM(aux, n, 0), A2D_ELEM(aux, n, 1));
        A2D_ELEM(aux, n, 2) = A2D_ELEM(hist, i, j);
        n++;
    }
    o << aux;
    return o;
}

/* Write to file ----------------------------------------------------------- */
void histogram2D::write(const FileName &fn)
{
    std::ofstream  fh;
    fh.open(fn.c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR(1, (std::string)"histogram2D::write: File " + fn + " cannot be openned for output");
    fh << *this;
    fh.close();
}
