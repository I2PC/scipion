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

#ifndef XMIPP_FUNCS_H
#define XMIPP_FUNCS_H

#include <fstream>
#include <typeinfo>
#include <complex>
#include <algorithm>
#include "xmipp_filename.h"
#include "xmipp_macros.h"
#include "xmipp_error.h"



// For timing functions
// Uncomment next line timing functions are giving problems in your system
//#define _NO_TIME
#ifndef _NO_TIME
#include <unistd.h>
#ifndef __MINGW32__
#include <sys/times.h>
#endif
#ifdef _IRIX65
#include <sys/types.h>
#include <time.h>
#endif
#endif



/// @defgroup GeneralFunctions General functions
/// @ingroup DataLibrary
//@{
/// @name Numerical functions
//@{
/** Tabulated Sinc = SIN(PI*X)/(PI*X)
 *
 * A lookup-table with the given sampling rate and range is created.
 *
 * @code
 * tabsinc TSINC(0.0001,64);
 * if (TSINC(1)==0)
 *    std::cout << "This is true!\n";
 * @endcode
 *
 * This class is not ported to Python.
 */
class Tabsinc
{
public:

    double sampl;
    double isampl;
    int xmax;
    int no_elem;
    double* tabulatedsinc;

public:
    /** Constructor with sampling rate and range */
    Tabsinc(const double dd, const int xx)
    {
        sampl = dd;
        isampl = 1.0/sampl;
        xmax = xx;
        filltable();
    }

    // Destructor
    virtual ~Tabsinc()
    {
        delete tabulatedsinc;
    }

#define TSINCVALUE(Tsinc, x,y) \
     { \
      int TSINCVALUEaux=(int)(x * Tsinc.isampl); \
         y=Tsinc.tabulatedsinc[ABS(TSINCVALUEaux)]; \
     }


    /** Value access. Tabulated sine in radians */
    double operator()(double val) const
    {
        int aux=(int)(val * isampl);
        return tabulatedsinc[ABS(aux)];
    }

    /** Actually fill the table */
    void filltable()
    {
        no_elem = (int)(xmax / sampl);
        tabulatedsinc = new double[no_elem];
        tabulatedsinc[0] = 1;
        for (int i = 1; i < no_elem; i++)
        {
            double xx = (double) i * sampl * PI;
            tabulatedsinc[i] = sin(xx) / xx;
        }
    }
};

/** Kaiser-Bessel function
 *
 *  This code was modified from SPARX and originally written by Pawel
 *  Penczek at the University of Texas - Houston Medical School
 *
 *  see P. A. Penczek, R. Renka, and H. Schomberg,
 *      J. Opt. Soc. Am. _21_, 449 (2004)
 *
 *  The I0 version can be tabulated and interpolated upon
 *  demand, but the max error needs to be checked.  The
 *  "vtable" parameter corresponds to the maximum value of x
 *  for which the I0 selfWindow is non-zero.  Setting "vtable"
 *  different from "v" corresponds to a change in units of x.
 *  In practice, it is often handy to replace x in some sort
 *  of absolute units with x described in terms of grid
 *  intervals.
 *
 *  The get_kbsinh_win and get_kbi0_win functions return
 *  single-argument function objects, which is what a
 *  generic routine is likely to want.
 *
 * @code
 *  kb = KaiserBessel(alpha, K, r, v , N);
 *  double wx = kb.sinhwin(32);
 *  double tablex1 = kb.i0win_tab(delx-inxold+3);
 * @endcode
 */
class KaiserBessel
{
protected:
    double alpha, v, r; /** Kaiser-Bessel parameters */
    int N; /** size in Ix-space */
    int K; /** I0 selfWindow size */
    double vtable; /** table I0 non-zero domain maximum */
    int ntable;
    std::vector<double> i0table;
    double dtable; /** table spastd::cing */
    double alphar; /** alpha*r */
    double fac; /** 2*pi*alpha*r*v */
    double vadjust;
    double facadj; /** 2*pi*alpha*r*vadjust */
    void build_I0table(); /** Tabulate I0 selfWindow for speed */
    double fltb;

public:
    /** Constructor with parameters */
    KaiserBessel(double alpha_, int K, double r_,
                 double v_, int N_, double vtable_=0.,
                 int ntable_ = 5999);

    /** Compute the maximum error in the table */
    double I0table_maxerror();
    std::vector<double> dump_table() const
    {
        return i0table;
    }

    /** Kaiser-Bessel Sinh selfWindow function */
    double sinhwin(double x) const;

    /** Kaiser-Bessel I0 selfWindow function */
    double i0win(double x) const;

    /** Kaiser-Bessel I0 selfWindow function (uses table lookup) */
    inline double i0win_tab(double x) const
    {
        double xt;
        if(x<0.)
            xt = -x*fltb+0.5;
        else
            xt = x*fltb+0.5;
        return i0table[ (int) xt];
    }

    /** Return the size of the I0 selfWindow */
    int get_window_size() const
    {
        return K;
    }
};

#if defined(__APPLE__) || defined(__MINGW32__)
/** Calculate sin and cos at the same time
 *
 */
void sincos(double angle, double * sine, double * cosine);
#endif

/** Solve second degree equation
 *
 *  ax^2+bx+c=0
 *
 * It returns the number of real solutions, 0 if the two roots are complex
 * and -1 if the equation is impossible to satisfy (Ex: 7=0). A number is said
 * to be 0 if its absolute magnitude is smaller than precision. This is used to
 * avoid dividing by 0
 */
int solve_2nd_degree_eq(double a,
                        double b,
                        double c,
                        double& x1,
                        double& x2,
                        double prec = XMIPP_EQUAL_ACCURACY);


/**
 *  Structure to define random generation mode
 */
enum RandomMode
{
    RND_UNIFORM = 0,
    RND_GAUSSIAN = 1
} ;

/** 1D gaussian value
 *
 * This function returns the value of a univariate gaussian function at the
 * point x.
 */
double gaussian1D(double x, double sigma, double mu = 0);

/** 1D t-student value
 *
 * This function returns the value of a univariate t-student function at the
 * point x, and with df degrees of freedom
 */
double tstudent1D(double x, double df, double sigma, double mu = 0);

/** Inverse Cumulative distribution function for a Gaussian
 *
 * This function returns the z of a N(0,1) such that the probability below z is p
 *
 * The function employs an fast approximation to z which is valid up to 1e-4.
 * See http://www.johndcook.com/normal_cdf_inverse.html
 */
double icdf_gauss(double p);

/** Cumulative distribution function for a Gaussian
 *
 * This function returns the value of the CDF of a univariate gaussian function at the
 * point x.
 */
double cdf_gauss(double x);

/** Cumulative distribution function for a t-distribution
 *
 * This function returns the value of the CDF of a univariate t-distribution
 * with k degrees of freedom  at the point t.
 *  Adapted by Sjors from: http://www.alglib.net/specialfunctions/distributions/student.php
 */
double cdf_tstudent(int k, double t);

/** Cumulative distribution function for a Snedecor's F-distribution.
 *
 * This function returns the value of the CDF of a univariate Snedecor's
 * F-distribution
 * with d1, d2 degrees of freedom  at the point x.
 */
double cdf_FSnedecor(int d1, int d2, double x);

/** Inverse Cumulative distribution function for a Snedecor's F-distribution.
 *
 * This function returns the value of the ICDF of a univariate Snedecor's
 * F-distribution
 * with d1, d2 degrees of freedom with probability p, i.e., it returns
 * x such that CDF(d1,d2,x)=p
 */
double icdf_FSnedecor(int d1, int d2, double p);

/** 2D gaussian value
 *
 * This function returns the value of a multivariate (2D) gaussian function at
 * the point (x,y) when the X axis of the gaussian is rotated ang
 * (counter-clockwise) radians (the angle is positive when measured from the
 * universal X to the gaussian X). X and Y are supposed to be independent.
 */
double gaussian2D(double x,
                  double y,
                  double sigmaX,
                  double sigmaY,
                  double ang,
                  double muX = 0,
                  double muY = 0);
//@}

/// @name Miscellaneous functions
//@{
/** Divides a number into most equally groups
 *
 * For example you want to distribute N jobs between M workers
 * so each worker will have N/M jobs and some of them(N % M first)
 * will have N/M + 1 jobs
 * So for the worker 'rank' will be computed the first and last job to do
 * Return the number of jobs assigned, that could be N/M + 1 or N/M
 *
 */
size_t divide_equally(size_t N, size_t size, size_t rank, size_t &first, size_t &last);

/** In which group (of divide_equally) is myself situated?
 */
size_t divide_equally_group(size_t N, size_t size, size_t myself);

/** Compute statistics of a std::vector
 */
template <class T>
void computeStats(const std::vector<T> &V, double& avg, double& stddev,
                  T& minval, T& maxval)
{
    if (V.size()<= 0)
        return;

    avg = 0;
    stddev = 0;

    minval = maxval = V[0];

    size_t nmax=V.size();
    const T* ptr=&V[0];
    for(size_t n=0; n<nmax; ++n, ++ptr)
    {
        double val=*ptr;
        avg += val;
        stddev += val * val;

        if (val > maxval)
            maxval = val;
        else if (val < minval)
            minval = val;
    }
    avg /= nmax;

    if (nmax > 1)
    {
        stddev = stddev / nmax - avg * avg;
        stddev *= nmax / (nmax - 1);

        // Foreseeing numerical instabilities
        stddev = sqrt(static_cast< double >(ABS(stddev)));
    }
    else
        stddev = 0;
}

/** Compute statistics of a std::vector
 */
template <class T>
void computeAvgStddev(const std::vector<T> &V, double& avg, double& stddev)
{
    if (V.size()<= 0)
        return;

    avg = 0;
    stddev = 0;

    size_t nmax=V.size();
    const T* ptr=&V[0];
    for(size_t n=0; n<nmax; ++n, ++ptr)
    {
        double val=*ptr;
        avg += val;
        stddev += val * val;
    }
    avg /= nmax;

    if (nmax > 1)
    {
        stddev = stddev / nmax - avg * avg;
        stddev *= nmax / (nmax - 1);

        // Foreseeing numerical instabilities
        stddev = sqrt(static_cast< double >(ABS(stddev)));
    }
    else
        stddev = 0;
}

/** Initialize std::vector with constant value */
template <class T>
void initConstant(std::vector<T> &V, T &value)
{
    const T* ptr=&V[0];
    size_t nmax=V.size();
    for(size_t n=0; n<nmax; ++n, ++ptr)
    	*ptr=value;
}
//@}

/// @name Complex functions
//@{
/** Real/Imaginary --> Complex
 *
 * The output array(s) must be already resized.
 *
 * This function is not ported to Python.
 */
template <typename T>
void RealImag2Complex(const T* _real,
                      const T* _imag,
                      std::complex< double >* _complex,
                      int length)
{
    T* aux_real = (T*) _real;
    T* aux_imag = (T*) _imag;
    double* aux_complex = (double*) _complex;

    for (int i = 0; i < length; i++)
    {
        *aux_complex++ = (double)(*aux_real++);
        *aux_complex++ = (double)(*aux_imag++);
    }
}

/** Amplitude/Phase --> Complex
 *
 * The output array(s) must be already resized.
 *
 * This function is not ported to Python.
 */
template <typename T>
void AmplPhase2Complex(const T* _ampl,
                       const T* _phase,
                       std::complex< double >* _complex,
                       int length)
{
    T* aux_ampl = (T*) _ampl;
    T* aux_phase = (T*) _phase;
    double* aux_complex = (double*) _complex;

    for (int i = 0; i < length; i++)
    {
        double ampl = (double)(*aux_ampl++);
        double phase = (double)(*aux_phase++);
        *aux_complex++ = ampl * cos(phase);
        *aux_complex++ = ampl * sin(phase);
    }
}

/** Complex --> Real/Imag
 *
 * The output array(s) must be already resized.
 *
 * This function is not ported to Python.
 */
template <typename T>
void Complex2RealImag(const std::complex< double >* _complex,
                      T* _real,
                      T* _imag,
                      int length)
{
    T* aux_real = (T*) _real;
    T* aux_imag = (T*) _imag;
    double* aux_complex = (double*) _complex;

    for (int i = 0; i < length; i++)
    {
        *aux_real++ = (T)(*aux_complex++);
        *aux_imag++ = (T)(*aux_complex++);
    }
}

/** Complex --> Amplitude/Phase
 *
 * The output array(s) must be already resized.
 *
 * This function is not ported to Python.
 */
template <typename T>
void Complex2AmplPhase(const std::complex< double >* _complex,
                       T* _ampl,
                       T* _phase,
                       int length)
{
    T* aux_ampl = (T*) _ampl;
    T* aux_phase = (T*) _phase;
    double* aux_complex = (double*) _complex;

    for (int i = 0; i < length; i++)
    {
        double re = *aux_complex++;
        double im = *aux_complex++;
        *aux_ampl++ = sqrt(re * re + im * im);
        *aux_phase++ = atan2(im, re);
    }
}
//@}

/** @name Random functions
 *
 * These functions allow you to work in an easier way with the random functions
 * of the Numerical Recipes. Only an uniform and a gaussian random number
 * generators have been implemented. In fact only a uniform generator exists and
 * the gaussian one is based on a call to it. For this reason, if you initialize
 * the gaussian random generator, you are also initialising the uniform one.
 *
 * Here goes an example for uniform random numbers to show how to use this set
 * of functions.
 *
 * @code
 * // Initialise according to the clock
 * randomize_random_generator();
 *
 * // Show 10 random numbers between -1 and 1
 * for (int i=0; i<10; i++)
 *     std::cout << rnd_unif(-1,1) << std::endl;
 * @endcode
 */
//@{
/** Reset uniform random generator to a known point
 *
 * If you initialize the random generator with this function each time, then the
 * same random sequence will be generated
 *
 * @code
 * init_rnd_unif();
 * init_rnd_unif(17891)
 * @endcode
 */
void init_random_generator(int seed = -1);

/** Reset random generator according to the clock.
 *
 * This time the initialisation itself assures a random sequence different each
 * time the program is run. Be careful not to run the program twice within the
 * same second as the initialisation will be the same for both runs.
 * Returns seed.
 */
unsigned int randomize_random_generator();

/** Produce a uniform random number between 0 and 1
 *
 * @code
 * std::cout << "This random number should be between 0 and 1: " << rnd_unif()
 * << std::endl;
 * @endcode
 */
double rnd_unif();

/** Produce a uniform random number between a and b
 *
 * @code
 * std::cout << "This random number should be between 0 and 10: " << rnd_unif(0,10)
 * << std::endl;
 * @endcode
 */
double rnd_unif(double a, double b);

/** Produce a t-distributed random number with mean 0 and standard deviation 1 and nu degrees of freedom
 *
 * @code
 * std::cout << "This random number should follow t(0,1) with 3 degrees of freedon: " << rnd_student_t(3.)
 * << std::endl;
 * @endcode
 */
double rnd_student_t(double nu);

/** Produce a gaussian random number with mean a and standard deviation b and nu degrees of freedom
 *
 * @code
 * std::cout << "This random number should follow t(1,4) with 3 d.o.f.: " << rnd_gaus(3,1,2)
 * << std::endl;
 * @endcode
 */
double rnd_student_t(double nu, double a, double b);

/** Produce a gaussian random number with mean 0 and standard deviation 1
 *
 * @code
 * std::cout << "This random number should follow N(0,1): " << rnd_gaus()
 * << std::endl;
 * @endcode
 */
double rnd_gaus();

/** Produce a gaussian random number with mean a and standard deviation b
 *
 * @code
 * std::cout << "This random number should follow N(1,4): " << rnd_gaus(1,2)
 * << std::endl;
 * @endcode
 */
double rnd_gaus(double a, double b);

/** Gaussian area from -x0 to x0
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * x0 must be positive
 */
double gaus_within_x0(double x0, double mean = 0, double stddev = 1);

/** Gaussian area outisde -x0 to x0
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * x0 must be positive
 */
double gaus_outside_x0(double x0, double mean = 0, double stddev = 1);

/** Gaussian area from -inf to x0
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * There is no restriction over the sign of x0
 */
double gaus_up_to_x0(double x0, double mean = 0, double stddev = 1);

/** Gaussian area from x0 to inf
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * There is no restriction over the sign of x0
 */
double gaus_from_x0(double x0, double mean = 0, double stddev = 1);

/** t0 for a given two-sided probability
 *
 * This function returns t0 such that the student probability outside t0 is
 * equal to p
 */
double student_outside_probb(double p, double degrees_of_freedom);

/** student area from -t0 to t0
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * t0 must be positive
 */
double student_within_t0(double t0, double degrees_of_freedom);

/** student area outisde -t0 to t0
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * t0 must be positive
 */
double student_outside_t0(double t0, double degrees_of_freedom);

/** student area from -inf to t0
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * There is no restriction over the sign of t0
 */
double student_up_to_t0(double t0, double degrees_of_freedom);

/** student area from t0 to inf
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * There is no restriction over the sign of t0
 */
double student_from_t0(double t0, double degrees_of_freedom);

/** chi2 area from -inf to t0
 *
 * By default the chi2 mean is 0 and the chi2 standard deviation is 1.
 * There is no restriction over the sign of t0
 */
double chi2_up_to_t0(double t0, double degrees_of_freedom);

/** chi2 area from t0 to inf
 *
 * By default the chi2 mean is 0 and the chi2 standard deviation is 1.
 * There is no restriction over the sign of t0
 */
double chi2_from_t0(double t0, double degrees_of_freedom);

/** Produce a log uniform random number between a and b
 *
 * Watch out that the following inequation must hold 0<a<=b.
 *
 * @code
 * std::cout << "This random number should be between 1 and 1000: "
 * << rnd_log(10,1000)<< std::endl;
 * @endcode
 */
double rnd_log(double a, double b);
//@}

/** @name Filename handling
 *
 * Filenames are something very common in any package, so it might be worthy to
 * devote some effort to make easier their manipulation. In filename conventions
 * you have a detailed explanation of the Filenames dealed. Here you are some
 * examples of what you can do with this class
 *
 * @code
 * FileName fn_vol;
 * FileName fn_root;
 * VolumeXmipp vol_voxels;
 *
 * // Read Volume name
 * fn_vol = argv[1]; // Suppose you give a more complex parameter reading.
 *                   // What we intend to have with this line, is something
 *                   // like fn_vol="art00001.vol"
 * fn_root = fn_vol.without_extension();
 *
 * // Read volume
 * vol_voxels.read(fn_vol);
 *
 * // Process volume
 * ...
 *
 * // Write volume at iteration 1, with a name in the fashion
 * // "art_it1_00001.voxels"
 * vol_voxels.write(fn_root.get_root() + "_it" + integerToString(1) + "_" +
 * fn_root.get_number() + ".voxels");
 * @endcode
 */

/** @name Time managing
 *
 * These functions are used to make time measures of the algorithms. If you know
 * the total amount of work to do then some estimation can be done about how
 * much time is left before finishing. The time functions are very machine
 * dependent, we've tried to accomodate the compilation for several machines,
 * but if still programs do not work, you may configure Xmipp to avoid these
 * time measurements, functions are then substituted by null functions doing
 * nothing.
 *
 * @code
 * // Variable declaration
 * TimeStamp t0;
 *
 * // Beginning of the program
 * time_config();
 * ...
 *
 * annotate_processor_time(&t0);
 * // Part to be measured
 * ...
 *
 * // End of part to be measured
 * print_elapsed_time(t0);
 * @endcode
 *
 * While for an estimation of time to go you can make it in two ways:  one
 * analytical, and another graphical.
 *
 * Analytical:
 *
 * @code
 * // Variable declaration
 * ProcessorTimeStamp t0;
 * double to_go;
 *
 * // Beginning of the program
 * time_config();
 * ...
 *
 * annotate_processor_time(&t0);
 * // Part to be measured
 * for (int i=0; i<60; i++)
 * {
 *     ...
 *     // Compute the time to go with the fraction of work already done
 *     to_go = time_to_go(t0, (double) (i + 1) / 60);
 *     std::cout << "I think you will be here " << to_go << "seconds more\n";
 * }
 * // End of part to be measured
 * print_elapsed_time(t0);
 * @endcode
 *
 * Alternatively:
 *
 * @code
 * // Variable declaration
 * TimeStamp t0;
 * double to_go;
 *
 * annotate_time(&t0);
 * // Part to be measured
 * for (int i=0; i<60; i++)
 * {
 *     ...
 *     // Compute the time to go with the fraction of work already done
 *     to_go = time_to_go(t0, (double) (i + 1) / 60);
 *     std::cout << "I think you will be here " << to_go << "seconds more\n";
 * }
 * // End of part to be measured
 * print_elapsed_time(t0);
 * @endcode
 *
 * Graphical:
 * @code
 * // Beginning of the program
 * time_config();
 * ...
 *
 * // Init the progress bar with the total amount of work to do
 * // It is very important that there is no print out to stdout but
 * // the progress bar
 * init_progress_bar(60);
 *
 * // Part to be measured
 * for (int i=0; i<60; i++)
 * {
 *     ...
 *     progress_bar(i+1);
 * }
 *
 * // In this case the following call is useless since it has been
 * // already done in the loop, but there are cases where a final call
 * // with the total amount of work is not performed and although the
 * // whole task has been finished it seems that it hasn't as the
 * // progress bar hasn't been called with the final work but with
 * // a quantity a little smaller.
 * progress_bar(60);
 * @endcode
 *
 * These functions is not ported to Python.
 */
//@{
#if defined _NO_TIME || defined __MINGW32__
typedef int ProcessorTimeStamp; // Any other kind of data will do
typedef int TimeStamp; // Any other kind of data will do
struct tm* localtime_r (const time_t *clock, struct tm *result);
#else
typedef struct tms ProcessorTimeStamp; // Renaming of the time structure
typedef size_t TimeStamp;              // Timestamp in miliseconds
#endif

/** Read the system clock frequency
 *
 * This operation is needed only once in a program always we want to have a time
 * measure, or an estimation of remaining time.
 *
 * @code
 * time_config();
 * @endcode
 *
 * This function is not ported to Python.
 */
void time_config();

#if !defined _NO_TIME && !defined __MINGW32__
/** Annotate actual time
 *
 * This annotation is used later to compute the elapsed time.
 *
 * @code
 * ProcessorTimeStamp t0;
 * annotate_processor_time(&t0);
 * @endcode
 *
 * This function is not ported to Python.
 */
void annotate_processor_time(ProcessorTimeStamp* time);
#endif

/** Annotate actual time
 *
 * This annotation is used later to compute the elapsed time.
 *
 * @code
 * TimeStamp t0;
 * annotate_time(&t0);
 * @endcode
 *
 * This function is not ported to Python.
 */
void annotate_time(TimeStamp* time);

#if !defined _NO_TIME && !defined __MINGW32__
/** Accumulate time
 *
 * Initially dest_time should be set to orig time. Then you acumulate succesive
 * times calling this function (Destination time=destination_time + (now -
 * original time)) and finally the elapsed time is the dest time minus the first
 * one (the one which initiliazed the dest time.
 *
 * This function is not ported to Python.
 */
void acum_time(ProcessorTimeStamp* orig, ProcessorTimeStamp* dest);
#endif

/** Compute elapsed time since a given annotation
 *
 * Given an annotation of time, this function computes the time elapsed since
 * then in seconds. The annotation is not modified. Usually the time is shown in
 * seconds, but you might specify to show it in clock ticks setting the variable
 * _IN_SECS to FALSE.
 *
 * @code
 * TimeStamp t0;
 * annotate_time(&t0);
 * ...;
 * double elapsed = elapsed_time(t0);
 *
 * TimeStamp t0;
 * annotate_time(&t0);
 * ...;
 * double elapsed = elapsed_time(t0, FALSE);
 * @endcode
 *
 * This function is not ported to Python.
 */
double elapsed_time(ProcessorTimeStamp& time, bool _IN_SECS = true);

#if !defined _NO_TIME && !defined __MINGW32__
/** Show on screen the elapsed time since a given annotation
 *
 * The format of the printing is "Elapsed time: User(13) System(1)" that means
 * that the user has used 13 seconds and the system 1, a total of 14 seconds
 * since the last annotation in this TimeStamp variable.
 *
 * @code
 * ProcessorTimeStamp t0;
 * annotate_processor_time(&t0);
 * ...;
 * print_elapsed_time(t0);
 * @endcode
 *
 * Usually the time is shown in seconds, but you might specify to show it in
 * clock ticks setting the variable _IN_SECS to FALSE.
 *
 * This function is not ported to Python.
 */
void print_elapsed_time(ProcessorTimeStamp& time, bool _IN_SECS = true);
#endif

/** Show on screen the elapsed time since a given annotation
 *
 * The format of the printing is "Elapsed time: User(13) System(1)" that means
 * that the user has used 13 seconds and the system 1, a total of 14 seconds
 * since the last annotation in this TimeStamp variable.
 *
 * @code
 * TimeStamp t0;
 * annotate_time(&t0);
 * ...;
 * print_elapsed_time(t0);
 * @endcode
 *
 * Usually the time is shown in seconds, but you might specify to show it in
 * clock ticks setting the variable _IN_SECS to FALSE.
 *
 * This function is not ported to Python.
 */
void print_elapsed_time(TimeStamp& time, bool _IN_SECS = true);

/** This class will encapsulate the logic to time printing.
 * Useful for debugging.
 */
class Timer
{
private:
  struct timeval  tv;
  struct tm tm;
  size_t tic_time;

public:
  size_t now();
  void tic();
  void toc(const char * msg=NULL);
};

#if !defined _NO_TIME && !defined __MINGW32__
/** Returns the estimated time left to finish
 *
 * To make this estimation the starting time must have been annotated before and
 * the fraction of the total amount of work must be estimated by the programmer.
 * See Time managing for an example.
 *
 * This function is not ported to Python.
 */
double time_to_go(ProcessorTimeStamp& time, double fraction_done);
#endif

/** Initialise the progress bar
 *
 * The progress bar is initialised to count for a total amount of work. For
 * instance, if we are to do something 60 times, the progress bar should be
 * initialised to that value. At the same time the bar is printed with the
 * initial guess of time left (ie, nothing "0000/????"). The number before the
 * slash is the elapsed time since initialisation of the progress bar, while the
 * second number is the estimation of total time that this task will take. See
 * Time managing for a more detailed example.
 *
 * @code
 * init_progress_bar(60);
 * @endcode
 */
void init_progress_bar(long total);

/** Update progress bar
 *
 * With this function you can change the already done amount of work, if
 * something is to be done 60 times and now we have already done 13 then we
 * could tell this to the progress bar with
 *
 * @code
 * progress_bar(13);
 * @endcode
 *
 * The information that this thing was to be done 60 times was given at the
 * initialisation of the progress bar. It is very important that during the use
 * of the progress bar, nobody prints anything to stdout as it is being used by
 * the progress bar. At the end you could make a call to progress_bar with the
 * total amount of work just to make sure that the printout is pretty enough.
 */
void progress_bar(long act_time);

/** Function to get current time as string
 *
 */
char * getCurrentTimeString();


/** BaseListener
 *
 * This class implements the xmipp listener class for notification of progress
 * status and other operations. It is an abstract class that contains base
 * functions useful for implementation of customer-design notification classes.
 *
 * This class is not ported to Python.
 */
class BaseListener
{
public :
    /** Default constructor */
    BaseListener(): verbosity(0), cancel(false)
    {}
    
    /** Destructor */
    virtual ~BaseListener() {}

    /** Initialize progress bar.
     *
     * _est_it: Defines the estimated number of iterations
     *
     * This method will initialize the progress bar with the number of estimated
     * iterations.
     */
    virtual void OnInitOperation(unsigned long _est_it) = 0;

    /** Shows a message indicating the operation in progress
     *
     * Class inheriting from this abstract class can output this message in the
     * way they want.
     *
     * _rsop: string message
     */
    virtual void OnReportOperation(const std::string& _rsOp) = 0;

    /** Show a bar with the progress in time
     *
     * When the input is negative then we are setting the progress bar, this
     * will be the total of elements to process. Afterwards the call to this
     * routine must be in ascending order, ie, 0, 1, 2, ... No. elements
     *
     * Class inheriting from this base class can define their own progress bar
     * here.
     *
     * _it: iteration number
     */
    virtual void OnProgress(unsigned long _it) = 0;

    /** This method will get the verbosity level
     */
    virtual const unsigned& getVerbosity() const
    {
        return verbosity;
    }

    /** This method will set the verbosity level to be used
     */
    virtual unsigned& setVerbosity()
    {
        return verbosity;
    }

    /** This method returns true if a cancel command was set
     *
     * It can be used to check whether an external event is trying to cancel any
     * operation. Inside an algorithm you can can this method to check if a
     * cancel operation was requested.
     */
    virtual const bool& OnUserCancel() const
    {
        return cancel;
    }

    /** This method is used to send a cancel command
     */
    virtual bool& setCancel()
    {
        return cancel;
    }

private:
    unsigned verbosity;
    bool cancel;
};

/** TextualListener
 *
 * This class implements the xmipp classical textual listener class for
 * notification of progress status. It inherits from BaseListener.
 *
 * This class is not ported to Python.
 */
class TextualListener: public BaseListener
{
public :

    /** Initialize progress bar
     *
     * _est_it: Defines the estimated number of iterations
     *
     * This method will initialize the progress bar with the number of estimated
     * iterations.
     */
    virtual void OnInitOperation(unsigned long _est_it);

    /** Show a textual bar with the progress in time
     *
     * When the input is negative then we are setting the progress bar, this
     * will be the total of elements to process. Afterwards the call to this
     * routine must be in ascending order, ie, 0, 1, 2, ... No. elements
     * The elapsed and estimation time is also shown in the output console.
     * _it: iteration number
     */
    virtual void OnProgress(unsigned long _it);

    /** Shows a message indicating the operation
     *
     * _rsop: string message
     */
    virtual void OnReportOperation(const std::string& _rsOp);
};

/** Shows a message and the time it was produced
 *
 * The format of the printing is "14:11:32 (12) => Hello, world", ie, The
 * message "Hello, world" was produced at 14:11:32 o'clock of the day 12. This
 * function needs not to read the time configuration (see time_config).
 *
 * @code
 * TimeMessage("Hello, world");
 * @endcode
 *
 * This function is not ported to Python.
 */
void TimeMessage(const std::string &message);
//@}

/** @name Little/Big endian
 *
 * These set of functions helps you to deal with the little/big endian
 * problems.
 */
//@{
/** Read from file
 *
 * This function is the same as fread from C, but at the end there is a flag
 * saying if data should be read in reverse order or not.
 *
 * @code
 * double f;
 * xmippFREAD(&f, sizeof(double), 1, fp, TRUE); // Reverse order
 *
 * double f;
 * xmippFREAD(&f, sizeof(double), 1, fp); // Normal order
 * @endcode
 *
 * This function is not ported to Python.
 */
size_t xmippFREAD(void* dest, size_t size, size_t nitems, FILE*& fp,
                  bool reverse = false);

/** Write to file
 *
 * This function is the same as fread from C, but at the end there is a flag
 * saying if data should be read in reverse order or not.
 *
 * This function is not ported to Python.
 */
size_t xmippFWRITE(const void* src,
                   size_t size,
                   size_t nitems,
                   FILE*& fp,
                   bool reverse = false);

/** Map file to memory.
 *
 * If size is less than 0, then the whole file is mapped, and the size
 * is correctly set to the file size.
 * */
void mapFile(const FileName &filename, char*&map,size_t &size, int &fileDescriptor, bool readOnly=true);

/** Unmap file*/
void unmapFile(char *&map, size_t &size, int& fileDescriptor);

/** Conversion little-big endian
 *
 * This function is not ported to Python.
 */
#define little22bigendian(x) swapbytes((unsigned char*)& x,sizeof(x))

/** Conversion little-big endian
 *
 * This function is not ported to Python.
 */
void swapbytes(char* v, unsigned long n);

/** Returns 1 if machine is big endian else 0
 */
bool IsBigEndian(void);

/** Returns 1 if machine is little endian else 0
 *  little-endian format (sometimes called the Intel format
 */
bool IsLittleEndian(void);
//@}

/** @name Marsaglia Marsaglia random functions
 *
 * These functions are not ported to Python.
 */
//@{
/** Marsaglia class.
 *
 * This class has been designed for those programs that  have to have a large
 * set of random numbers but do not have time to generate them properly on the
 * fly. The idea is to have  a large file with lots of random numbers in it and
 * to store some of them in memory and retrieve (from memory), as many times as
 * needed, in a random fashion.
 *
 * As source of the numbers I recommend the "Marsaglia's Random Number CDROM"
 * available, for free, at your favorite web site. (Search for it in any search
 * engine  and you will get tons of hits.) The following is from the CD
 * description:
 *
 * "This CDROM will contain 4.8 billion random bits. They were produced by a
 * combination of several of the best deterministic random number generators
 * (RNG's), together with three sources of white noise, as well as black noise
 * (from a rap music digital recording). My intent is to provide an unassailable
 * source for those who absolutely positively have to have a large, reliable set
 * of random numbers for serious simulation (Monte Carlo) studies."
 *
 * When the random numbers are doubles, by default they are in the interval [0,1]
 *
 * This class is not ported to Python.
 */
template <typename T>
class Marsaglia
{
private:
    char* random_vector; // read the data right here
    T* T_random_vector;
    long pointer_in_memory;
    FileName fn_rand;
    long vector_size;
    long Number_of_Numbers;

public:
    /** Constructor
     * @ingroup Marsaglia
     *
     * @code
     * Marsaglia rodalnino("masaglia", 1000000, 34);
     * @endcode
     *
     * M_max (optional) is the magnitude of the maximum value of the random
     * number (exclusive), therefore must be positive
     */
    Marsaglia(const FileName& fn_in, int No_Numbers)
    {
        Init(fn_in, No_Numbers);
    }

    /// Empty constructor
    Marsaglia()
    {}

    /** You may use init for reading another set of random numbers
     */
    void Init(const FileName& fn_in, int No_Numbers)
    {
        int Type_size; // sizeof(type)

        pointer_in_memory = 0;
        Number_of_Numbers = No_Numbers; // initialize class variable
        Type_size = sizeof(T);

        std::ifstream in(fn_in.c_str());
        in.seekg(0, std::ios::end); // End of file
        std::streampos sp = in.tellg(); // Size of file
        if (sp < Number_of_Numbers * Type_size)
            REPORT_ERROR(ERR_IO_SIZE, (std::string) "Marsaglia::Init: File " + fn_in +
                         "is too small");
        else
        {
            // get a random number to set the file pointer at a random position
            randomize_random_generator();

            random_vector = new char[(Number_of_Numbers * Type_size)];
            T_random_vector = (T*) random_vector;
            in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (double)(sp -
                                            (std::streamoff)(Number_of_Numbers * Type_size)))), std::ios::beg);
            in.read(random_vector, (Number_of_Numbers * Type_size));

            in.close();
        }
        if (typeid(double) == typeid(T))
            Verify_double();
    }

    /** Get a random number from the memory
     *
     * If you are at the end of the stream the pointer will be radomly moved
     * before stracting the number.
     */
    T Get_One_Number()
    {
        if (pointer_in_memory >= Number_of_Numbers)
            pointer_in_memory = (int) FLOOR(rnd_unif(0.f, (double)
                                            (Number_of_Numbers - 1)));
        return (T_random_vector[pointer_in_memory++]);
    }

    /** Calculate random vector log (use only with doubles)
     */
    void Marsaglia_log()
    {
        if (typeid(double) != typeid(T) && typeid(double) != typeid(T))
            REPORT_ERROR(ERR_TYPE_INCORRECT,
                         "Marsaglia: I do not know how to calculate integer logs");

        for (int hh = 0; hh < Number_of_Numbers; hh++)
            if (T_random_vector[hh] == 0.)
                T_random_vector[hh] = -1e+20f;
            else
                T_random_vector[hh] = log(T_random_vector[hh]);
    }

    /** Multiply random vector by constant
     */
    void mul(T mul_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] *= mul_cte;
    }

    /** Calculate mod of random vector, only make sense with integers
     */
    void operator&= (T mod_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] &= mod_cte;
    }

    /** Add a constant
     */
    void add(T add_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] += add_cte;
    }

    /** Set Maximum value (only valid for integers)
     */
    void M_max(const FileName& fn_in, T m_max)
    {
        int Type_size; // sizeof(type)
        Type_size = sizeof(T);

        std::ifstream in(fn_in.c_str());
        in.seekg(0, std::ios::end);              // End of file
        std::streampos sp = in.tellg();     // Size of file
        T power_of_2 = (T) NEXT_POWER_OF_2(m_max);
        if (power_of_2 == m_max)
            power_of_2 = (T) NEXT_POWER_OF_2(m_max + 1);
        T mask = power_of_2 - 1;
        T aux_number;
        m_max;

        // get a random number to set the file pointer at a random position
        in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (double)(sp -
                                        (std::streamoff)(Number_of_Numbers*Type_size)))), std::ios::beg);
        for (int ii = 0; ii < Number_of_Numbers;)
        {
            aux_number = T_random_vector[ii];
            aux_number &= mask;
            if (aux_number > m_max ||
                (T_random_vector[ii] <= 0) && (aux_number == 0))
            {
                if (in.eof())
                    in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (double)
                                                    (sp - (std::streamoff)(Number_of_Numbers*Type_size)))),
                             std::ios::beg);
                in.read((char*) &(T_random_vector[ii]), Type_size);
            }
            else
            {
                T_random_vector[ii] = aux_number * (T) SGN(T_random_vector[ii]);
                ii++;
            }
        }
        in.close();
    }

private:
    /** Verify double
     *
     * Be aware that Marsaglia reads blindly the data, therefore if the type
     * double is selected several of the "random" numbers may not be valid (the
     * number are created from a source of random bits and although 4 random
     * bits are one random integer, four random bits may not be a valid double).
     * If Marsaglia is "double" The constructor will run the following function
     * that will fix the problem
     */
    void Verify_double()
    {
        unsigned int* int_random_vector;
        long long MaxInteger;
        if (sizeof(double) != sizeof(int))
            REPORT_ERROR(ERR_TYPE_INCORRECT,
                         "Marsaglia: I do not know how to make the double correction");
        MaxInteger = (long long) pow(2.0, sizeof(unsigned int) * 8.0);
        int_random_vector = (unsigned int*) random_vector;
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] = (T)((double) int_random_vector[hh] /
                                      (double) MaxInteger);
    }
};
//@}
/// binary comparison of two files skipping first "offset" bytes
bool compareTwoFiles(const FileName &fn1, const FileName &fn2, size_t offset = 0);



//@}
#endif
