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

#ifndef FUNCS_H
#define FUNCS_H

// TODO remove this, never open namespaces in header files
using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <climits>
#include <algorithm>
#include <vector>

// For timing functions
#include <unistd.h>
#include <sys/times.h>
#ifdef _IRIX65
#include <sys/types.h>
#include <time.h>
#endif

#include "numerical_recipes.h"
#include "macros.h"
#include "error.h"

// FIXME local paths in includes is not recommended. The file could be removed.
#include "../../xmippConfiguration.inc"

/// @defgroup GeneralFunctions General functions
/// @ingroup DataLibrary

/// @defgroup NumericalFunctions Numerical functions
/// @ingroup GeneralFunctions

/** Tabulated Sinc = SIN(PI*X)/(PI*X)
 * @ingroup NumericalFunctions
 *
 * A lookup-table with the given sampling rate and range is created.
 *
 * @code
 * tabsinc TSINC(0.0001,64);
 * if (TSINC(1)==0)
 *    cout << "This is true!\n";
 * @endcode
 */
class tabsinc
{
private:

    float sampl;
    int xmax;
    int no_elem;
    float* tabulatedsinc;

public:
    /** Constructor with sampling rate and range */
    tabsinc(const float dd, const int xx)
    {
        sampl = dd;
        xmax = xx;
        filltable();
    }

    // Destructor
    ~tabsinc()
    {
        free(tabulatedsinc);
    }

    /** Value access. Tabulated sine in radians */
    float operator()(float val) const
    {
        return tabulatedsinc[ABS((int)(val / sampl))];
    }

    /** Actually fill the table */
    void filltable()
    {
        no_elem = (int)(xmax / sampl);
        tabulatedsinc = (float*) malloc(no_elem * sizeof(float));
        tabulatedsinc[0] = 1;
        for (int i = 1; i < no_elem; i++)
        {
            float xx = (float) i * sampl * PI;
            tabulatedsinc[i] = sin(xx) / xx;
        }
    }
};

/** Kaiser-Bessel function
 * @ingroup NumericalFunctions
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
 *  for which the I0 window is non-zero.  Setting "vtable"
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
    int K; /** I0 window size */
    double vtable; /** table I0 non-zero domain maximum */
    int ntable;
    vector<double> i0table;
    double dtable; /** table spacing */
    double alphar; /** alpha*r */
    double fac; /** 2*pi*alpha*r*v */
    double vadjust; 
    double facadj; /** 2*pi*alpha*r*vadjust */
    virtual void build_I0table(); /** Tabulate I0 window for speed */
    double fltb;

public:
    /** Empty constructor */
    KaiserBessel() { }

    /** Constructor with parameters */
    KaiserBessel(double alpha_, int K, double r_,
		 double v_, int N_, double vtable_=0., 
		 int ntable_ = 5999);

    /** Destructor */
    virtual ~KaiserBessel() {};

    /** Compute the maximum error in the table */
    double I0table_maxerror();
    vector<double> dump_table() {
	return i0table;
    }

    /** Kaiser-Bessel Sinh window function */
    virtual double sinhwin(double x) const;

    /** Kaiser-Bessel I0 window function */
    virtual double i0win(double x) const;

    /** Kaiser-Bessel I0 window function (uses table lookup) */
    inline double i0win_tab(double x) const {
	double xt;
	if(x<0.) xt = -x*fltb+0.5; else xt = x*fltb+0.5;
	return i0table[ (int) xt];
    }

    /** Return the size of the I0 window */
    int get_window_size() const { return K; }

};

/** Solve second degree equation
 * @ingroup NumericalFunctions
 *
 *  ax^2+bx+c=0
 *
 * It returns the number of real solutions, 0 if the two roots are complex
 * and -1 if the equation is impossible to satisfy (Ex: 7=0). A number is said
 * to be 0 if its absolute magnitude is smaller than precision. This is used to
 * avoid dividing by 0
 */
int solve_2nd_degree_eq(float a,
                        float b,
                        float c,
                        float& x1,
                        float& x2,
                        float prec = XMIPP_EQUAL_ACCURACY);

/** 1D gaussian value
 * @ingroup NumericalFunctions
 *
 * This function returns the value of a univariate gaussian function at the
 * point x.
 */
double gaussian1D(double x, double sigma, double mu = 0);

/** 2D gaussian value
 * @ingroup NumericalFunctions
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

/// @defgroup MiscellaneousFunctions Miscellaneous functions
/// @ingroup GeneralFunctions

/** Print a boolean value
 * @ingroup MiscellaneousFunctions
 * */
void print(ostream& o, const bool b);

/** Print a value in binary
 * @ingroup MiscellaneousFunctions
 *
 * So far not instatiate for float/double number
 */
template <typename T>
void printb(ostream& o, T value)
{
    char buf[CHAR_BIT * sizeof(T) + 1];
    size_t i;

    for (i = 0; i < CHAR_BIT * sizeof(T); ++i)
    {
        buf[i] = '0' + (value & 1);
        value >>= 1;
    }
    buf[i] = 0;

    o << buf;
}

/// @defgroup ComplexFunctions Complex functions
/// @ingroup GeneralFunctions

/** Real/Imaginary --> Complex
 * @ingroup ComplexFunctions
 *
 * The output array(s) must be already resized.
 */
template <typename T>
void RealImag2Complex(const T* _real,
                      const T* _imag,
                      complex< double >* _complex,
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
 * @ingroup ComplexFunctions
 *
 * The output array(s) must be already resized.
 */
template <typename T>
void AmplPhase2Complex(const T* _ampl,
                       const T* _phase,
                       complex< double >* _complex,
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
 * @ingroup ComplexFunctions
 *
 * The output array(s) must be already resized.
 */
template <typename T>
void Complex2RealImag(const complex< double >* _complex,
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
 * @ingroup ComplexFunctions
 *
 * The output array(s) must be already resized.
 */
template <typename T>
void Complex2AmplPhase(const complex< double >* _complex,
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

/** @defgroup RandomFunctions Random functions
 * @ingroup GeneralFunctions
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
 *     cout << rnd_unif(-1,1) << endl;
 * @endcode
 */

/** Reset uniform random generator to a known point
 * @ingroup RandomFunctions
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

/** Reset uniform random generator according to the clock
 * @ingroup RandomFunctions
 *
 * This time the initialisation itself assures a random sequence different each
 * time the program is run. Be careful not to run the program twice within the
 * same second as the initialisation will be the same for both runs.
 *
 * @code
 * rnd_rnd_unif();
 * @endcode
 */
void randomize_random_generator();

/** Produce a uniform random number between 0 and 1
 * @ingroup RandomFunctions
 *
 * @code
 * cout << "This random number should be between 0 and 1: " << rnd_unif()
 * << endl;
 * @endcode
 */
float rnd_unif();

/** Produce a uniform random number between a and b
 * @ingroup RandomFunctions
 *
 * @code
 * cout << "This random number should be between 0 and 10: " << rnd_unif(0,10)
 * << endl;
 * @endcode
 */
float rnd_unif(float a, float b);

/** Produce a gaussian random number with mean 0 and standard deviation 1
 * @ingroup RandomFunctions
 *
 * @code
 * cout << "This random number should follow N(0,1): " << rnd_gaus()
 * << endl;
 * @endcode
 */
float rnd_gaus();

/** Produce a gaussian random number with mean a and standard deviation b
 * @ingroup RandomFunctions
 *
 * @code
 * cout << "This random number should follow N(1,4): " << rnd_gaus(1,2)
 * << endl;
 * @endcode
 */
float rnd_gaus(float a, float b);

/** Gaussian area from -x0 to x0
 * @ingroup RandomFunctions
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * x0 must be positive
 */
float gaus_within_x0(float x0, float mean = 0, float stddev = 1);

/** Gaussian area outisde -x0 to x0
 * @ingroup RandomFunctions
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * x0 must be positive
 */
float gaus_outside_x0(float x0, float mean = 0, float stddev = 1);

/** Gaussian area from -inf to x0
 * @ingroup RandomFunctions
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * There is no restriction over the sign of x0
 */
float gaus_up_to_x0(float x0, float mean = 0, float stddev = 1);

/** Gaussian area from x0 to inf
 * @ingroup RandomFunctions
 *
 * By default the gaussian mean is 0 and the gaussian standard deviation is 1.
 * There is no restriction over the sign of x0
 */
float gaus_from_x0(float x0, float mean = 0, float stddev = 1);

/** t0 for a given two-sided probability
 * @ingroup RandomFunctions
 *
 * This function returns t0 such that the student probability outside t0 is
 * equal to p
 */
float student_outside_probb(float p, float degrees_of_freedom);

/** student area from -t0 to t0
 * @ingroup RandomFunctions
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * t0 must be positive
 */
float student_within_t0(float t0, float degrees_of_freedom);

/** student area outisde -t0 to t0
 * @ingroup RandomFunctions
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * t0 must be positive
 */
float student_outside_t0(float t0, float degrees_of_freedom);

/** student area from -inf to t0
 * @ingroup RandomFunctions
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * There is no restriction over the sign of t0
 */
float student_up_to_t0(float t0, float degrees_of_freedom);

/** student area from t0 to inf
 * @ingroup RandomFunctions
 *
 * By default the student mean is 0 and the student standard deviation is 1.
 * There is no restriction over the sign of t0
 */
float student_from_t0(float t0, float degrees_of_freedom);

/** chi2 area from -inf to t0
 * @ingroup RandomFunctions
 *
 * By default the chi2 mean is 0 and the chi2 standard deviation is 1.
 * There is no restriction over the sign of t0
 */
float chi2_up_to_t0(float t0, float degrees_of_freedom);

/** chi2 area from t0 to inf
 * @ingroup RandomFunctions
 *
 * By default the chi2 mean is 0 and the chi2 standard deviation is 1.
 * There is no restriction over the sign of t0
 */
float chi2_from_t0(float t0, float degrees_of_freedom);

/** Produce a log uniform random number between a and b
 * @ingroup RandomFunctions
 *
 * Watch out that the following inequation must hold 0<a<=b.
 *
 * @code
 * cout << "This random number should be between 1 and 1000: "
 * << rnd_log(10,1000)<< endl;
 * @endcode
 */
float rnd_log(float a, float b);

/** @defgroup Filenames Filename handling
 * @ingroup GeneralFunctions
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

/** Filenames
 * @ingroup Filenames
 *
 * This class allows you a lot of usual and common manipulations with filenames.
 * See filename conventions for a detailed explanation of the Filenames dealed
 * here, although most of the functions work with the more general model
 * "name.extension"
 */
class FileName: public string
{
public:
    /// @defgroup FilenameConstructors Filename constructors
    /// @ingroup Filenames

    /** Empty constructor
     * @ingroup FilenameConstructors
     *
     * The empty constructor is inherited from the string class, so an empty
     * FileName is equal to "".
     *
     * @code
     * FileName fn_blobs;
     * @endcode
     */
    FileName(): string("")
    {}

    /** Constructor from string
     * @ingroup FilenameConstructors
     *
     * The constructor from a string allows building complex expressions based
     * on the string class. Notice that in the following example the type
     * casting to string is very important, if not, the operation is just a
     * pointer movement instead of a string concatenation.
     *
     * @code
     * FileName fn_blobs((string) "art00001" + ".blobs");
     * @endcode
     */
    FileName(const string& str): string(str)
    {}

    /** Constructor from char*
     * @ingroup FilenameConstructors
     */
    FileName(const char* str) : string(str)
    {}

    /** Copy constructor
     * @ingroup FilenameConstructors
     */
    FileName(const FileName& fn) : string(fn)
    {}

    // FIXME What is this?
#ifdef NEVER_DEFINED
    /** Constructor from root, number and extension
     * @ingroup FilenameConstructors
     *
     * The number and extension are optional.
     *
     * @code
     * FileName fn_proj("g1ta00001.xmp"); // fn_proj = "g1ta00001.xmp"
     * FileName fn_proj("g1ta",1,"xmp"); // fn_proj = "g1ta00001.xmp"
     * FileName fn_proj("g1ta",1); // fn_proj = "g1ta00001"
     * @endcode
     */
    FileName(const char* str, int no = -1, const string& ext = "")
    {
        compose(str, no, ext);
    }
#endif

    /** Constructor from root and extension
     * @ingroup FilenameConstructors
     *
     * None of the parameters is optional
     *
     * @code
     * FileName fn_proj("g1ta00001", "xmp"); // fn_proj = "g1ta00001.xmp"
     * @endcode
     */
    FileName(const char* str, const string& ext): string(str + ext)
    {}

    /// @defgroup FilenameComposing Composing/Decomposing the filename
    /// @ingroup Filenames

    /** Compose from root, number and extension
     * @ingroup FilenameComposing
     *
     * @code
     * fn_proj.compose("g1ta", 1, "xmp");  // fn_proj = "g1ta00001.xmp"
     * @endcode
     */
    void compose(const string& str, int no, const string& ext);

    /** Get the root from a filename
     * @ingroup FilenameComposing
     *
     * @code
     * FileName fn_root, fn_proj("g1ta00001");
     * fn_root = fn_proj.get_root();
     * @endcode
     */
    FileName get_root() const;

    /** Get the base name from a filename
     * @ingroup FilenameComposing
     */
    string get_baseName() const;

    /** Get the number from a filename
     * @ingroup FilenameComposing
     *
     * If there is no number a -1 is returned.
     *
     * @code
     * FileName proj("g1ta00001");
     * int num = proj.get_number();
     * @endcode
     */
    int get_number() const;

    /** Get the last extension from filename
     * @ingroup FilenameComposing
     *
     * The extension is returned without the dot. If there is no extension "" is
     * returned.
     *
     * @code
     * string ext = fn_proj.get_extension();
     * @endcode
     */
    string get_extension() const;

    /** Random name
     * @ingroup FilenameComposing
     *
     * Generate a random name of the desired length.
     */
    void init_random(int length);

    /// @defgroup FilenameManipulators Filename manipulators
    /// @ingroup Filenames

    /** Add string at the beginning
     * @ingroup FilenameManipulators
     *
     * If there is a path then the prefix is added after the path.
     *
     * @code
     * fn_proj = "imgs/g1ta00001";
     * fn_proj.add_prefix("h"); // fn_proj == "imgs/hg1ta00001"
     *
     * fn_proj = "g1ta00001";
     * fn_proj.add_prefix("h"); // fn_proj == "hg1ta00001"
     * @endcode
     */
    FileName add_prefix(const string& prefix) const;

    /** Add extension at the end.
     * @ingroup FilenameManipulators
     *
     * The "." is added. If teh input extension is "" then the same name is
     * returned, with nothing added.
     *
     * @code
     * fn_proj = "g1ta00001";
     * fn_proj.add_extension("xmp"); // fn_proj == "g1ta00001.xmp"
     * @endcode
     */
    FileName add_extension(const string& ext) const;

    /** Remove last extension, if any
     * @ingroup FilenameManipulators
     *
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.without_extension(); // fn_proj == "g1ta00001"
     *
     * fn_proj = "g1ta00001";
     * fn_proj = fn_proj.without_extension(); // fn_proj == "g1ta00001"
     * @endcode
     */
    FileName without_extension() const;

    /** Remove the root
     * @ingroup FilenameManipulators
     *
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.without_root(); // fn_proj == "00001.xmp"
     *
     * fn_proj = "g1ta00001";
     * fn_proj = fn_proj.without_root(); // fn_proj == "00001"
     * @endcode
     */
    FileName without_root() const;

    /** Insert before first extension
     * @ingroup FilenameManipulators
     *
     * If there is no extension, the insertion is performed at the end.
     *
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.insert_before_extension("pp");
     * // fn_proj == "g1ta00001pp.xmp"
     *
     * fn_proj = "g1ta00001";
     * fn_proj = fn_proj.insert_before_extension("pp");
     * // fn_proj=="g1ta00001pp"
     * @endcode
     */
    FileName insert_before_extension(const string& str) const;

    /** Remove a certain extension
     * @ingroup FilenameManipulators
     *
     * It doesn't matter if there are several extensions and the one to be
     * removed is in the middle. If the given extension is not present in the
     * filename nothing is done.
     *
     * @code
     * fn_proj = "g1ta00001.xmp.bak";
     * fn_proj = fn_proj.remove_extension("xmp");
     * // fn_proj == "g1ta00001.bak"
     * @endcode
     */
    FileName remove_extension(const string& ext) const;

    /** Remove all extensions
     * @ingroup FilenameManipulators
     */
    FileName remove_all_extensions() const;

    /** Substitute ext1 by ext2
     * @ingroup FilenameManipulators
     *
     * It doesn't matter if ext1 is in the middle of several extensions. If ext1
     * is not present in the filename nothing is done.
     *
     * @code
     * fn_proj = "g1ta00001.xmp.bak";
     * fn_proj = fn_proj.substitute_extension("xmp", "bor");
     * // fn_proj == "g1ta00001.bor.bak"
     *
     * fn_proj = "g1ta00001.xmp.bak";
     * fn_proj = fn_proj.substitute_extension("tob", "bor");
     * // fn_proj=="g1ta00001.xmp.bak"
     * @endcode
     */
    FileName substitute_extension(const string& ext1, const string& ext2) const;

    /** Without a substring
     * @ingroup FilenameManipulators
     *
     * If the substring is not present the same FileName is returned, if it is
     * there the substring is removed.
     */
    FileName without(const string& str) const;

    /** Remove until prefix
     * @ingroup FilenameManipulators
     *
     * Remove the starting string until the given prefix, inclusively. For
     * instance /usr/local/data/ctf-image00001.fft with ctf- yields
     * image00001.fft. If the prefix is not found nothing is done.
     */
    FileName remove_until_prefix(const string& str) const;

    /** Remove all directories
     * @ingroup FilenameManipulators
     *
     * Or if keep>0, then keep the lowest keep directories
     */
    FileName remove_directories(int keep = 0) const;
};

/** This class is used for comparin filenames.
 * @ingroup Filenames
 *
 * Example: "g0ta00001.xmp" is less than "g0ta00002.xmp"
 *
 * This class is needed to define a std::map as
 * map<FileName,FileName,FileNameComparison> myMap;
 */
class FileNameComparison {
public:
    bool operator ()(const FileName &fn1, const FileName &fn2)
    {
    	return fn1<fn2;
    }
};

/** True if the file exists in the current directory
 * @ingroup Filenames
 *
 * @code
 * if (exists("g1ta00001"))
 *     cout << "The file exists" << endl;
 * @endcode
 */
int exists(const FileName& fn);

/** This function raised an ERROR if the filename if not empty and if
 * the corresponding file does not exist. 
 * This may be useful to have a better (killing) control on (mpi-regulated) jobs
 * @ingroup Filenames
 *
 * @code
 *   exit_if_not_exists("control_file.txt");
 * @endcode
 */
void exit_if_not_exists(const FileName &fn);

/** Waits until the given filename has a stable size
 * @ingroup Filenames
 *
 * The stable size is defined as having the same size within two samples
 * separated by time_step (microsecs).
 *
 * An exception is throw if the file exists but its size cannot be stated.
 */
void wait_until_stable_size(const FileName& fn,
                            unsigned long time_step = 250000);

/** Write a zero filled file with the desired size.
 * @ingroup Filenames
 *
 * The file is written by blocks to speed up, you can modify the block size.
 * An exception is thrown if any error happens
 */
void create_empty_file(const FileName& fn,
                       unsigned long long size,
                       unsigned long long block_size = 102400);

/** Returns the base directory of the Xmipp installation
 * @ingroup Filenames
 */
FileName xmippBaseDir();

/** @defgroup TimeManaging Time managing
 * @ingroup GeneralFunctions
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
 * annotate_time(&t0);
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
 * TimeStamp t0;
 * float to_go;
 *
 * // Beginning of the program
 * time_config();
 * ...
 *
 * annotate_time(&t0);
 * // Part to be measured
 * for (int i=0; i<60; i++)
 * {
 *     ...
 *     // Compute the time to go with the fraction of work already done
 *     to_go = time_to_go(t0, (float) (i + 1) / 60);
 *     cout << "I think you will be here " << to_go << "seconds more\n";
 * }
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
 */

#ifdef _NO_TIME
typedef int TimeStamp; // Any other kind of data will do
#else
typedef struct tms TimeStamp; // Renaming of the time structure
#endif

/** Read the system clock frequency
 * @ingroup TimeManaging
 *
 * This operation is needed only once in a program always we want to have a time
 * measure, or an estimation of remaining time.
 *
 * @code
 * time_config();
 * @endcode
 */
void time_config();

/** Annotate actual time
 * @ingroup TimeManaging
 *
 * This annotation is used later to compute the elapsed time.
 *
 * @code
 * TimeStamp t0;
 * annotate_time(&t0);
 * @endcode
 */
void annotate_time(TimeStamp* time);

/** Acumulate time
 * @ingroup TimeManaging
 *
 * Initially dest_time should be set to orig time. Then you acumulate succesive
 * times calling this function (Destination time=destination_time + (now -
 * original time)) and finally the elapsed time is the dest time minus the first
 * one (the one which initiliazed the dest time.
 */
void acum_time(TimeStamp* orig, TimeStamp* dest);

/** Compute elapsed time since a given annotation
 * @ingroup TimeManaging
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
 * float elapsed = elapsed_time(t0);
 *
 * TimeStamp t0;
 * annotate_time(&t0);
 * ...;
 * float elapsed = elapsed_time(t0, FALSE);
 * @endcode
 */
float elapsed_time(TimeStamp& time, bool _IN_SECS = true);

/** Show on screen the elapsed time since a given annotation
 * @ingroup TimeManaging
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
 */
void print_elapsed_time(TimeStamp& time, bool _IN_SECS = true);

/** Returns the estimated time left to finish
 * @ingroup TimeManagement
 *
 * To make this estimation the starting time must have been annotated before and
 * the fraction of the total amount of work must be estimated by the programmer.
 * See Time managing for an example.
 */
float time_to_go(TimeStamp& time, float fraction_done);

/** Initialise the progress bar
 * @ingroup TimeManagement
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
#define init_progress_bar(total) progress_bar(-(total));

/** Update progress bar
 * @ingroup TimeManaging
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

/** xmippBaseListener
 * @ingroup TimeManaging
 *
 * This class implements the xmipp listener class for notification of progress
 * status and other operations. It is an abstract class that contains base
 * functions useful for implementation of customer-design notification classes.
 */
class xmippBaseListener
{
public :
    /** Default constructor */
    xmippBaseListener(): verbosity(0), cancel(false)
    {}

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
    virtual void OnReportOperation(const string& _rsOp) = 0;

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


/** xmippTextualListener
 * @ingroup TimeManaging
 *
 * This class implements the xmipp classical textual listener class for
 * notification of progress status. It inherits from xmippBaseListener.
 */
class xmippTextualListener: public xmippBaseListener
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
    virtual void OnReportOperation(const string& _rsOp);
};

/** Shows a message and the time it was produced
 * @ingroup TimeManaging
 *
 * The format of the printing is "14:11:32 (12) => Hello, world", ie, The
 * message "Hello, world" was produced at 14:11:32 o'clock of the day 12. This
 * function needs not to read the time configuration (see time_config).
 *
 * @code
 * TimeMessage("Hello, world");
 * @endcode
 */
void TimeMessage(string message);

/** @defgroup LittleBigEndian Little/Big endian
 * @ingroup GeneralFunctions
 *
 * These set of functions helps you to deal with the little/big endian
 * problems.
 */

/** Read from file
 * @ingroup LittleBigEndian
 *
 * This function is the same as fread from C, but at the end there is a flag
 * saying if data should be read in reverse order or not.
 *
 * @code
 * float f;
 * FREAD(&f, sizeof(float), 1, fp, TRUE); // Reverse order
 *
 * float f;
 * FREAD(&f, sizeof(float), 1, fp); // Normal order
 * @endcode
 */
size_t FREAD(void* dest, size_t size, size_t nitems, FILE*& fp,
             bool reverse = false);

/** Write to file
 * @ingroup LittleBigEndian
 *
 * This function is the same as fread from C, but at the end there is a flag
 * saying if data should be read in reverse order or not.
 */
size_t FWRITE(const void* src,
              size_t size,
              size_t nitems,
              FILE*& fp,
              bool reverse = false);

/** Conversion little-big endian
 * @ingroup LittleBigEndian
 */
#define little22bigendian(x) ByteSwap((unsigned char*)& x,sizeof(x))

void ByteSwap(unsigned char* b, int n);

/** Returns 1 if machine is big endian else 0
 * @ingroup LittleBigEndian
 */
int IsBigEndian(void);

/** Returns 1 if machine is little endian else 0
 * @ingroup LittleBigEndian
 *  little-endian format (sometimes called the Intel format
 */
int IsLittleEndian(void);

/** @defgroup Marsaglia Marsaglia random functions
 * @ingroup RandomFunctions
 */

/** Marsaglia
 * @ingroup Marsaglia
 *
 * This class has been designed for those programs that  have to have a large
 * set of random numbers but do not have time to generate them properly on the
 * fly. The idea is to have  a large file with lots of random numbers in it and
 * to store some of them in memory and retrieve (from memory), as many times as
 * needed, in a random fashion.
 *
 * As source of the numbers I recomend the "Marsaglia's Random Number CDROM"
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
 * When the random numbers are floats, by default they are in the interval [0,1]
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

    // TODO Document
    Marsaglia()
    {}

    /** You may use init for reading another set of random numbers
     * @ingroup Marsaglia
     */
    void Init(const FileName& fn_in, int No_Numbers)
    {
        int Type_size; // sizeof(type)

        pointer_in_memory = 0;
        Number_of_Numbers = No_Numbers; // initialize class variable
        Type_size = sizeof(T);

        ifstream in(fn_in.c_str());
        in.seekg(0, ios::end); // End of file
        std::streampos sp = in.tellg(); // Size of file
        if (sp < Number_of_Numbers * Type_size)
            REPORT_ERROR(1, (string) "Marsaglia::Init: File " + fn_in +
                         "is too small");
        else
        {
            // get a random number to set the file pointer at a random position
            randomize_random_generator();

            random_vector = new char[(Number_of_Numbers * Type_size)];
            T_random_vector = (T*) random_vector;
            in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (float)(sp -
                                            (std::streamoff)(Number_of_Numbers * Type_size)))), ios::beg);
            in.read(random_vector, (Number_of_Numbers * Type_size));

            in.close();
        }
        if (typeid(float) == typeid(T))
            Verify_float();
    }

    /** Get a random number from the memory
     * @ingroup Marsaglia
     *
     * If you are at the end of the stream the pointer will be radomly moved
     * before stracting the number.
     */
    T Get_One_Number()
    {
        if (pointer_in_memory >= Number_of_Numbers)
            pointer_in_memory = (int) FLOOR(rnd_unif(0.f, (float)
                                            (Number_of_Numbers - 1)));
        return (T_random_vector[pointer_in_memory++]);
    }

    /** Calculate random vector log (use only with flloats)
     * @ingroup Marsaglia
     */
    void Marsaglia_log()
    {
        if (typeid(float) != typeid(T) && typeid(double) != typeid(T))
            REPORT_ERROR(1,
                         "Marsaglia: I do not know how to calculate integer logs");

        for (int hh = 0; hh < Number_of_Numbers; hh++)
            if (T_random_vector[hh] == 0.)
                T_random_vector[hh] = -1e+20f;
            else
                T_random_vector[hh] = log(T_random_vector[hh]);
    }

    /** Multiply random vector by constant
     * @ingroup Marsaglia
     */
    void mul(T mul_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] *= mul_cte;
    }

    /** Calculate mod of random vector, only make sense with integers
     * @ingroup Marsaglia
     */
    void operator&= (T mod_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] &= mod_cte;
    }

    /** Add a constant
     * @ingroup Marsaglia
     */
    void add(T add_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] += add_cte;
    }

    /** Set Maximun value (only valid for integers)
     * @ingroup Marsaglia
     */
    void M_max(const FileName& fn_in, T m_max)
    {
        int Type_size; // sizeof(type)
        Type_size = sizeof(T);

        ifstream in(fn_in.c_str());
        in.seekg(0, ios::end);              // End of file
        std::streampos sp = in.tellg();     // Size of file
        T power_of_2 = (T) NEXT_POWER_OF_2(m_max);
        if (power_of_2 == m_max)
            power_of_2 = (T) NEXT_POWER_OF_2(m_max + 1);
        T mask = power_of_2 - 1;
        T aux_number;
        m_max;

        // get a random number to set the file pointer at a random position
        in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (float)(sp -
                                        (std::streamoff)(Number_of_Numbers*Type_size)))), ios::beg);
        for (int ii = 0; ii < Number_of_Numbers;)
        {
            aux_number = T_random_vector[ii];
            aux_number &= mask;
            if (aux_number > m_max ||
                (T_random_vector[ii] <= 0) && (aux_number == 0))
            {
                if (in.eof())
                    in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (float)
                                                    (sp - (std::streamoff)(Number_of_Numbers*Type_size)))),
                             ios::beg);
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
    /** Verify float
     * @ingroup Marsaglia
     *
     * Be aware that Marsaglia reads blindly the data, therefore if the type
     * float is selected several of the "random" numbers may not be valid (the
     * number are created from a source of random bits and although 4 random
     * bits are one random integer, four random bits may not be a valid float).
     * If Marsaglia is "float" The constructor will run the following function
     * that will fix the problem
     */
    void Verify_float()
    {
        unsigned int* int_random_vector;
        long long MaxInteger;
        if (sizeof(float) != sizeof(int))
            REPORT_ERROR(1,
                         "Marsaglia: I do not know how to make the float correction");
        MaxInteger = (long long) pow(2.0, sizeof(unsigned int) * 8.0);
        int_random_vector = (unsigned int*) random_vector;
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] = (T)((double) int_random_vector[hh] /
                                      (double) MaxInteger);
    }
};

#endif
