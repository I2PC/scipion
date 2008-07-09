
// class based in previous fftw wrapper by  Anna Kreshuk   07/4/2006

//////////////////////////////////////////////////////////////////////////
//                                                                      
// xmippFftw                                                       
//                                                                      
// An interface classes to the FFTW package
// Only the basic interface of FFTW is implemented.
//
// Computes a real input/complex output discrete Fourier transform in 1 or more
// dimensions. However, only out-of-place transforms are now supported for transforms
// in more than 1 dimension. For detailed information about the computed transforms,
// please refer to the FFTW manual
//
// How to use it:
// 1) Create an instance of xmippFftw - this will allocate input and output
//    arrays (unless an in-place transform is specified)
// 2) Run the Init() function with the desired flags and settings (see function
//    comments for possible kind parameters)
// 3) Set the data (via SetPoints()or SetPoint() functions)
// 4) Run the Transform() function
// 5) Get the output (via GetPoints() or GetPoint() functions)
// 6) Repeat steps 3)-5) as needed
// For a transform of the same size, but with different flags, 
// rerun the Init() function and continue with steps 3)-5)
//
// NOTE: 1) running Init() function will overwrite the input array! Don't set any data
//          before running the Init() function
//       2) FFTW computes unnormalized transform, so doing a transform followed by 
//          its inverse will lead to the original array scaled by the transform size
// 
//
//////////////////////////////////////////////////////////////////////////

#ifndef XMIPPFFTW
#define XMIPPFFTW

#include <string>
#include <iostream>
#include <complex.h> 
#include "../../external/fftw-3.1.2/api/fftw3.h"
#include <climits>
#include "matrix3d.h"

class xmippFftw{
 public:
   /** I grant direct access to the arrays, hope everything will be ok */
 
   /** input array */
   double     *fIn;      
   /**output array  */
   double     *fOut;  
   /* size of fft_output array */
   int sizeout ;
   /* total size of the transform */
   int     fTotalSize;  
       
   unsigned MapFlag(std::string flag);
 protected:
   /** delete input array when done */
   bool destroy_fIn; 
   /* fftw plan (the plan how to compute the transform) */
   void     *fPlan;      
   /*transform sizes in each dimension */
   int    *fN;         
   /** number of dimensions */
   int     fNdim;       
   /** transform flags */
   std::string fFlags;    
   /**sign of the exponent of the transform 
   (-1 is FFTW_FORWARD and +1 FFTW_BACKWARD) */
   int     fSign; 
   
   /** Use the input array to store the output data */
   bool inPlace;    
   
   /** wisdown file */
   std::string wisdom_name;

 public:
 /** Default constructor */
   xmippFftw();

  /** Initialization */
  void myxmippFftw(void);

   /**For 1d transforms Allocates memory for the input array, and, 
   if inPlace = false, for the output array. 
   If already_reserved != NULL no allocation is done

   Since we do not use SIND we are not required to use fftw_malloc..
   */ 
   xmippFftw(int n, bool inPlace, double * already_reserved=NULL);
   /** Initialization 1D constructor */

   void myxmippFftw(int n, bool my_inPlace, double * already_reserved);
      
   /* For ndim-dimensional transforms, Second argurment contains sizes of the
   transform in each dimension */
   xmippFftw(int ndim, int *n, bool inPlace=false, double * already_reserved=NULL);
   /** Initialization ND constructor */
   void myxmippFftw(int ndim, int *n, bool my_inPlace, 
                     double * already_reserved);
   /** Destroy... */
   ~xmippFftw();
   
   /** Destroys the data arrays and the plan. However, some plan 
   information stays around until the root session is over, 
   and is reused if other plans of the same size are created
   */

   void deleteXmippFftw(void);
   /**Creates the fftw-plan
   
   NOTE:  input and output arrays are overwritten during initialisation,
          so don't set any points, before running this function!!!!!
          unless "estimate" is used.
          
   Arguments  kind is dummy and not need to be specified
   Possible flag_options:
   "ES" (from "estimate") - no time in preparing the transform, but probably sub-optimal
      performanc
   "M" (from "measure") - some time spend in finding the optimal way to do the transform
   "P" (from "patient") - more time spend in finding the optimal way to do the transform
   "EX" (from "exhaustive") - the most optimal way is found
   This option should be chosen depending on how many transforms of the same size and
   type are going to be done. Planning is only done once, for the first transform of this
   size and type.
   
   flags is a boolean OR (`|') of zero or more of the following:

    * FFTW_MEASURE: this flag tells FFTW to find the optimal plan by actually
computing several FFTs and measuring their execution time. Depending on the
installation, this can take some time. (2)

    * FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan (for a
RISC processor with many registers). If neither FFTW_ESTIMATE nor FFTW_MEASURE
is provided, the default is FFTW_ESTIMATE.

    * FFTW_OUT_OF_PLACE: produce a plan assuming that the input and output
arrays will be distinct (this is the default).

    * FFTW_IN_PLACE: produce a plan assuming that you want the output in the
input array. The algorithm used is not necessarily in place: FFTW is able to
compute true in-place transforms only for small values of n. If FFTW is not
able to compute the transform in-place, it will allocate a temporary array
(unless you provide one yourself), compute the transform out of place, and copy
the result back. Warning: This option changes the meaning of some parameters of
fftw (see section Computing the One-dimensional Transform).

      The in-place option is mainly provided for people who want to write their
own in-place multi-dimensional Fourier transform, using FFTW as a base. For
example, consider a three-dimensional n * n * n transform. An out-of-place
algorithm will need another array (which may be huge). However, FFTW can
compute the in-place transform along each dimension using only a temporary
array of size n. Moreover, if FFTW happens to be able to compute the transform
truly in-place, no temporary array and no copying are needed. As distributed,
FFTW `knows' how to compute in-place transforms of size 1, 2, 3, 4, 5, 6, 7, 8,
9, 10, 11, 12, 13, 14, 15, 16, 32 and 64.

      The default mode of operation is FFTW_OUT_OF_PLACE.

    * FFTW_USE_WISDOM: use any wisdom that is available to help in the creation
of the plan. (See section 2.6 Words of Wisdom.) This can greatly speed the
creation of plans, especially with the FFTW_MEASURE option. FFTW_ESTIMATE plans
can also take advantage of wisdom to produce a more optimal plan (based on past
measurements) than the estimation heuristic would normally generate. When the
FFTW_MEASURE option is used, new wisdom will also be generated if the current
transform size is not completely understood by existing wisdom.

Algorithm-restriction flags

    * FFTW_DESTROY_INPUT specifies that an out-of-place transform is allowed to
overwrite its input array with arbitrary data; this can sometimes allow more
efficient algorithms to be employed.

    * FFTW_PRESERVE_INPUT specifies that an out-of-place transform must not
change its input array. This is ordinarily the default, except for c2r and hc2r
(i.e. complex-to-real) transforms for which FFTW_DESTROY_INPUT is the default.
In the latter cases, passing FFTW_PRESERVE_INPUT will attempt to use algorithms
that do not destroy the input, at the expense of worse performance; for
multi-dimensional c2r transforms, however, no input-preserving algorithms are
implemented and the planner will return NULL if one is requested.

    * FFTW_UNALIGNED specifies that the algorithm may not impose any unusual
alignment requirements on the input/output arrays (i.e. no SIMD may be used).
This flag is normally not necessary, since the planner automatically detects
misaligned arrays. The only use for this flag is if you want to use the guru
interface to execute a given plan on a different array that may not be aligned
like the original. (Using fftw_malloc makes this flag unnecessary even then.) 

By the way I have not find much any adventage using the wisdom mechanism or
The MEASURE flag 
       */
   void       Init( std::string flags, int sign, bool wisdom_flag=false);
   
   /**Computes the transform, specified in Init() function*/
   virtual void       Transform();

   /** Fills the array data with the computed transform.
       or the initial points.
       IMPORTANT note that fftw and xmipp use different access order.
       first coordiante is contigous in memory
       in fftw last coordinate is contigous in memory.
       This routine performs a blind copy and will not
       take care of this difference
       Only (roughly) a half of the transform is copied (exactly the output of FFTW),
       the rest being Hermitian symmetric with the first half
   */
   void       GetPoints(double * data, bool fromInput = false) const;

   /** Fills the array data with the computed transform as complex double.
       or the initial points.
       IMPORTANT note that fftw and xmipp use different access order.
       first coordiante is contigous in memory
       in fftw last coordinate is contigous in memory.
       This routine performs a blind copy and will not
       take care of this difference
       Only (roughly) a half of the transform is copied (exactly the output of FFTW),
       the rest being Hermitian symmetric with the first half
   */
   void       GetPoints(std::complex<double> * data, bool fromInput = false) const;

   /**Set all input points.
   The values are copied. If inverse points should be ordered as follows:
   [re_0, im_0, re_1, im_1, ..., re_n, im_n) in the input array
       IMPORTANT note that fftw and xmipp use different access order.
       first coordiante is contigous in memory
       in fftw last coordinate is contigous in memory.
       This routine performs a blind copy and will not
       take care of this difference
   */   
   void       SetPoints(const double *data);

   /**Set all input points from complex double array.
   The values are copied. If inverse points should be ordered as follows:
   [re_0, im_0, re_1, im_1, ..., re_n, im_n) in the input array
       IMPORTANT note that fftw and xmipp use different access order.
       first coordiante is contigous in memory
       in fftw last coordinate is contigous in memory.
       This routine performs a blind copy and will not
       take care of this difference
   */   
   void       SetPoints(const std::complex<double> *data);

   int      GetSize() const {return fTotalSize;}
   int     *GetN()    const {return fN;}
   int      GetNdim() const {return fNdim;}
   int      GetSign() const {return fSign;}
   std::string  GetTransformFlag() const {return fFlags;}
   bool     IsInplace() const {if (fOut) return true; else return false;};

   /** Normalize transform dividing by the number of points. Computing the
       forward transform followed by the backward transform (or vice versa) yields the
       original array scaled by the size of the array. (For multi-dimensional
       transforms, the size of the array is the product of the dimensions.) We could,
       instead, have chosen a normalization that would have returned the unscaled
       array. Or, to accomodate the many conventions in this matter, the transform
       routines could have accepted a "scale factor" parameter. We did not do this,
       however, for two reasons. First, we didn't want to sacrifice performance in the
       common case where the scale factor is 1. Second, in real applications the FFT
       is followed or preceded by some computation on the data, into which the scale
       factor can typically be absorbed at little or no cost.
       */
   void Normalize(void);
   /**
   
 This function instructs FFTW to
spend at most `seconds' seconds (approximately) in the planner. If `seconds ==
FFTW_NO_TIMELIMIT' (the default value, which is negative), then planning time
is unbounded. Otherwise, FFTW plans with a progressively wider range of
algorithms until the the given time limit is reached or the given range of
algorithms is explored, returning the best available plan. For example,
specifying `FFTW_PATIENT' first plans in `FFTW_ESTIMATE' mode, then in
`FFTW_MEASURE' mode, then finally (time permitting) in `FFTW_PATIENT'. If
`FFTW_EXHAUSTIVE' is specified instead, the planner will further progress to
`FFTW_EXHAUSTIVE' mode.std::cerr << " flags " << opt << "\n";
   
   */
   void TimeLimit(double seconds);

   /** Wrapper for Matrix 1D, calls to constructor
       with appropiate values. Does NOT call to plan
       or Transform */
   xmippFftw(Matrix1D<double> &img, bool my_inPlace=false, 
                                            bool already_reserved=true);
   /** Wrapper for Matrix 2D, calls to constructor
       with appropiate values and perform the Fourier transform 
       if init_and_do_transform = true. For 1 single image is OK for
       a sel file reusing old fftw structure may be faster
        
       If you have an array stored in column-major order (as xmipp does) and wish to
       transform it using FFTW, it is quite easy to do. When creating the object, simply
       pass the dimensions of the array to the planner in reverse order. For example,
       if your array is a rank three N x M x L matrix in column-major order, you
       should pass the dimensions of the array as if it were an L x M x N matrix
       (which it is, from the perspective of FFTW). Of course the FFT
       will be transposed
       Not tested for non square images
        */

    xmippFftw(Matrix2D<double> &img, bool init_and_do_transform=true);

   /** Wrapper for Matrix 3D, calls to constructor
       with appropiate values and perform the Fourier transform 
       if init_and_do_transform = true. For I single image is OK for
       a sel file reusing old aatw structure may be faster
        
        
       If you have an array stored in column-major order (as xmipp does) and wish to
       transform it using FFTW, it is quite easy to do. When creating the object, simply
       pass the dimensions of the array to the planner in reverse order. For example,
       if your array is a rank three N x M x L matrix in column-major order, you
       should pass the dimensions of the array as if it were an L x M x N matrix
       (which it is, from the perspective of FFTW). Of course the FFT
       will be transposed
       Not tested for non square images
        */

    xmippFftw(Matrix3D<double> &img, bool init_and_do_transform=true);
                                            
   /** Modify Real Data so the Fourier trasform is at the center
       CenterRealDataAfterTransform function must be applied after Fourier inversion
    */
   void CenterRealDataBeforeTransform(void);
    /** Modify Output Real Data so the Fourier transform is at the center
       CenterRealDataBeforeTransform function must be applied before Fourier inversion
    */
   void CenterRealDataAfterTransform(void);
  /** Delete fIn vector while keeping fOut. This is usefull if you
       want to save memory */
   void delete_fIn(void);

   /** Fill the fOut array with the frequency values (0-0.5) for each pixel */
   void getResolutionAllPoints(void);

   /** Applies a bandpass filter to an image (either  2 or 3 D). 
       frecuencies in range 0-0.5 */
   
   void img_bandpass_filter(double res_hi, double width);

/** Fourier-Ring-Correlation between two 2D-matrices using FFT
 * @ingroup FourierOperations
 */
    void fourier_ring_correlation(xmippFftw& fft_m2,
                              double sampling_rate,
                              Matrix1D< double >& freq,
                              Matrix1D< double >& frc,
                              Matrix1D< double >& frc_noise);
                              
/** Radial average for Fourier transforms*/
void fftwRadialAverage(double *AUX,
                                  Matrix1D< double >& radial_mean,
                                  Matrix1D< int >& radial_count,
                                  bool rounding =true );

/** Read Fourier transform from disk*/
void read(FileName fn);

/** Write Fourier transform to disk*/
void write(FileName fn);

};
#include <sys/stat.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h> 
#include <iostream>
#include <fstream>
class Lock
{
    std::string filename;
    public:
        /** If not looked this is true */
        bool non_locked;
    public:
        Lock  (std::string basefilename);
        ~Lock (void);

};

#endif
