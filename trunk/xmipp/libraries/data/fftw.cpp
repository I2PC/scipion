
//////////////////////////////////////////////////////////////////////////
//                                                                      
// xmippFftw                                                       
//                                                                      
//  interface classes to the FFTW package. Only the basic interface of FFTW is implemented.
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

#include "fftw.h"


//_____________________________________________________________________________
xmippFftw::xmippFftw()
{
//default constructor
myxmippFftw();

}
void xmippFftw::myxmippFftw(void)
{
//almost the default constructor

   fIn   = NULL;
   fOut  = NULL;
   fPlan = NULL;
   fN    = NULL;
   sizeout=0;
   destroy_fIn=true;   
   fTotalSize=1;

   wisdom_name="/tmp/myfftw";
}

//_____________________________________________________________________________
xmippFftw::xmippFftw(int n, bool my_inPlace, double * already_reserved)
{
   myxmippFftw();
   myxmippFftw(n, my_inPlace, already_reserved);
}

void xmippFftw::myxmippFftw(int n, bool my_inPlace, double * already_reserved)
{
//For 1d transforms
//Allocates memory for the input array, and, if inPlace = false, 
//for the output array. If already_reserved != NULL no allocation is
//done

/**
You are not required to use fftw_malloc. You can allocate your data in any way
that you like, from malloc to new (in C++) to a fixed-size array declaration. If
the array happens not to be properly aligned, FFTW will not use the SIMD
extensions.
*/ // we do not use sind extensions...
   //check for overflows, maximum allocation size cannot be greater than INT_MAX
   double check_size;
   check_size = (double) n * sizeof(double);
   if (check_size > (double)INT_MAX)
      {
        std::cout << "Error allocating memory." << std::endl;
        std::cout << "size requested" 
                  << check_size  
                  << " maximun size "
                  <<  INT_MAX
                  << std::endl;
        exit(1);
      }    
       
   inPlace = my_inPlace;
   sizeout=n;
   if (!inPlace){
      //fIn = fftw_malloc(sizeof(double)*n);
      try
      {
          fOut = new double [2*(n/2+1)];
      }
      catch (std::bad_alloc&)
      {
        std::cout << "Error allocating memory." << std::endl;
        exit(1);
      }    
   } 
   if (already_reserved == NULL){
      try
      {
          fIn  = new double [2*(n/2+1)];
      }
      catch (std::bad_alloc&)
      {
        std::cout << "Error allocating memory." << std::endl;
        exit(1);
      }    
   } else {
      destroy_fIn=false;
      fIn  = already_reserved;
   }

   fN = new int [1];
   fN[0] = n;
   fTotalSize = n;
   fNdim = 1;
}

//_____________________________________________________________________________
xmippFftw::xmippFftw(int ndim, int *n, bool my_inPlace, 
                     double * already_reserved)
{
   myxmippFftw();
   myxmippFftw(ndim, n, my_inPlace, already_reserved);
}


void xmippFftw::myxmippFftw(int ndim, int *n, bool my_inPlace, 
                     double * already_reserved)
{
//For ndim-dimensional transforms
//Second argurment contains sizes of the transform in each dimension

   inPlace = my_inPlace;

   if (ndim>1 && inPlace==true){
      std::cerr << "xmippFftw, multidimensional in-place r2c transforms are not implemented yet"
                << std::endl;
      exit(1);
   }
   fNdim = ndim;
   fTotalSize = 1;
   fN = new int[fNdim];

   //check for overflows, maximum allocation size cannot be greater than INT_MAX
   double check_size;
   check_size=1.;
   for (int i=0; i<fNdim; i++){
      fN[i] = n[i];
      check_size*= (double)n[i];
   }
   check_size *= (double) sizeof(double);
   if (check_size > (double)INT_MAX)
      {
        std::cout << "Error allocating memory." << std::endl;
        std::cout << "size requested" 
                  << check_size  
                  << " maximun size "
                  <<  INT_MAX
                  << std::endl;
        exit(1);
      }
      
   //compute vector size
   for (int i=0; i<fNdim; i++){
      fN[i] = n[i];
      fTotalSize*=n[i];
   } 
   sizeout = int(double(fTotalSize)*(int)(n[ndim-1]/2+1)/n[ndim-1]);
   if (!inPlace) {
      //fIn =malloc(sizeof(double)*fTotalSize);
      try
      {
          fOut = new double [2*sizeout];
       }
      catch (std::bad_alloc&)
      {
        std::cout << "Error allocating memory." << std::endl;
        exit(1);
      }    
   }
   if (already_reserved == NULL){
      try
      {
          fIn  = new double [2*sizeout];
      }
      catch (std::bad_alloc&)
      {
        std::cout << "Error allocating memory." << std::endl;
        exit(1);
      }    
   } 
   else {
      destroy_fIn=false;
      fIn  = already_reserved;
   }
}

//_____________________________________________________________________________
xmippFftw::xmippFftw(Matrix1D<double> &img, bool my_inPlace, 
                                            bool already_reserved)
{
    int * myfN; 
    int ndim    = 1;
    myfN    = new int [ndim];     
    myfN[0] = XSIZE(img);

    if(already_reserved)
       myxmippFftw(myfN[0], my_inPlace,MULTIDIM_ARRAY(img));
    else
       myxmippFftw(myfN[0], my_inPlace,NULL);
    
}
//_____________________________________________________________________________
/*

If you have an array stored in column-major order and wish to transform it using FFTW,
it is quite easy to do. When creating the plan, simply pass the dimensions of the
array to the planner in reverse order. For example, if your array is a rank three N x
M x L matrix in column-major order, you should pass the dimensions of the array as if
it were an L x M x N matrix (which it is, from the perspective of FFTW).

Be Aware that you will get the transpose Fourier transform,
*/
xmippFftw::xmippFftw(Matrix2D<double> &img, bool init_and_do_transform)
{
    int * myfN; 
    int ndim    = 2; 
    myfN    = new int [ndim];     
    /* map y to x dim */
    myfN[0] = YSIZE(img);
    myfN[1] = XSIZE(img);
    //#define xmippFftw_Matrix2D
    #ifdef  xmippFftw_Matrix2D
    std::cerr << "XSIZE(img) YSIZE(img) " << XSIZE(img) << " " 
                                          << YSIZE(img) << std::endl;
    #endif
    #undef xmippFftw_Matrix2D
    bool my_inPlace = false;//no implace transform for 2D
    myxmippFftw(ndim, myfN, my_inPlace,MULTIDIM_ARRAY(img));
    // I presume if you call this routine with a matrix you 
    // want to do a forward fourier transform
    if(init_and_do_transform)
        {
        Init("ES",FFTW_FORWARD,false);
        Transform();
        }
}

/*

If you have an array stored in column-major order and wish to transform it using FFTW,
it is quite easy to do. When creating the plan, simply pass the dimensions of the
array to the planner in reverse order. For example, if your array is a rank three N x
M x L matrix in column-major order, you should pass the dimensions of the array as if
it were an L x M x N matrix (which it is, from the perspective of FFTW).

*/
xmippFftw::xmippFftw(Matrix3D<double> &img, bool init_and_do_transform)
{
    int * myfN; 
    int ndim    = 3; 
    myfN    = new int [ndim];     
    //map x to z dim
    myfN[0] = ZSIZE(img);
    myfN[1] = YSIZE(img);
    myfN[2] = XSIZE(img);

    bool my_inPlace = false;//no implace transform for 2D
    myxmippFftw(ndim, myfN, my_inPlace,MULTIDIM_ARRAY(img));
    // I presume if you call this routine with a matrix you 
    // want to do a forward fourier transform
    if(init_and_do_transform)
        {
        Init("ES",FFTW_FORWARD,false);
        Transform();
        }
}

//_____________________________________________________________________________
xmippFftw::~xmippFftw()
{
//Destroys the data arrays and the plan. However, some plan information stays around
//until the root session is over, and is reused if other plans of the same size are
//created

   fftw_destroy_plan((fftw_plan)fPlan);
   fPlan = NULL;
   if(fIn!=NULL && destroy_fIn==true)
   {
       delete [] fIn;
       fIn = NULL;
   }
   else
       fIn = NULL;
   if (fOut!=NULL)
   {
       delete [] fOut;
   }
   else
       fOut = NULL;

   if (fN!=NULL)
   {
      delete[] fN;
      fN = NULL;
   }
   else
      fN = NULL;
}

//_____________________________________________________________________________
//sometimes memory is importan and you want to free fIn while keeping fOut
void xmippFftw::delete_fIn()
{
   if(fIn!=NULL)
   {
       delete [] fIn;
       fIn = 0;
   }
}


//_____________________________________________________________________________
void xmippFftw::Init(std::string flags,int sign,bool wisdom_flag)
{
//Creates the fftw-plan
//
//NOTE:  input and output arrays are overwritten during initialisation,
//       so don't set any points, before running this function!!!!!
//
//Possible flag_options:
//"ES" (from "estimate") - no time in preparing the transform, but probably sub-optimal
//   performanc
//"M" (from "measure") - some time spend in finding the optimal way to do the transform
//"P" (from "patient") - more time spend in finding the optimal way to do the transform
//"EX" (from "exhaustive") - the most optimal way is found
//This option should be chosen depending on how many transforms of the same size and
//type are going to be done. Planning is only done once, for the first transform of this
//size and type.
// After trying wisdom I think is totally useless ROBERTO


//     void fftw_export_wisdom_to_file(FILE *output_file);
//     int fftw_import_wisdom_from_file(FILE *input_file);

//save wishdon using a lock mechanism and 
//delete when xmipp recompile -> ¿/tmp/fftw?
    /* Get any accumulated wisdom. */
/*
std::cerr << "Init_-3" << "plan= " << fPlan << std::endl;   
    if (fPlan != NULL)
    {
          std::cerr << "Init_-2" << "plan= " << fPlan << std::endl;   
         fftw_destroy_plan((fftw_plan)fPlan);
    }
*/
    FILE *wisdom;
    if(wisdom_flag)
    {
        wisdom = fopen((wisdom_name+".data").c_str(), "r");
		if (wisdom) {
			fftw_import_wisdom_from_file(wisdom);
			fclose(wisdom);
		}
    }    
    fFlags = flags;
    fSign  = sign;
    if(fSign == FFTW_FORWARD)
    {
       if (fOut==NULL)
       {
          fPlan = (void*)fftw_plan_dft_r2c(fNdim, fN, (double*)fIn, (fftw_complex*)fIn,MapFlag(flags));
       }
       else
       {
          fPlan = (void*)fftw_plan_dft_r2c(fNdim, fN, (double*)fIn, (fftw_complex*)fOut, MapFlag(flags));
       }
    }
    else if (fSign == FFTW_BACKWARD)
       {   
       if (fOut==NULL)
          fPlan = (void*)fftw_plan_dft_c2r(fNdim, fN,(fftw_complex*)fIn,(double*)fIn, MapFlag(flags));
       else
          fPlan = (void*)fftw_plan_dft_c2r(fNdim, fN, (fftw_complex*)fIn, (double*)fOut, MapFlag(flags));
       }
    else
       {
       std::cerr << "Invalid sign value:" 
                 << fSign 
                 << " in xmippFftw::Init "
                 <<std::endl;
       }   
    /* Save the wisdom. */
    if(wisdom_flag)
     {
		wisdom = fopen((wisdom_name+".data").c_str(), "w");
		if (wisdom) {
			fftw_export_wisdom_to_file(wisdom);
			fclose(wisdom);
		}
     }
    // Check plan was created 
    if(fPlan==NULL){
       std::cerr << "Error in xmippFftw::Init. Plan cannot be created\n"; 
    }    
}

//_____________________________________________________________________________
void xmippFftw::Transform()
{
//Computes the transform, specified in Init() function
#define MYSIZE 8

   if (fPlan){
      fftw_execute((fftw_plan)fPlan);
   }
   else {
      std::cerr << "xmippFftw:Transform, transform hasn't been initialised"
                << std::endl;
      exit(1);
   }
}

//_____________________________________________________________________________
void xmippFftw::GetPoints(double * data, bool fromInput) const
{
//Fills the array data with the computed transform.
//or the initial points
//Only (roughly) a half of the transform is copied (exactly the output of FFTW),
//the rest being Hermitian symmetric with the first half
// from input mens get original input data
   if(fSign == FFTW_FORWARD)
   {
       if (fromInput){
          for (int i=0; i<fTotalSize; i++){
             data[i] = fIn[i];
          }
       } else {
          int realN = 2*int(double(fTotalSize)*(fN[fNdim-1]/2+1)/fN[fNdim-1]);
          if (fOut){
             for (int i=0; i<realN; i+=2){
                data[i]    = ((fftw_complex*)fOut)[i/2][0];
                data[i+1]  = ((fftw_complex*)fOut)[i/2][1];
             }
          }
          else {
             for (int i=0; i<realN; i++)
                data[i] = (double) fIn[i];
          }
       }
    }   
    else if (fSign == FFTW_BACKWARD)
    {
        if (fromInput){
             int realN = 2*int(double(fTotalSize)*(fN[fNdim-1]/2+1)/fN[fNdim-1]);
             for (int i=0; i<realN; i+=2){
                data[i]    = ((fftw_complex*)fIn)[i/2][0];
                data[i+1]  = ((fftw_complex*)fIn)[i/2][1];
             }
        }
        if (fOut){
           for (int i=0; i<fTotalSize; i++)
              data[i] = ((double*)fOut)[i];
        }
        else{
           for (int i=0; i<fTotalSize; i++)
              data[i] = ((double*)fIn)[i];
        }
    }   
    else
    {
       std::cerr << "Invalid sign value:" 
                 << fSign 
                 << " in xmippFftw::GetPoints "
                 <<std::endl;
    }

}

//_____________________________________________________________________________
void xmippFftw::GetPoints(std::complex<double> * data, bool fromInput) const
{
//Fills the array data with the computed transform.
//or the initial points
//Only (roughly) a half of the transform is copied (exactly the output of FFTW),
//the rest being Hermitian symmetric with the first half
// from input mens get original input data
   if(fSign == FFTW_FORWARD)
   {
       if (fromInput){
           REPORT_ERROR(1,"FFTW: Not implemented 1");
       } else {
          int realN = 2*int(double(fTotalSize)*(fN[fNdim-1]/2+1)/fN[fNdim-1]);
          if (fOut){
             for (int i=0; i<realN; i+=2){
                 ((fftw_complex*)data)[i/2][0]=fOut[i];
                 ((fftw_complex*)data)[i/2][1]=fOut[i+1];

             }
          }
          else {
             REPORT_ERROR(1,"FFTW: Not implemented 2");
          }
       }
    }   
    else if (fSign == FFTW_BACKWARD)
    {
        REPORT_ERROR(1,"FFTW: Not implemented 3");
    }   
    else
    {
       std::cerr << "Invalid sign value:" 
                 << fSign 
                 << " in xmippFftw::GetPoints "
                 <<std::endl;
    }

}

//_____________________________________________________________________________
void xmippFftw::SetPoints(const double *data)
{
//Set all input points

   if(fSign == FFTW_FORWARD)
   {
       for (int i=0; i<fTotalSize; i++){
          ((double*)fIn)[i]=data[i];
       }
   }    
//set all points. the values are copied. points should be ordered as follows:
//[re_0, im_0, re_1, im_1, ..., re_n, im_n)
    else if (fSign == FFTW_BACKWARD)
    {

       int sizein = int(double(fTotalSize)*(fN[fNdim-1]/2+1)/fN[fNdim-1]);

       for (int i=0; i<2*(sizein); i+=2){
          ((fftw_complex*)fIn)[i/2][0]=data[i];
          ((fftw_complex*)fIn)[i/2][1]=data[i+1];
       }
    }       
}

//_____________________________________________________________________________
void xmippFftw::SetPoints(const std::complex<double> *data)
{
//Set all input points

   if(fSign == FFTW_FORWARD)
   {
       REPORT_ERROR(1,"FFTW: Not implemented 4");
   }    
//set all points. the values are copied. points should be ordered as follows:
//[re_0, im_0, re_1, im_1, ..., re_n, im_n)
    else if (fSign == FFTW_BACKWARD)
    {

       int sizein = int(double(fTotalSize)*(fN[fNdim-1]/2+1)/fN[fNdim-1]);
       for (int i=0; i<(sizein); i++){
           ((std::complex<double>*)fIn)[i]=data[i];
       }
    }       
}

//_____________________________________________________________________________
void xmippFftw::Normalize(void)
{
   if (inPlace)   
   {
      for (int i=0; i<fTotalSize; i++)
         ((double*)fIn)[i]  /= fTotalSize;
   }
   else   
   {
       for (int i=0; i<fTotalSize; i++)
         ((double*)fOut)[i]   /= fTotalSize;
   }
}

//_____________________________________________________________________________
unsigned xmippFftw::MapFlag(std::string opt)
{
//allowed options:
//"ES"
//"M"
//"P"
//"EX"
   unsigned myflag;
   myflag=FFTW_ESTIMATE;
   if (opt =="ES")
      myflag=FFTW_ESTIMATE;
   else if (opt == "M")
      myflag=FFTW_MEASURE;
   else if (opt == "P")
      myflag=FFTW_PATIENT;
   else if (opt == "EX")
      myflag=FFTW_EXHAUSTIVE;
   else
      std::cerr << "No valid Plan opt using FFTW_ESTIMATE (MapFlag)";   
   return myflag; 
}


void xmippFftw::TimeLimit(double seconds)
{
/*

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

fftw_set_timelimit( seconds);
}

void xmippFftw::CenterRealDataBeforeTransform(void)
{
  int ii=0;
  for(int i=0; i<fN[1]; i++)
     for(int j=0; j<fN[0]; j++)
       {
       if(((i+j)%2) == 1) fIn[ii] *= (-1.);
       ii++;
       }
}

/* Applies a bandpass filter to an image. */
/* frecuencies in range 0-0.5 */
void xmippFftw::img_bandpass_filter(double res_hi, double width)
{
    int Xdim,Xsize;
    int Ydim=1,Ysize=1;
    int Zdim=1,Zsize=1;
    /* rememebr to transpose de matrix */
    if(fNdim==1)
    {   
        Zdim=Zsize=1;
        Ydim=Ysize=1;
        Xsize=fN[0];
        Xdim = (int) fN[0]/2 +1;        
    }
    else if (fNdim==2)
    {
        Zdim=Zsize=1;
        Ydim=Ysize=fN[0];
        Xsize=fN[1];
        Xdim = (int) fN[1]/2 +1;
    
    }    
    else if (fNdim==3)
    {
        Zdim=Zsize=fN[0];
        Ydim=Ysize=fN[1];
        Xsize=fN[2];
        Xdim = (int)  fN[2]/2 +1;
    }
    else
    {
        std::cerr << "Error in img_bandpass_filter\n";
        exit(0);
    }
    int zz, yy, xx;
    double sz2,sy2,sx2,s;   

	double res_hi2= res_hi*res_hi;
    std::complex<double> * cfOut, *cfIn;
    cfOut = (std::complex<double> *)fOut; //fOut is double *
    cfIn  = (std::complex<double> *)fIn;  //fIn is double *
    for ( int z=0, i=0; z<Zdim; z++ ) {
		if ( z > (Zsize - 1)/2 ) zz = Zsize-z;
		else zz = z;
        sz2 = (double)zz/Zsize;
		sz2 *= sz2;
		for ( int y=0; y<Ydim; y++ ) {
			if ( y > (Ysize - 1)/2 ) yy = Ysize-y;
			else yy = y;
			sy2 = (double)yy/Ysize;
			sy2 *= sy2;
			for ( int x=0; x<Xdim; x++, i++ ) {
				if ( x > (Xsize - 1)/2 ) xx = Xsize-x;
				else xx = x;
				sx2 = (double)xx/Xsize;
				sx2 *= sx2;
				s = (sx2 + sy2 + sz2);
                if (s > res_hi2) 
				    cfOut[i] = 0.;
/*                else    
				    cfOut[i] = 1.;

std::cerr << "i xx yy s res_hi2 z y x " 
          << i << " " << sy2 << " " << sx2 
          << " " << s << " " 
          << res_hi2 << " "
          << z << " "
          << y << " "
          << x << " "
          << fOut[2*i] << "\n";
*/
			}
		}
	}

}

/** Fourier-Ring-Correlation between two 2D-matrices using FFT
 * @ingroup FourierOperations
 */
void xmippFftw::fourier_ring_correlation(xmippFftw & fft_m2,
                              double sampling_rate,
                              Matrix1D< double >& freq,
                              Matrix1D< double >& frc,
                              Matrix1D< double >& frc_noise)
{
    int dim;
    /* remember to transpose de matrix */
    if(fNdim==1)
    {
        dim = (int)  fN[0]/2;
    }
    else if (fNdim==2)
    {
        dim = (int)  fN[1]/2;
    }
    else if (fNdim==3)
    {
        dim = (int)  fN[2]/2;
    }
    else
    {
        std::cerr << "Error in fftwRadialAverage\n";
        exit(0);
    }

    double     *aux, *realFT1;
    int sizeout = int(double(fTotalSize)*(int)(fN[fNdim-1]/2+1)/fN[fNdim-1]);
    try
    {
       aux = new double [2*sizeout];
       realFT1 = new double [2*sizeout];
    }
    catch (std::bad_alloc&)
    {
        std::cout << "Error allocating memory." << std::endl;
        exit(1);
    }

    Matrix1D< int >    radial_count;
    Matrix1D< double > tmp1, tmp2;
    //xmipp y contigous, fftw x contigous
    std::complex<double> * IMG1, * IMG2;//,*AUX;
    IMG1 = (std::complex<double> *)fOut;
    //AUX = (std::complex<double> *)aux;
    for (int ii=0;ii<sizeout;ii++)
    {
       aux[ii] = abs(IMG1[ii])*abs(IMG1[ii]);
    }
    tmp1.initZeros();
    fftwRadialAverage(aux, tmp1, radial_count,true);
    IMG2 = (std::complex<double> *)fft_m2.fOut;
    for (int ii=0;ii<sizeout;ii++)
    {
       aux[ii] = abs(IMG2[ii])*abs(IMG2[ii]);
    }
    tmp2.initZeros();
    fftwRadialAverage(aux, tmp2, radial_count, true);
    //radialAverageHalf

    //std::complex<double> * REALFT1;
    //REALFT1 = (std::complex<double> *)realFT1;
    for (int ii=0;ii<sizeout;ii++)
    {
       realFT1[ii] = real(conj(IMG1[ii]) * IMG2[ii]);
    }
    
    frc.initZeros();
    fftwRadialAverage(realFT1, frc, radial_count, true);
    frc.resize(dim);
    frc_noise.resize(dim);
    freq.resize(dim);
    
    FOR_ALL_ELEMENTS_IN_MATRIX1D(freq)
    {
        int j = i;
        VEC_ELEM(freq, i) = (double) j / (dim * 2 * sampling_rate);
        
        VEC_ELEM(frc, i) = VEC_ELEM(frc, i) / sqrt(VEC_ELEM(tmp1, i)
                           * VEC_ELEM(tmp2, i));
        VEC_ELEM(frc_noise, i) = 2 / sqrt((double) VEC_ELEM(radial_count, i));
    }
}


/* Radial average for Fourier transforms*/
void xmippFftw::fftwRadialAverage(double * AUX,
                                  Matrix1D< double >& radial_mean,
                                  Matrix1D< int >& radial_count,
                                  bool rounding /*=true*/ )
{
    int Xdim,Xsize;
    int Ydim=1,Ysize=1;
    int Zdim=1,Zsize=1;
    int maxDistance;
    /* remember to transpose de matrix */
    if(fNdim==1)
    {   
        Zdim=Zsize=1;
        Ydim=Ysize=1;
        Xsize=fN[0];
        Xdim = (int) fN[0]/2 +1;
        maxDistance=Xdim;        
    }
    else if (fNdim==2)
    {
        Zdim=Zsize=1;
        Ydim=Ysize=fN[0];
        Xsize=fN[1];
        Xdim = (int) fN[1]/2 +1;
        maxDistance = (int) FLOOR(sqrt(Xdim*Xdim+(Ydim/2)*(Ydim/2)));
    }    
    else if (fNdim==3)
    {
        Zdim=Zsize=fN[0];
        Ydim=Ysize=fN[1];
        Xsize=fN[2];
        Xdim = (int)  fN[2]/2 +1;
        maxDistance = (int) FLOOR(sqrt(Xdim*Xdim+
                                  (Ydim/2)*(Ydim/2)+
                                  (Zdim/2)*(Zdim/2)));
    }
    else
    {
        std::cerr << "Error in fftwRadialAverage\n";
        exit(0);
    }
    // Define the vectors
    radial_mean.resize(maxDistance+1);//I think +1 will be needed for rounding
                                      //case. ROB
    radial_mean.initZeros();
    radial_count.resize(maxDistance+1);
    radial_count.initZeros();

    int zz, yy, xx;
    double sz2,sy2,sx2,s;   
    /* maximun distance*/
    std::complex<double> * cfOut;
    cfOut = (std::complex<double> *)fOut; //fOut is double *
    for ( int z=0, ii=0; z<Zdim; z++ ) {
		if ( z > (Zsize - 1)/2 ) zz = Zsize-z;
		else zz = z;
		for ( int y=0; y<Ydim; y++ ) {
			if ( y > (Ysize - 1)/2 ) yy = Ysize-y;
			else yy = y;
			for ( int x=0; x<Xdim; x++, ii++ ) {
				if ( x > (Xsize - 1)/2 ) xx = Xsize-x;
				else xx = x;
				s = sqrt(xx*xx+yy*yy+zz*zz);
                // Determine distance to the center
                int distance;
                if (rounding)
                    distance = (int) ROUND(s);
                else
                    distance = (int) FLOOR(s);
                // Sum the value to the pixels with the same distance
                radial_mean(distance) += AUX[ii];

                // Count the pixel
                radial_count(distance)++;
                //each point counts by two
                //except those in the redundant part
                if (xx!=0)
                {
                    radial_mean(distance) += AUX[ii];
                    radial_count(distance)++;
                }
			}
		}
	}
    // Perfor the mean
    FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
    {
        radial_mean(i) /= (double) radial_count(i);
    }

}

//
// lockFile is the filename to use to represent the lock.  Basically, if the
// file exists, then someone has the lock (in this case to prtsuffix.data
// If the lockFile does not exist, then there is not currently a lock
//


Lock::Lock (std::string basefilename)
{
    // ------------------------------
    // Attempt to create a new lockFile...
    // ------------------------------
    int   count = 0;
    time_t delta, current_time;
    struct stat statptr;
    filename = basefilename + "LOCK";

   /*
    * Get the time of last modification of the /etc/boottime file, which is
    * created each time the system boots.
    */

    non_locked=false;
    std::ofstream      l;
    if  (stat(filename.c_str(),&statptr) )
    {
        time(&current_time); /* get the current time */
        ctime(&statptr.st_mtime);
        delta=current_time-statptr.st_atime; /* calculate the delta time in seconds*/
        if (delta > 3600)
            unlink (filename.c_str());
    }
    else //if file does not exists open it
       l.open (filename.c_str(), std::ios::out);
    //
    // If we were unable to create the file, then sleep and try again.
    // In this particular example, I know that the reason I am locking the
    // file is pretty quick, so we should be able to get a lock in a
    // reasonable amount of time.  If 11 seconds go by and we still don't
    // have a lock, then another program probably died without freeing the
    // lock, so just assume responsibility for the lock and go on as if it is ours.
    //
    while (!l && (++count < 11))
   {
       sleep (1);
       if( !stat(filename.c_str(),&statptr))
          l.open (filename.c_str(), std::ios::out);
   }  
   if (l)
   {
       l << "locked.\n";
       l.close();
       chmod(filename.c_str(),0666);
       non_locked=true;
   }

}

Lock::~Lock (void)
{
 // ------------------------------
 // Remove the lock file.
 // ------------------------------
 if(non_locked)
   unlink (filename.c_str());

}
/** Example
int main (void)
{
 Lock lock("/tmp/fftw"); // we now have an exclusive lock
 if(lock.non_locked)
     use the file, is yours...
 // do something

} // lock gets reliquished
*/
