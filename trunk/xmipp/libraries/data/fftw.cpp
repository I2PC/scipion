
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
   std::cerr << "check_size INT_MAX " << check_size << " " << INT_MAX << std::endl;
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
   sizeout = int(double(fTotalSize)*(n[ndim-1]/2+1)/n[ndim-1]);
   if (!inPlace){
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
    int ndim    = 1; //useless for one dimension
                 //the funtion getdimension resturns size for 1D!!
    myfN    = new int [ndim];     
    myfN[0] = XSIZE(img);

    if(already_reserved)
       myxmippFftw(myfN[0], my_inPlace,MULTIDIM_ARRAY(img));
    else
       myxmippFftw(myfN[0], my_inPlace,NULL);
    
}
//_____________________________________________________________________________
xmippFftw::xmippFftw(Matrix2D<double> &img, bool my_inPlace, 
                                            bool already_reserved)
{
    int * myfN; 
    int ndim    = 2; //useless for one dimension
                 //the funtion getdimension resturns size for 1D!!
    myfN    = new int [ndim];     
    myfN[0] = XSIZE(img);
    myfN[1] = YSIZE(img);

    if(already_reserved)
        myxmippFftw(ndim, myfN, my_inPlace,MULTIDIM_ARRAY(img));
    else
        myxmippFftw(ndim, myfN, my_inPlace,NULL);
}

xmippFftw::xmippFftw(Matrix3D<double> &img, bool my_inPlace, 
                                            bool already_reserved)
{
    int * myfN; 
    int ndim    = 3; //useless for one dimension
                 //the funtion getdimension resturns size for 1D!!
    myfN    = new int [ndim];     
    myfN[0] = XSIZE(img);
    myfN[1] = YSIZE(img);
    myfN[2] = ZSIZE(img);

    if(already_reserved)
        myxmippFftw(ndim, myfN, my_inPlace,MULTIDIM_ARRAY(img));
    else
        myxmippFftw(ndim, myfN, my_inPlace,NULL);
}

//_____________________________________________________________________________
xmippFftw::~xmippFftw()
{
//Destroys the data arrays and the plan. However, some plan information stays around
//until the root session is over, and is reused if other plans of the same size are
//created

   fftw_destroy_plan((fftw_plan)fPlan);
   fPlan = 0;
   if(destroy_fIn)
   {
       delete [] fIn;
       fIn = 0;
   }
   else
       fIn = 0;
   if (!inPlace)
   {
       delete [] fOut;
   }
   else
       fOut = 0;

   if (fN!=NULL)
   {
      delete[] fN;
      fN = 0;
   }
   else
      fN = 0;

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
          std::cerr << "Creating plan with: \n" 
                    << "fNdim " <<  fNdim
                    << " fN    " <<  fN[0] << " " << fN[1]    
                    << " flags " <<  flags
                    << "\n" ;
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
    long delta, current_time;
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
