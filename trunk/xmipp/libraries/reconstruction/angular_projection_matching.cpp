/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.uam.es (2004)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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

#include "angular_projection_matching.h"

//#define DEBUG
//#define TIMING

// For blocking of threads
pthread_mutex_t update_refs_in_memory_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t debug_mutex = PTHREAD_MUTEX_INITIALIZER;

// Read arguments ==========================================================
void Prog_angular_projection_matching_prm::read(int argc, char **argv)  {

  // Read command line
  if (checkParameter(argc,argv,"-more_options")) { usage(); extendedUsage();}
  fn_ref  = getParameter(argc,argv,"-ref","ref");
  fn_exp  = getParameter(argc,argv,"-i");
  fn_root = getParameter(argc,argv,"-o","out");

  // Additional commands
  pad=XMIPP_MAX(1.,textToFloat(getParameter(argc,argv,"-pad","1.")));
  Ri=textToInteger(getParameter(argc,argv,"-Ri","-1"));
  Ro=textToInteger(getParameter(argc,argv,"-Ro","-1"));
  search5d_shift  = textToInteger(getParameter(argc,argv,"-search5d_shift","0"));
  search5d_step = textToInteger(getParameter(argc,argv,"-search5d_step","2"));
  max_shift = textToFloat(getParameter(argc,argv,"-max_shift","-1"));
  avail_memory = textToFloat(getParameter(argc,argv,"-mem","1"));
  fn_ctf  = getParameter(argc,argv,"-ctf","");
  phase_flipped = checkParameter(argc, argv, "-phase_flipped");
  // Number of threads
  threads = textToInteger(getParameter(argc, argv, "-thr","1"));

  // Hidden stuff
  verb=textToInteger(getParameter(argc,argv,"-verb","1"));

}

// Show ====================================================================
void Prog_angular_projection_matching_prm::show() {

  if (verb>0) {
    std::cerr << "  Input images            : "<< fn_exp << std::endl;
    std::cerr << "  Output rootname         : "<< fn_root << std::endl;
    if (Ri>0)
	std::cerr << "  Inner radius rot-search : "<<Ri<<std::endl;
    if (Ro>0)
	std::cerr << "  Outer radius rot-search : "<<Ro<<std::endl;
    if (max_nr_refs_in_memory<total_nr_refs)
    {
	std::cerr << "  Number of references    : "<<total_nr_refs<<std::endl;
	std::cerr << "  Nr. refs in memory      : "<<max_nr_refs_in_memory<< " (using "
		  <<avail_memory<<" Gb)"<<std::endl;
    }
    else
    {
	std::cerr << "  Number of references    : "<<total_nr_refs<<" (all stored in memory)"<<std::endl;	
    }
    std::cerr << "  Max. allowed shift      : +/- " <<max_shift<<" pixels"<<std::endl;
    if (search5d_shift > 0)
    {
	std::cerr << "  5D-search shift range   : "<<search5d_shift<<" pixels (sampled "<<nr_trans<<" times)"<<std::endl;
    }
    if (fn_ctf!="")
    {
        if (Is_ImageXmipp(fn_ctf))
        {
            std::cerr << "  CTF image               :  " <<fn_ctf<<std::endl;
            if (pad > 1.)
                std::cerr << "  Padding factor          : "<< pad << std::endl;
        }
        else
        {
            std::cerr << "  CTF parameter file      :  " <<fn_ctf<<std::endl;
            if (phase_flipped)
                std::cerr << "    + Assuming images have been phase flipped " << std::endl;
            else
                std::cerr << "    + Assuming images have not been phase flipped " << std::endl;
        }
    }
    if (threads>1)
    {
        std::cerr << "  -> Using "<<threads<<" parallel threads"<<std::endl;
    }
    std::cerr << " ================================================================="<<std::endl;
  }
}

// Usage ===================================================================
void Prog_angular_projection_matching_prm::usage() {
  std::cerr << "Usage:  projection_matching [options] "<<std::endl;
  std::cerr << "   -i <docfile>            : Docfile with input images \n"
	    << "   -ref <ref rootname=\"ref\"> : Rootname for reference projection files \n"
            << "   -o [rootname=\"out\"]       : Output rootname \n"
            << " [ -search5d_shift <int=0> ]   : Search range (in +/- pix) for 5D shift search\n"
            << " [ -search5d_step  <int=2> ]   : Step size for 5D shift search (in pix) \n"
	    << " [ -Ri <float=1> ]             : Inner radius to limit rotational search \n"
	    << " [ -Ro <float=dim/2 - 1> ]     : Outer radius to limit rotational search \n"
	    << " [ -more_options ]             : Show all program options\n";
}

// Extended usage ===================================================================
void Prog_angular_projection_matching_prm::extendedUsage() {
  std::cerr << "Additional options: \n"
            << " [ -mem <float=1> ]            : Available memory for reference library (Gb)\n"
	    << " [ -max_shift <float=-1> ]     : Max. change in origin offset (+/- pixels; neg= no limit) \n"
	    << " [ -ctf <filename> ]           : CTF to apply to the reference projections, either a  \n"
            << "                                 CTF parameter file or a 2D image with the CTF amplitudes\n"
            << " [-pad <float=1>]              : Padding factor (for CTF correction only)\n"
	    << " [ -phase_flipped ]            : Use this if the experimental images have been phase flipped\n";
  exit(1);
}

// Side info stuff ===================================================================
void Prog_angular_projection_matching_prm::produceSideInfo() {

    ImageXmipp       img,empty;
    Projection       proj;
    DocFile          DF;
    SelFile          SFr,emptySF;
    SymList          SL;
    FileName         fn_img;
    double           mean,stddev,psi=0.;
    Matrix2D<double> Maux;
    Matrix1D<double> dataline(3);
    int              nl;
    Polar<double>    P;
    Polar<std::complex <double> > fP;

    // Read Selfile and get dimensions
    DFexp.read(fn_exp);

    // Thread barrier
    barrier_init(&thread_barrier, threads);

    // Read one image to get dim
    DFexp.go_first_data_line();
    DFexp.previous();
    if (DFexp.get_current_line().Is_comment()) 
	{
	fn_img = ((DFexp.get_current_line()).get_text()).erase(0, 3);
	}
    else
	REPORT_ERROR(1,"BUG: no comment in DFexp where expected....");
    img.read(fn_img);
    
    dim = XSIZE(img());
    
    // Set padding dimension
    paddim=ROUND(pad*dim);

    // Set max_shift
    if (max_shift<0) max_shift = dim/2;

    // Set ring defaults
    if (Ri<1) Ri=1;
    if (Ro<0) Ro=(dim/2)-1;

    // Calculate necessary memory per image
    img().produceSplineCoefficients(Maux,3);
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro);
    P.calculateFftwPlans(global_plans);
    fourierTransformRings(P,fP,global_plans,false);
    double memory_per_ref = 0.;
    for (int i = 0; i < fP.getRingNo(); i++)
    {
	memory_per_ref += (double) fP.getSampleNo(i) * 2 * sizeof(double);
    }
    memory_per_ref += dim * dim * sizeof(double);
    max_nr_imgs_in_memory = ROUND( 1024 * 1024 * 1024 * avail_memory / memory_per_ref);

    // Set up angular sampling
    mysampling.read_sampling_file(fn_ref,false);
    total_nr_refs = mysampling.no_redundant_sampling_points_angles.size();

    // Don't reserve more memory than necessary
    max_nr_refs_in_memory = XMIPP_MIN(max_nr_imgs_in_memory, total_nr_refs);

    // Initialize pointers for reference retrieval 
    pointer_allrefs2refsinmem.resize(total_nr_refs,-1);
    pointer_refsinmem2allrefs.resize(max_nr_refs_in_memory,-1);
    counter_refs_in_memory = 0;
    loop_forward_refs=true;

    // Initialize 5D search vectors
    search5d_xoff.clear();
    search5d_yoff.clear();
    // Make sure origin is included
    int myfinal=search5d_shift + search5d_shift%search5d_step;
    nr_trans = 0;
    for (int xoff = -myfinal; xoff <= myfinal; xoff+= search5d_step)
    {
	for (int yoff = -myfinal; yoff <= myfinal; yoff+= search5d_step)
	{
	    // Only take a circle (not a square)
	    if ( xoff*xoff + yoff*yoff <= search5d_shift*search5d_shift)
	    {
		search5d_xoff.push_back(xoff);
		search5d_yoff.push_back(yoff);
                nr_trans++;
	    }
	}
    }

    // Initialize all arrays
    try
    {
        fP_ref = new Polar<std::complex<double> >[max_nr_refs_in_memory];
        proj_ref = new Matrix2D<double>[max_nr_refs_in_memory];
        fP_img = new Polar<std::complex<double> >[nr_trans];
        fPm_img = new Polar<std::complex<double> >[nr_trans];
        
        stddev_ref = new double[max_nr_refs_in_memory];
        stddev_img = new double[nr_trans];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in produceSideInfo");
    }
  
    // CTF stuff
    if (fn_ctf != "")
    {
        if (Is_ImageXmipp(fn_ctf))
        {
            ImageXmipp img;
            img.read(fn_ctf);
	    Mctf=img();
            if (XSIZE(Mctf) != paddim)
            {
                std::cerr<<"image size= "<<dim<<" padding factor= "<<pad<<" padded image size= "<<paddim<<" Wiener filter size= "<<XSIZE(Mctf)<<std::endl;
                REPORT_ERROR(1,"Incompatible padding factor for this CTF filter");
            }
        }
        else
        {
            XmippCTF ctf;
            Matrix2D<std::complex<double> >  ctfmask;
            ctf.read(fn_ctf);
            if (ABS(ctf.DeltafV - ctf.DeltafU) >1.) 
            {
                REPORT_ERROR(1, "ERROR%% Only non-astigmatic CTFs are allowed!");
            }
            ctf.enable_CTF = true;
            ctf.Produce_Side_Info();
            ctf.Generate_CTF(paddim, paddim, ctfmask);
            Mctf.resize(paddim,paddim);
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mctf)
            {
                if (phase_flipped) dMij(Mctf, i, j) = fabs(dMij(ctfmask, i, j).real());
                else dMij(Mctf, i, j) = dMij(ctfmask, i, j).real();
            }
        }
    }

}

int Prog_angular_projection_matching_prm::getCurrentReference(int refno,
                                                              Polar_fftw_plans &local_plans)
{

    FileName                      fnt;
    ImageXmipp                    img;
    double                        mean,stddev;
    Matrix2D<double>              Maux;
    Polar<double>                 P;
    Polar<std::complex <double> > fP;
    XmippFftw                     local_transformer;

    // Image was not stored yet: read it from disc and store
    fnt.compose(fn_ref,refno+1,"xmp");
    img.read(fnt);
    img().setXmippOrigin();

    // Apply CTF
    if (fn_ctf!="")
    {
	Matrix2D<std::complex<double> > Faux;
        if (paddim > dim)
        {
            // pad real-space image
            int x0 = FIRST_XMIPP_INDEX(paddim);
            int xF = LAST_XMIPP_INDEX(paddim);
            img().window(x0, x0, xF,xF, 0.);
        }
        local_transformer.FourierTransform(img(),Faux);
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
	{
	    dMij(Faux,i,j) *= dMij(Mctf,i,j);
	}
	local_transformer.inverseFourierTransform(Faux,img());
        if (paddim > dim)
        {
            // de-pad real-space image
            int x0 = FIRST_XMIPP_INDEX(dim);
            int xF = LAST_XMIPP_INDEX(dim);
            img().window(x0, x0, xF,xF, 0.);
        }
    }

    // Calculate FTs of polar rings and its stddev
    img().produceSplineCoefficients(Maux,3);
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro);
    mean = P.computeSum(true);
    stddev = P.computeSum2(true);
    stddev = sqrt(stddev - mean * mean);
    P -= mean;
    fourierTransformRings(P,fP,local_plans,true);

    pthread_mutex_lock(  &update_refs_in_memory_mutex );

    int counter = counter_refs_in_memory % max_nr_refs_in_memory;
    pointer_allrefs2refsinmem[refno] = counter;
    if (pointer_refsinmem2allrefs[counter] != -1)
    {
	// This position was in use already
	// Images will be overwritten, so reset the
	// pointer_allrefs2refsinmem of the old images to -1
	pointer_allrefs2refsinmem[pointer_refsinmem2allrefs[counter]] = -1;
    }
    pointer_refsinmem2allrefs[counter] = refno;
    fP_ref[counter] = fP;
    stddev_ref[counter] = stddev;
    proj_ref[counter] = img();
 #ifdef DEBUG
    std::cerr<<"counter= "<<counter<<"refno= "<<refno<<" stddev = "<<stddev;
    std::cerr<<" refsinmem2allrefs= "<<pointer_refsinmem2allrefs[counter];
    std::cerr<<" allrefs2refsinmem= "<<pointer_allrefs2refsinmem[pointer_refsinmem2allrefs[counter]] <<std::endl;

#endif
   
    counter_refs_in_memory++;
    pthread_mutex_unlock(  &update_refs_in_memory_mutex );

}

void * threadRotationallyAlignOneImage( void * data )
{
    structThreadRotationallyAlignOneImage * thread_data = (structThreadRotationallyAlignOneImage *) data;

    // Variables from above
    int thread_id = thread_data->thread_id;
    int thread_num = thread_data->thread_num;
    Prog_angular_projection_matching_prm *prm = thread_data->prm;
    Matrix2D<double> *img = thread_data->img;
    int *this_image = thread_data->this_image;
    int *opt_refno = thread_data->opt_refno;
    double *opt_psi = thread_data->opt_psi;
    double *opt_flip = thread_data->opt_flip; 
    double *maxcorr = thread_data->maxcorr;

    // Local variables
    Matrix2D<double>            Maux;
    Matrix1D<double>            ang, corr;
    int                         max_index, refno, myinit, myfinal, myincr;
    bool                        done_once=false;
    double                      mean, stddev;
    Polar<double>               P;
    Polar<std::complex <double> > fP,fPm;
    XmippFftw                   local_transformer;
    Polar_fftw_plans            local_plans;
    int                         imgno = *this_image;

#ifdef TIMING
    TimeStamp t0,t1,t2; 
    time_config();
    annotate_time(&t0);
    annotate_time(&t2);
#endif

    *maxcorr = -99.e99;
    (*img).produceSplineCoefficients(Maux,3);
    // Precalculate polar transform of each translation
    // This loop is also threaded
    myinit = thread_id;
    myincr = thread_num;
    for (int itrans = myinit; itrans < prm->nr_trans; itrans+=myincr)
    {
        P.getPolarFromCartesianBSpline(Maux,prm->Ri,prm->Ro,
                                       (double)prm->search5d_xoff[itrans],
                                       (double)prm->search5d_yoff[itrans]);
	mean = P.computeSum(true);
	stddev = P.computeSum2(true);
	stddev = sqrt(stddev - mean * mean);
	P -= mean; // for normalized cross-correlation coefficient
        if (itrans == myinit)
            P.calculateFftwPlans(local_plans);
        fourierTransformRings(P,prm->fP_img[itrans],local_plans,false);
        fourierTransformRings(P,prm->fPm_img[itrans],local_plans,true);
	prm->stddev_img[itrans] = stddev;
        done_once=true;
    }
    // If thread did not have to do any itrans, initialize fftw plans
    if (!done_once)
    {
        P.getPolarFromCartesianBSpline(Maux,prm->Ri,prm->Ro);
        P.calculateFftwPlans(local_plans);
    }
    // Prepare FFTW plan for rotational correlation
    corr.resize(P.getSampleNoOuterRing());
    local_transformer.setReal(corr);
    local_transformer.FourierTransform();

    // All threads have to wait until the itrans loop is done
    barrier_wait(&(prm->thread_barrier));

#ifdef TIMING
    float prepare_img = elapsed_time(t0);
    float get_refs = 0.;
    annotate_time(&t0);
#endif

    //pthread_mutex_lock(  &debug_mutex );
    // Switch the order of looping through the references every time.
    // That way, in case max_nr_refs_in_memory<total_nr_refs
    // the references read in memory for the previous image
    // will still be there when processing the next image
    // Every thread processes a part of the references
    if (prm->loop_forward_refs)
    {
        myinit = 0;
        myfinal = prm->mysampling.my_neighbors[imgno].size();
        myincr = +1;
    }
    else
    {
	myinit = prm->mysampling.my_neighbors[imgno].size() - 1;
	myfinal = -1;
	myincr = -1;
    }

    // Loop over all relevant "neighbours" (i.e. directions within the search range)
    //for (int i = myinit; i != myfinal; i+=myincr)
    for (int i = myinit; i != myfinal; i += myincr)
    {
        if (i%thread_num == thread_id)
        {

#ifdef DEBUG_THREADS
            pthread_mutex_lock(  &debug_mutex );
            std::cerr<<" thread_id= "<<thread_id<<" i= "<<i<<" "<<myinit<<" "<<myfinal<<" "<<myincr<<std::endl;
            pthread_mutex_unlock(  &debug_mutex );
#endif

#ifdef TIMING
            annotate_time(&t1);
#endif
            // Get pointer to the current reference image
            refno = prm->pointer_allrefs2refsinmem[prm->mysampling.my_neighbors[imgno][i]];
            if (refno == -1)
            {
                // Reference is not stored in memory (anymore): (re-)read from disc
                prm->getCurrentReference(prm->mysampling.my_neighbors[imgno][i],local_plans);
                refno = prm->pointer_allrefs2refsinmem[prm->mysampling.my_neighbors[imgno][i]];
            }
#ifdef TIMING
            get_refs += elapsed_time(t1);
#endif
#ifdef DEBUG
            std::cerr<<"Got refno= "<<refno<<" pointer= "<<prm->mysampling.my_neighbors[imgno][i]<<std::endl;
#endif	    
        
            // Loop over all 5D-search translations
            for (int itrans = 0; itrans < prm->nr_trans; itrans++)
            {
                // A. Check straight image
                rotationalCorrelation(prm->fP_img[itrans],prm->fP_ref[refno],ang,local_transformer);
                corr /= prm->stddev_ref[refno] * prm->stddev_img[itrans]; // for normalized ccf
                for (int k = 0; k < XSIZE(corr); k++)
                {
                    if (corr(k)> *maxcorr)
                    {
                        *maxcorr = corr(k);
                        *opt_psi = ang(k);
                        *opt_refno = prm->mysampling.my_neighbors[imgno][i];
                        *opt_flip = 0.;
                    }
                }
#ifdef DEBUG
                std::cerr<<"straight: corr "<<*maxcorr<<std::endl;
#endif	    

                // B. Check mirrored image
                rotationalCorrelation(prm->fPm_img[itrans],prm->fP_ref[refno],ang,local_transformer);
                corr /= prm->stddev_ref[refno] * prm->stddev_img[itrans]; // for normalized ccf
                for (int k = 0; k < XSIZE(corr); k++)
                {
                    if (corr(k)> *maxcorr)
                    {
                        *maxcorr = corr(k);
                        *opt_psi = ang(k);
                        *opt_refno = prm->mysampling.my_neighbors[imgno][i];
                        *opt_flip = 1.;
                    }
                }
                
#ifdef DEBUG
                std::cerr<<"mirror: corr "<<*maxcorr;
                if (*opt_flip==1.) std::cerr<<"**";
                std::cerr<<std::endl;
#endif	    
            }
        }
    }
#ifdef TIMING
    float all_rot_align = elapsed_time(t0);
    float total_rot = elapsed_time(t2);
    std::cerr<<" rotal%% "<<total_rot
	     <<" => prep: "<<prepare_img
	     <<" all_refs: "<<all_rot_align
	     <<" (of which "<<get_refs
	     <<" to get "<< prm->mysampling.my_neighbors[imgno].size()
	     <<" refs for imgno "<<imgno<<" )"
	     <<std::endl;
#endif
    //pthread_mutex_unlock(  &debug_mutex );

}

void Prog_angular_projection_matching_prm::translationallyAlignOneImage(Matrix2D<double> &img,
								    const int &opt_refno, 
								    const double &opt_psi, 
								    const double &opt_flip, 
								    double &opt_xoff, 
								    double &opt_yoff, 
								    double &maxcorr)
{

    Matrix2D<double> Mtrans,Mimg,Mref;
    int refno;
    Mtrans.setXmippOrigin();
    Mimg.setXmippOrigin();
    Mref.setXmippOrigin();

#ifdef TIMING
    TimeStamp t0,t1,t2; 
    time_config();
    annotate_time(&t0);
#endif

#ifdef DEBUG
    std::cerr<<"start trans: opt_refno= "<<opt_refno<<" pointer= "<<pointer_allrefs2refsinmem[opt_refno]<<" opt_psi= "<<opt_psi<<"opt_flip= "<<opt_flip<<std::endl;
#endif	    

    // Get pointer to the correct reference image in memory
    refno = pointer_allrefs2refsinmem[opt_refno];
    if (refno == -1)
    {
	// Reference is not stored in memory (anymore): (re-)read from disc
	getCurrentReference(opt_refno,global_plans);
	refno = pointer_allrefs2refsinmem[opt_refno];
    }

    // Rotate stored reference projection by phi degrees
    proj_ref[refno].rotateBSpline(3,-opt_psi,Mref,DONT_WRAP);

#ifdef DEBUG
    std::cerr<<"rotated ref "<<std::endl;
#endif	    

    if (opt_flip > 0.)
    {
	// Flip experimental image
	Matrix2D<double> A(3,3);
	A.initIdentity();
	A(0, 0) *= -1.;
	A(0, 1) *= -1.;
	applyGeometry(Mimg,A, img, IS_INV, DONT_WRAP);
    }
    else
	Mimg = img;

    // Perform the actual search for the optimal shift
    if (max_shift>0) 
	best_shift(Mref,Mimg,opt_xoff,opt_yoff);
    else 
	opt_xoff = opt_yoff = 0.;
    if (opt_xoff * opt_xoff + opt_yoff * opt_yoff > max_shift * max_shift) 
	opt_xoff = opt_yoff = 0.;

#ifdef DEBUG
    std::cerr<<"optimal shift "<<opt_xoff<<" "<<opt_yoff<<std::endl;
#endif	 
   
    // Calculate standard cross-correlation coefficient
    Mimg.translate(vectorR2(opt_xoff,opt_yoff),Mtrans,true);
    maxcorr = correlation_index(Mref,Mtrans);

#ifdef DEBUG
    std::cerr<<"optimal shift corr "<<maxcorr<<std::endl;
#endif	 

    // Correct X-shift for mirrored images
    if (opt_flip>0.)
	opt_xoff *= -1.;	

#ifdef TIMING
    float total_trans = elapsed_time(t0);
    std::cerr<<" trans%% "<<total_trans <<std::endl;
#endif

}

void Prog_angular_projection_matching_prm::processSomeImages(int * my_images, double * my_output) 
{


  ImageXmipp img;
  double opt_rot, opt_tilt, opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr=-99.e99;
  int opt_refno;

  int nr_images = my_images[0];
  // Prepare my_output
  my_output[0] = nr_images * MY_OUPUT_SIZE;

  // Loop over all images
  for (int imgno = 0; imgno < nr_images; imgno++)
  {

      // add one because first number is number of elements in the array
      int this_image = my_images[imgno + 1]; 

      // read experimental image in memory
      //std::cerr << "get_image(this_image,img,true " <<  this_image << std::endl;
      //DFexp.get_image IS NOT EQUIVALENT TO getCurrentImage
      //DFexp.get_image(this_image,img,true); 
      getCurrentImage(this_image,img);

      //rotationallyAlignOneImage(img(),this_image, opt_refno, opt_psi, opt_flip, maxcorr);

      // Call threads to calculate the rotational alignment of each image in the selfile
      pthread_t * th_ids = (pthread_t *)malloc( threads * sizeof( pthread_t));
      structThreadRotationallyAlignOneImage * threads_d = (structThreadRotationallyAlignOneImage *) 
          malloc ( threads * sizeof( structThreadRotationallyAlignOneImage ) );
      for( int c = 0 ; c < threads ; c++ )
      {
          threads_d[c].thread_id = c;
          threads_d[c].thread_num = threads;
          threads_d[c].prm = this;
          threads_d[c].img=&img();
          threads_d[c].this_image=&this_image;
          threads_d[c].opt_refno=&opt_refno;
          threads_d[c].opt_psi=&opt_psi;
          threads_d[c].opt_flip=&opt_flip; 
          threads_d[c].maxcorr=&maxcorr;
          pthread_create( (th_ids+c), NULL, threadRotationallyAlignOneImage, (void *)(threads_d+c) );
      }

      // Wait for threads to finish and get optimal refno, psi, flip and maxcorr
      for( int c = 0 ; c < threads ; c++ )
      {
          if (*threads_d[c].maxcorr > maxcorr)
          {
              maxcorr = *threads_d[c].maxcorr;
              opt_refno = *threads_d[c].opt_refno;
              opt_psi = *threads_d[c].opt_psi;
              opt_flip = *threads_d[c].opt_flip;
          }
          pthread_join(*(th_ids+c),NULL);
      }

      // Flip order to loop through references
      loop_forward_refs = !loop_forward_refs;

      opt_rot  = XX(mysampling.no_redundant_sampling_points_angles[opt_refno]);
      opt_tilt = YY(mysampling.no_redundant_sampling_points_angles[opt_refno]);
      translationallyAlignOneImage(img(), opt_refno, opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr);

      // Add previously applied translation to the newly found one
      opt_xoff += img.Xoff();
      opt_yoff += img.Yoff();

      // Output
      my_output[imgno * MY_OUPUT_SIZE + 1] = this_image;
      my_output[imgno * MY_OUPUT_SIZE + 2] = opt_rot;
      my_output[imgno * MY_OUPUT_SIZE + 3] = opt_tilt;
      my_output[imgno * MY_OUPUT_SIZE + 4] = opt_psi;
      my_output[imgno * MY_OUPUT_SIZE + 5] = opt_xoff;
      my_output[imgno * MY_OUPUT_SIZE + 6] = opt_yoff;
      my_output[imgno * MY_OUPUT_SIZE + 7] = opt_refno; 
      my_output[imgno * MY_OUPUT_SIZE + 8] = opt_flip;
      my_output[imgno * MY_OUPUT_SIZE + 9] = maxcorr;

  }

}

void Prog_angular_projection_matching_prm::getCurrentImage(int imgno, ImageXmipp &img)
{

    FileName	     fn_img;
    DocLine	     DL;
    Matrix2D<double> A;

    // jump to line imgno+1 in DFexp, get data and filename
    DFexp.locate(imgno+1);
    DL = DFexp.get_current_line();
    DFexp.previous();
    if (DFexp.get_current_line().Is_comment()) 
    {
       fn_img = ((DFexp.get_current_line()).get_text()).erase(0, 3);
    }
    else
    {
       REPORT_ERROR(1,"BUG: no comment in DFexp where expected....");
    }

    // Read actual image
    img.read(fn_img);
    img().setXmippOrigin();

    // Store translation in header and apply it to the actual image
    img.set_Xoff(DL[3]);
    img.set_Yoff(DL[4]);
    img.set_rot(0.);
    img.set_tilt(0.);
    img.set_psi(0.);
    img.set_flip(0.);

    A = img.get_transformation_matrix(true);
    if (!A.isIdentity())
       img().selfApplyGeometryBSpline(A, 3, IS_INV, WRAP);

}

