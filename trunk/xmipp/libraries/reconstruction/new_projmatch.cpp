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

#include "new_projmatch.h"
//#define DEBUG
//#define TIMING

// Read arguments ==========================================================
void Prog_new_projection_matching_prm::read(int argc, char **argv)  {

  // Read command line
  if (checkParameter(argc,argv,"-more_options")) { usage(); extendedUsage();}
  fn_ref  = getParameter(argc,argv,"-ref","ref");
  fn_exp  = getParameter(argc,argv,"-i");
  fn_root = getParameter(argc,argv,"-o","out");

  // Additional commands
  Ri=textToInteger(getParameter(argc,argv,"-Ri","-1"));
  Ro=textToInteger(getParameter(argc,argv,"-Ro","-1"));
  search5d_shift  = textToInteger(getParameter(argc,argv,"-search5d_shift","0"));
  search5d_step = textToInteger(getParameter(argc,argv,"-search5d_step","2"));
  max_shift = textToFloat(getParameter(argc,argv,"-max_shift","-1"));
  avail_memory = textToFloat(getParameter(argc,argv,"-mem","1"));
  fn_ctf  = getParameter(argc,argv,"-ctf","");
  phase_flipped = checkParameter(argc, argv, "-phase_flipped");

  // Hidden stuff
  verb=textToInteger(getParameter(argc,argv,"-verb","1"));

}

// Show ====================================================================
void Prog_new_projection_matching_prm::show() {

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
	std::cerr << "  5D-search shift range   : "<<search5d_shift<<" pixels (sampled "<<search5d_xoff.size()<<" times)"<<std::endl;
    }
    if (fn_ctf!="")
    {
	std::cerr << "  CTF parameter file      :  " <<fn_ctf<<std::endl;
	if (phase_flipped)
	    std::cerr << "    + Assuming images have been phase flipped " << std::endl;
	else
	    std::cerr << "    + Assuming images have not been phase flipped " << std::endl;
    }
    std::cerr << " ================================================================="<<std::endl;
  }
}

// Usage ===================================================================
void Prog_new_projection_matching_prm::usage() {
  std::cerr << "Usage:  projection_matching [options] "<<std::endl;
  std::cerr << "   -i <docfile>            : Docfile with input images \n"
	    << "   -ref <ref rootname=\"ref\"> : Rootname for reference projection files \n"
            << " [ -search5d_shift <int=0> ]   : Search range (in +/- pix) for 5D shift search\n"
            << " [ -search5d_step  <int=2> ]   : Step size for 5D shift search (in pix) \n"
	    << " [ -Ri <float=1> ]             : Inner radius to limit rotational search \n"
	    << " [ -Ro <float=dim/2 - 1> ]     : Outer radius to limit rotational search \n"
	    << " [ -more_options ]             : Show all program options\n";
}

// Extended usage ===================================================================
void Prog_new_projection_matching_prm::extendedUsage() {
  std::cerr << "Additional options: \n"
            << " [ -mem <float=1> ]            : Available memory for reference library (Gb)\n"
	    << " [ -max_shift <float=-1> ]     : Max. change in origin offset (+/- pixels; neg= no limit) \n"
	    << " [ -ctf <ctfparam-file> ]      : Apply this CTF to the reference projections \n"
	    << " [ -phase_flipped ]            : Use this if the experimental images have been phase flipped\n";
  exit(1);
}

// Side info stuff ===================================================================
void Prog_new_projection_matching_prm::produceSideInfo() {

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

    // Read one image to get dim
    DFexp.go_first_data_line();
    DFexp.previous();
    if (DFexp.get_current_line().Is_comment()) 
	fn_img = ((DFexp.get_current_line()).get_text()).erase(0, 3);
    else
	REPORT_ERROR(1,"BUG: no comment in DFexp where expected....");
    img.read(fn_img);
    dim = XSIZE(img());

    // Set max_shift
    if (max_shift<0) max_shift = dim/2;

    // Set ring defaults
    if (Ri<1) Ri=1;
    if (Ro<0) Ro=(dim/2)-1;

    // Calculate necessary memory per image
    img().produceSplineCoefficients(Maux,3);
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro);
    P.fourierTransformRings(true);
    double memory_per_ref = 0.;
    for (int i = 0; i < fP.getRingNo(); i++)
    {
	memory_per_ref += (double) fP.getSampleNo(i) * 2 * sizeof(double);
    }
    memory_per_ref += dim * dim * sizeof(double);
    max_nr_refs_in_memory = ROUND( 1024 * 1024 * 1024 * avail_memory / memory_per_ref);

    // Set up angular sampling
    mysampling.read_sampling_file(fn_ref,false);
    total_nr_refs = mysampling.no_redundant_sampling_points_angles.size();

    // Don't reserve more memory than necessary
    max_nr_refs_in_memory = XMIPP_MIN(max_nr_refs_in_memory, total_nr_refs);

    // Initialize pointers for reference retrieval 
    pointer_allrefs2refsinmem.resize(total_nr_refs,-1);
    pointer_refsinmem2allrefs.resize(max_nr_refs_in_memory,-1);
    counter_refs_in_memory = 0;

    // Initialize vectors with references
    Polar<std::complex <double> > fP_dum;
    fP_ref.resize(max_nr_refs_in_memory,fP_dum);
    Matrix2D<double> Mdum;
    proj_ref.resize(max_nr_refs_in_memory,Mdum);
    double dum;
    stddev_ref.resize(max_nr_refs_in_memory,dum);
    loop_forward_refs=true;

    // Initialize 5D search vectors
    search5d_xoff.clear();
    search5d_yoff.clear();
    // Make sure origin is included
    int myfinal=search5d_shift + search5d_shift%search5d_step;
    for (int xoff = -myfinal; xoff <= myfinal; xoff+= search5d_step)
    {
	for (int yoff = -myfinal; yoff <= myfinal; yoff+= search5d_step)
	{
	    // Only take a circle (not a square)
	    if ( xoff*xoff + yoff*yoff <= search5d_shift*search5d_shift)
	    {
		search5d_xoff.push_back(xoff);
		search5d_yoff.push_back(yoff);
	    }
	}
    }

    // CTF stuff
    if (fn_ctf != "")
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
	ctf.Generate_CTF(dim, dim, ctfmask);
	Mctf.resize(dim,dim);
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mctf)
	{
	    if (phase_flipped) dMij(Mctf, i, j) = fabs(dMij(ctfmask, i, j).real());
	    else dMij(Mctf, i, j) = dMij(ctfmask, i, j).real();
	}
    }

}

void Prog_new_projection_matching_prm::getCurrentImage(int imgno, ImageXmipp &img)
{

    FileName         fn_img;
    DocLine          DL;
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
    img.Xoff() = DL[3];
    img.Yoff() = DL[4];
    img.rot()  = 0.;
    img.tilt() = 0.;
    img.psi()  = 0.;
    img.flip() = 0.;

    A = img.get_transformation_matrix(true);
    if (!A.isIdentity())
	img().selfApplyGeometryBSpline(A, 3, IS_INV, WRAP);

}

int Prog_new_projection_matching_prm::getCurrentReference(int refno)
{

    FileName                      fnt;
    ImageXmipp                    img;
    double                        mean,stddev;
    Matrix2D<double>              Maux;
    Polar<double>                 P;
    Polar<std::complex <double> > fP;

    // Image was not stored yet: read it from disc and store
    fnt.compose(fn_ref,refno+1,"xmp");
    img.read(fnt);
    img().setXmippOrigin();

    // Apply CTF (this takes approx as long as calculating the polar
    // transform etc.)
    if (fn_ctf!="")
    {
	Matrix2D<std::complex<double> > Faux;
	FourierTransform(img(),Faux);
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mctf)
	{
	    dMij(Faux,i,j) *= dMij(Mctf,i,j);
	}
	InverseFourierTransform(Faux,img());
    }

    // Calculate FTs of polar rings and its stddev
    img().produceSplineCoefficients(Maux,3);
    P.getPolarFromCartesianBSpline(Maux,Ri,Ro);
    mean = P.computeSum(true);
    stddev = P.computeSum2(true);
    stddev = sqrt(stddev - mean * mean);
    P -= mean;
    fP = P.fourierTransformRings(true);

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

}

void Prog_new_projection_matching_prm::rotationallyAlignOneImage(Matrix2D<double> &img,
								 int imgno, 
								 int &opt_refno,
								 double &opt_psi, 
								 double &opt_flip, 
								 double &maxcorr)
{
    Matrix2D<double>         Maux;
    Polar<double>            P;
    Matrix1D<double>         ang,corr;
    int                      max_index, refno, myinit, myfinal, myincr, nr_trans;
    double                   mean, stddev;
    std::vector<double>      stddev_img;
    Polar<std::complex <double> > fP,fPm;
    std::vector< Polar <std::complex <double> > > fP_img,fPm_img;

#ifdef TIMING
    TimeStamp t0,t1,t2; 
    time_config();
    annotate_time(&t0);
    annotate_time(&t2);
#endif

    maxcorr = -99.e99;
    img.produceSplineCoefficients(Maux,3);
    // Precalculate polar transform of each translation
    nr_trans = search5d_xoff.size();
    for (int itrans = 0; itrans < nr_trans; itrans++)
    {
	P.getPolarFromCartesianBSpline(Maux,Ri,Ro,(double)search5d_xoff[itrans],(double)search5d_yoff[itrans]);
	mean = P.computeSum(true);
	stddev = P.computeSum2(true);
	stddev = sqrt(stddev - mean * mean);
	P -= mean; // for normalized cross-correlation coefficient
	fP = P.fourierTransformRings(false);
	fPm = P.fourierTransformRings(true);
	fP_img.push_back(fP);
	fPm_img.push_back(fPm);
	stddev_img.push_back(stddev);
    }

#ifdef TIMING
    float prepare_img = elapsed_time(t0);
    float get_refs = 0.;
    annotate_time(&t0);
#endif

    // Switch the order of looping through the references every time.
    // That way, in case max_nr_refs_in_memory<total_nr_refs
    // the references read in memory for the previous image
    // will still be there when processing the next image
    if (loop_forward_refs)
    {
	myinit = 0;
	myfinal = mysampling.my_neighbors[imgno].size();
	myincr = 1;
    }
    else
    {
	myinit = mysampling.my_neighbors[imgno].size()-1;
	myfinal = -1;
	myincr = -1;
    }

    // Loop over all relevant "neighbours" (i.e. directions within the search range)
    for (int i = myinit; i != myfinal; i+=myincr)
    {

#ifdef TIMING
	annotate_time(&t1);
#endif

	// Get pointer to the current reference image
	refno = pointer_allrefs2refsinmem[mysampling.my_neighbors[imgno][i]];
	if (refno == -1)
	{
	    // Reference is not stored in memory (anymore): (re-)read from disc
	    getCurrentReference(mysampling.my_neighbors[imgno][i]);
	    refno = pointer_allrefs2refsinmem[mysampling.my_neighbors[imgno][i]];
	}

#ifdef TIMING
	get_refs += elapsed_time(t1);
#endif

#ifdef DEBUG
	std::cerr<<"Got refno= "<<refno<<" pointer= "<<mysampling.my_neighbors[imgno][i]<<std::endl;
#endif	    


	// Loop over all 5D-search translations
	for (int itrans = 0; itrans < nr_trans; itrans++)
	{

	    // A. Check straight image
	    rotationalCorrelation(fP_img[itrans],fP_ref[refno],ang,corr);
	    corr /= stddev_ref[refno] * stddev_img[itrans]; // for normalized ccf
	    for (int k = 0; k < XSIZE(corr); k++)
	    {
		if (corr(k)> maxcorr)
		{
		    maxcorr = corr(k);
		    opt_psi = ang(k);
		    opt_refno = mysampling.my_neighbors[imgno][i];
		    opt_flip = 0.;
		}
	    }

#ifdef DEBUG
	    std::cerr<<"straight: corr "<<maxcorr<<std::endl;
#endif	    

	    // B. Check mirrored image
	    rotationalCorrelation(fPm_img[itrans],fP_ref[refno],ang,corr);
	    corr /= stddev_ref[refno] * stddev_img[itrans]; // for normalized ccf
	    for (int k = 0; k < XSIZE(corr); k++)
	    {
		if (corr(k)> maxcorr)
		{
		    maxcorr = corr(k);
		    opt_psi = ang(k);
		    opt_refno = mysampling.my_neighbors[imgno][i];
		    opt_flip = 1.;
		}
	    }

#ifdef DEBUG
	    std::cerr<<"mirror: corr "<<maxcorr;
	    if (opt_flip==1.) std::cerr<<"**";
	    std::cerr<<std::endl;
#endif	    

	}
    }

    // Flip order to loop through references
    loop_forward_refs = !loop_forward_refs;

#ifdef TIMING
    float all_rot_align = elapsed_time(t0);
    float total_rot = elapsed_time(t2);
    std::cerr<<" rotal%% "<<total_rot
	     <<" => prep: "<<prepare_img
	     <<" all_refs: "<<all_rot_align
	     <<" (of which "<<get_refs
	     <<" to get "<< mysampling.my_neighbors[imgno].size()
	     <<" refs for imgno "<<imgno<<" )"
	     <<std::endl;
#endif

}

void Prog_new_projection_matching_prm::translationallyAlignOneImage(Matrix2D<double> &img,
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
	getCurrentReference(opt_refno);
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



void Prog_new_projection_matching_prm::processSomeImages(int * my_images, double * my_output) 
{


  ImageXmipp img;
  double opt_rot, opt_tilt, opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr;
  int opt_refno;

  int nr_images = my_images[0];
  // Prepare my_output
  my_output[0] = nr_images * MY_OUPUT_SIZE;

  // Loop over all images
  for (int imgno = 0; imgno < nr_images; imgno++)
  {

      // add one because first number is number of elements in the array
      int this_image = my_images[imgno + 1]; 

      // read experimental image and all corresponding references in memory
      getCurrentImage(this_image,img);

      // Align the image (for now 3D+2D search only)
      rotationallyAlignOneImage(img(),this_image, opt_refno, opt_psi, opt_flip, maxcorr);
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

