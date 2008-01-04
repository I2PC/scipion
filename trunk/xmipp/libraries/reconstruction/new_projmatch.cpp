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

// Read arguments ==========================================================
void Prog_new_projection_matching_prm::read(int argc, char **argv)  {

  // Read command line
  if (checkParameter(argc,argv,"-more_options")) { usage(); extended_usage();}
  fn_vol=getParameter(argc,argv,"-vol");
  fn_img=getParameter(argc,argv,"-i");
  fn_root=getParameter(argc,argv,"-o","out");
  sampling=textToFloat(getParameter(argc,argv,"-sam","10"));
  max_shift=textToFloat(getParameter(argc,argv,"-max_shift","5"));

  // Additional commands
  ang_search=textToFloat(getParameter(argc,argv,"-ang_search","1000."));
  tilt_range0=textToFloat(getParameter(argc,argv,"-tilt0","0."));
  tilt_rangeF=textToFloat(getParameter(argc,argv,"-tiltF","90."));
  Ri=textToInteger(getParameter(argc,argv,"-Ri","-1"));
  Ro=textToInteger(getParameter(argc,argv,"-Ro","-1"));
  symmetry = getParameter(argc, argv, "-symmetry", "cn");
  sym_order = textToInteger(getParameter(argc, argv, "-sym_order", "1"));
  output_refs=checkParameter(argc,argv,"-output_refs");
  modify_header=!checkParameter(argc,argv,"-dont_modify_header");
  output_classes=checkParameter(argc,argv,"-output_classes");
  do_reproject=checkParameter(argc,argv,"-save_memA");

  // Hidden stuff
  verb=textToInteger(getParameter(argc,argv,"-verb","1"));

  // TODO!
  if (tilt_range0>0. || tilt_rangeF<90.) 
      REPORT_ERROR(1,"Limited tilt search not implemented yet in this version...");

}

// Show ====================================================================
void Prog_new_projection_matching_prm::show() {

  if (verb>0) {
    std::cerr << "  Input images            : "<< SF.name()<<" ("<<SF.ImgNo()<<")"<<std::endl;
    std::cerr << "  Reference volume        : "<< fn_vol<<std::endl;
    std::cerr << "  Output rootname         : "<< fn_root<<std::endl;
    std::cerr << "  Angular sampling rate   : "<< sampling <<std::endl;
    if (Ri>0)
    std::cerr << "  Inner radius rot-search : "<<Ri<<std::endl;
    if (Ro>0)
    std::cerr << "  Outer radius rot-search : "<<Ro<<std::endl;
    std::cerr << "  -> Limit search of origin offsets to  +/- "<<max_shift<<" pixels"<<std::endl;
    if (ang_search>0) {
      std::cerr << "  -> Limit search of rot and tilt angle to  +/- "<<ang_search<<" degrees"<<std::endl;
    }
    if (tilt_range0>0. || tilt_rangeF<180.)
    {
	std::cerr << "  -> Limited tilt range       : "<<tilt_range0<<"  "<<tilt_rangeF<<std::endl;
    }
    if (!modify_header)
    {
	std::cerr << "  -> Do not modify the image headers (only output docfile)"<<std::endl;
    }
    if (output_refs)
    {
	std::cerr << "  -> Output library projections, sel and docfile"<<std::endl;
    }
    if (output_classes)
    {
	std::cerr << "  -> Output class averages and selfiles for each projection direction "<<std::endl;
    }
    if (do_reproject)
    {
	std::cerr << "  -> save_memA: re-project volume for each translational search"<<std::endl;
    }

    std::cerr << " ================================================================="<<std::endl;
  }
}

// Usage ===================================================================
void Prog_new_projection_matching_prm::usage() {
  std::cerr << "Usage:  projection_matching [options] "<<std::endl;
  std::cerr << "   -i <selfile>                : Selfile with input images \n"
       << "   -vol <volume>               : Reference volume \n"
       << " [ -o <rootname=\"out\"> ]       : Output rootname \n"
       << " [ -sam <float=10> ]           : Sampling rate for rot, tilt & psi (degrees) \n"
       << " [ -max_shift <float=5> ]      : Maximum change in origin offset (+/- pixels) \n"
       << " [ -more_options ]             : Show all program options\n";
}

// Extended usage ===================================================================
void Prog_new_projection_matching_prm::extended_usage() {
  std::cerr << "Additional options: \n"
       << " [ -ang_search <float=-1> ]    : Maximum change in rot & tilt  (+/- degrees) \n"
       << " [ -tilt0 <float=0.> ]         : Lower-value for restricted tilt angle search \n"
       << " [ -tiltF <float=90.> ]        : Higher-value for restricted tilt angle search \n"
       << " [ -Ri <float=1> ]             : Inner radius to limit rotational search \n"
       << " [ -Ro <float=dim/2 - 1> ]     : Outer radius to limit rotational search \n"
       << " [ -symmetry <cn> ]            :  One of the 17 possible symmetries in\n"
       << "                                  single particle electronmicroscopy\n"
       << "                                  i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv, dnh, t, td, th, o, oh, i, ih\n"
       << " [ -sym_order <1> ]            : For infinite groups symmetry order\n" 
       << " [ -output_refs ]              : Output reference projections, sel and docfile \n"
       << " [ -output_classes ]           : Output averages and selfiles for all projection directions\n"
       << " [ -save_memA ]                : Save memory by re-projecting volume for translational searched\n"
       << " [ -dont_modify_header ]       : Do not store alignment parameters in the image headers \n";
  exit(1);
}

// Side info stuff ===================================================================
void Prog_new_projection_matching_prm::produce_Side_info() {

    ImageXmipp       img,empty;
    Projection       proj;
    DocFile          DF;
    SelFile          SFr,emptySF;
    SymList          SL;
    FileName         fn_tmp, fn_refs;
    double           mean,stddev,psi=0.;
    Matrix1D<double> dataline(3);
    int              nl;

    // Read Selfile and get dimensions
    SF.read(fn_img);
    SF.ImgSize(dim,dim);

    // Set ring defaults
    if (Ri<0) Ri=1;
    if (Ro<0) Ro=(dim/2)-1;

    // Initialize empty image
    if (output_classes)
    {
	empty().resize(dim,dim);
	empty().setXmippOrigin();
	empty.clear_header();
    }

    // Set up angular sampling
    mysampling.SetSampling(sampling);
    mysampling.SetNeighborhoodRadius(ang_search);
    mysampling.Compute_sampling_points(true); // STILL ADD MIN&MAX TILT!
    mysampling.create_sym_file(symmetry, sym_order);
    mysampling.remove_redundant_points(symmetry, sym_order);
    mysampling.compute_neighbors();    

    // Read the reference volume
    vol.read(fn_vol);
    vol().setXmippOrigin();

    // Generate reference projections from sampling
    if (verb>0) std::cerr << "--> Projecting the reference volume ..."<<std::endl;

    nl = mysampling.no_redundant_sampling_points_vector.size();
    if (verb>0) init_progress_bar(nl);
    for (int i = 0; i < nl; i++)
    {
	double rot,  tilt,  psi = 0.;
	rot = XX(mysampling.no_redundant_sampling_points_angles[i]);
	tilt = YY(mysampling.no_redundant_sampling_points_angles[i]);
	psi = ZZ(mysampling.no_redundant_sampling_points_angles[i]);

	proj.clear_header();
	project_Volume(vol(),proj,dim,dim,rot,tilt,psi);
	proj().setXmippOrigin();
	if (output_refs)
	{
	    DF.append_angles(rot,tilt,psi,"rot","tilt","psi");
	    fn_tmp.compose(fn_root+"_lib",i+1,"proj");
	    proj.write(fn_tmp);
	    SFr.insert(fn_tmp);
	}
	if (output_classes)
	{
	    empty.rot()=rot;
	    empty.tilt()=tilt;
	    class_avgs.push_back(empty);
	    class_selfiles.push_back(emptySF);
	}

	// Calculate FTs of polar rings (store its complex conjugated!)
	Matrix2D<double>         Maux;
	Polar<double>            P;
	Polar<std::complex <double> > fP;
	produceReverseGriddingMatrix2D(proj(),Maux,kb);
	P.getPolarFromCartesian(Maux,kb,Ri,Ro);
	mean = P.computeSum(true);
	stddev = P.computeSum2(true);
	stddev = sqrt(stddev - mean * mean);
	P -= mean;
	fP = P.fourierTransformRings(true);
	fP_ref.push_back(fP);
	stddev_ref.push_back(stddev);
	// If not save_memA: also store reference projections in cartesian coordinates 
	if (!do_reproject)
	    proj_ref.push_back(proj());
	
	if (verb>0 && (i%MAX(1,nl/60)==0)) progress_bar(i);
    }

    if (output_refs)
    {
	fn_tmp=fn_root+"_lib.doc";
	DF.write(fn_tmp);
	fn_tmp=fn_root+"_lib.sel";
	SFr.write(fn_tmp);
    }

    // If not reprojecting, remove the volume from memory
    if (!do_reproject)
	vol.clear();

    if (verb>0) progress_bar(nl);
    if (verb>0) std::cerr << " ================================================================="<<std::endl;

}

void Prog_new_projection_matching_prm::rotationally_align_one_image(const Matrix2D<double> &img,
								    const int &samplenr, int &opt_samplenr,
								    double &opt_psi, double &opt_flip, double &maxcorr)
{
    Matrix2D<double>         Maux;
    Polar<double>            P;
    Polar<std::complex <double> > fP,fPm;
    Matrix1D<double>         ang,corr;
    int                      max_index;
    double                   mean;

    // Calculate polar coordinates using gridding
    produceReverseGriddingMatrix2D(img,Maux,kb);
    P.getPolarFromCartesian(Maux,kb,Ri,Ro);
    mean = P.computeSum(true);
    P -= mean; // for normalized cross-correlation coefficient
    fP = P.fourierTransformRings(false);
    fPm = P.fourierTransformRings(true);

    // Loop over all relevant "neighbors" (i.e. directions within the search range)
    maxcorr = -99.e19;
    for (int j = 0; j < mysampling.my_neighbors[samplenr].size();j++)
    {
	int i = mysampling.my_neighbors[samplenr][j];

	// A. Check straight image
	rotationalCorrelation(fP,fP_ref[i],ang,corr);
	corr /= stddev_ref[i]; // for normalized ccf, forget about constant stddev_img
	for (int k = 0; k < XSIZE(corr); k++)
	    if (corr(k)> maxcorr)
	    {
		maxcorr = corr(k);
		opt_psi = ang(k);
		opt_samplenr = i;
		opt_flip = 0.;
	    }

	// B. Check mirrored image
	rotationalCorrelation(fPm,fP_ref[i],ang,corr);
	corr /= stddev_ref[i]; // for normalized ccf, forget about constant stddev_img
	for (int k = 0; k < XSIZE(corr); k++)
	    if (corr(k)> maxcorr)
	    {
		maxcorr = corr(k);
		opt_psi = ang(k);
		opt_samplenr = i;
		opt_flip = 1.;
	    }

    }

}

void Prog_new_projection_matching_prm::translationally_align_one_image(const Matrix2D<double> &img,
								       const int &samplenr, const double &psi, 
								       const double &opt_flip, double &opt_xoff, 
								       double &opt_yoff, double &maxcorr)
{
    Projection       proj;
    Matrix2D<double> Mtrans,Mimg,Mref;

    if (do_reproject)
    {
	// Re project volume according to optimal rot, tilt and psi
	double rot = XX(mysampling.no_redundant_sampling_points_angles[samplenr]);
	double tilt = YY(mysampling.no_redundant_sampling_points_angles[samplenr]);
	project_Volume(vol(),proj,dim,dim,rot,tilt,psi);
	Mref=proj();
    }	
    else
	// Rotate stored reference projection by phi degrees
	Mref=proj_ref[samplenr].rotate(-psi,DONT_WRAP);

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

    // Calculate standard cross-correlation coefficient
    Mtrans = Mimg.translate(vectorR2(opt_xoff,opt_yoff));
    maxcorr = correlation_index(Mref,Mtrans);

    // Correct X-shift for mirrored images
    if (opt_flip>0.)
	opt_xoff *= -1.;	

}

void Prog_new_projection_matching_prm::PM_loop_over_all_images(SelFile &SF, DocFile &DFo,
							   double &sumCC) {


  ImageXmipp img;
  FileName fn_img;
  Matrix1D<double> dataline(8);
  double opt_rot,opt_tilt,opt_psi,opt_flip,opt_xoff,opt_yoff,maxcorr;
  int c,nn,imgno,samplenr,opt_samplenr;

  if (verb>0) std::cerr << "--> Projection matching ... "<<std::endl;

  // Initialize
  nn=SF.ImgNo();
  if (verb>0) init_progress_bar(nn);
  c=MAX(1,nn/60);

  // Loop over all images
  sumCC=0.;
  imgno=0;
  SF.go_beginning();
  while ((!SF.eof())) {
    fn_img=SF.NextImg();
    // TODO: maybe instead of reading from the headers I could get the
    // angles from a docfile? Or allow both options?
    // That would allow always using the same images, without touching
    // their headers...
    img.read(fn_img,false,false,false,true);
    img().setXmippOrigin();


    // TODO: I WILL NEED A FUNCTION THAT RELATES THE CURRENT ANGLES TO THE
    // CLOSEST SAMPLE POINT IN MYSAMPLING (FOR NEIGHBOURHOOD SEARCHES!)
    samplenr = 0; // For now no limited searches!

    rotationally_align_one_image(img(), samplenr, opt_samplenr, opt_psi, opt_flip, maxcorr);

    opt_rot  = XX(mysampling.no_redundant_sampling_points_angles[opt_samplenr]);
    opt_tilt = YY(mysampling.no_redundant_sampling_points_angles[opt_samplenr]);
    opt_psi += ZZ(mysampling.no_redundant_sampling_points_angles[opt_samplenr]);

    translationally_align_one_image(img(), opt_samplenr, opt_psi, opt_flip, opt_xoff, opt_yoff, maxcorr);

    opt_xoff += img.Xoff();
    opt_yoff += img.Yoff();

    sumCC+=maxcorr;
    dataline(0) = opt_rot;                 // rot
    dataline(1) = opt_tilt;                // tilt
    dataline(2) = opt_psi;                 // psi
    dataline(3) = opt_xoff;                // Xoff
    dataline(4) = opt_yoff;                // Yoff
    dataline(5) = opt_samplenr+1;          // optimal direction number
    dataline(6) = opt_flip;                // Mirror
    dataline(7) = maxcorr;                 // maxCC
    DFo.append_comment(img.name());
    DFo.append_data_line(dataline);


    if (modify_header || output_classes)
    {
	// Re-read image to get the untransformed image matrix again
	// FIXME: THIS IS RATHER STUPID AND VERY DISK-ACCESS EXPENSIVE!
	// CANT I JUST MODIFY THE HEADER IF (modify_header)
	// AND DO A SECOND GRIDDING INTERPOLATION IF (output_classes)??

	// MAYBE THE OUTPUT_CLASSES OPTION SHOULDNT BE IN THIS PROGRAM
	// ANYWAY. MAYBE A STAND-ALONE PROGRAM THAT MAKES THESE FROM
	// THE OUTPUT DOCFILE WOULD BE A BETTER IDEA...
	img.read(fn_img);
	img.set_eulerAngles(opt_rot,opt_tilt,opt_psi);
	img.set_originOffsets(opt_xoff,opt_yoff);
	img.flip()=opt_flip;
	if (modify_header) 
	    img.write(fn_img);
	if (output_classes)
	{
	    img().selfApplyGeometryBSpline(img.get_transformation_matrix(),3,IS_INV,WRAP);
	    class_avgs[opt_samplenr]()+=img();
	    class_avgs[opt_samplenr].weight()+=1.;
	    class_selfiles[opt_samplenr].insert(img.name());
	}
    }

    if (verb>0) if (imgno%c==0) progress_bar(imgno);
    imgno++;
  }

  if (verb>0) progress_bar(nn);
  if (verb>0) std::cerr << " ================================================================="<<std::endl;

}

void Prog_new_projection_matching_prm::write_classes()
{

    FileName fn_base,fn_img,fn_sel;
    SelFile SF,SF2;

    fn_base = fn_root + "_class";
    SF.clear();
    SF2.clear();
    for (int i = 0; i < mysampling.no_redundant_sampling_points_vector.size(); i++)
    {
	fn_img.compose(fn_base, i + 1, "xmp");
	SF.insert(fn_img);
	fn_sel.compose(fn_base, i + 1, "sel");
	class_avgs[i]() /= class_avgs[i].weight();
	class_avgs[i].write(fn_img);
	class_selfiles[i].write(fn_sel);
    }
    fn_base += "es.sel";
    SF.write(fn_base);
}
