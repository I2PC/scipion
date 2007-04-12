/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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
#include "../Prog_projection_matching.hh"

// Read arguments ==========================================================
void Prog_projection_matching_prm::read(int argc, char **argv)  {

  // Read command line
  if (check_param(argc,argv,"-show_all_options")) { usage(); extended_usage();}
  fn_vol=get_param(argc,argv,"-vol","");
  SF.read(get_param(argc,argv,"-i"));
  SF.ImgSize(dim,dim);
  fn_root=get_param(argc,argv,"-o","out");
  sampling=AtoF(get_param(argc,argv,"-sam","10"));
  max_shift=AtoF(get_param(argc,argv,"-max_shift","5"));

  // Additional commands
  ang_search=AtoF(get_param(argc,argv,"-ang_search","-1"));
  Ri=AtoF(get_param(argc,argv,"-Ri","-1"));
  Ro=AtoF(get_param(argc,argv,"-Ro","-1"));
  fn_sym=get_param(argc,argv,"-sym","");
  fn_ang=get_param(argc,argv,"-ang","");
  output_refs=check_param(argc,argv,"-output_refs");
  modify_header=!check_param(argc,argv,"-dont_modify_header");
  fn_ref=get_param(argc,argv,"-ref","");
  output_classes=check_param(argc,argv,"-output_classes");
  tilt_range0=AtoF(get_param(argc,argv,"-tilt0","0."));
  tilt_rangeF=AtoF(get_param(argc,argv,"-tiltF","180."));

  // Checks
  if (fn_ref=="" && fn_vol=="") 
    REPORT_ERROR(1," Provide either -vol or -ref!");
  if (fn_ang!="" && ang_search>0)
    REPORT_ERROR(1," option -ang and -ang_search are incompatible!");

  // Hidden stuff
  verb=AtoI(get_param(argc,argv,"-verb","1"));

}

// Show ====================================================================
void Prog_projection_matching_prm::show() {

  if (verb>0) {
    cerr << "  Input images            : "<< SF.name()<<" ("<<SF.ImgNo()<<")"<<endl;
    cerr << "  Reference volume        : "<< fn_vol<<endl;
    cerr << "  Output rootname         : "<< fn_root<<endl;
    cerr << "  Angular sampling rate   : "<< sampling <<endl;
    if (Ri>0) 
    cerr << "  Inner radius rot-search : "<<Ri<<endl;
    if (Ro>0) 
    cerr << "  Outer radius rot-search : "<<Ro<<endl;
    cerr << "  -> Limit search of origin offsets to  +/- "<<max_shift<<" pixels"<<endl;
    if (ang_search>0) {
      cerr << "  -> Limit search of rot and tilt angle to  +/- "<<ang_search<<" degrees"<<endl;
    }
    if (tilt_range0>0. || tilt_rangeF<180.)
    {
	cerr << "  -> Limited tilt range       : "<<tilt_range0<<"  "<<tilt_rangeF<<endl;
    }
    if (fn_sym!="") 
    {
	cerr << "  -> Limit angular search to asymmetric part, as defined by: "<<fn_sym<<endl;
    }
    if (fn_ang!="") 
    {
	cerr << "  -> Document file with angles for projection library: "<<fn_ang<<endl;
    }
    if (!modify_header) 
    {
	cerr << "  -> Do not modify the image headers (only output docfile)"<<endl;
    }
    if (output_refs) 
    {
	cerr << "  -> Output library projections, sel and docfile"<<endl;
    }
    if (output_classes) 
    {
	cerr << "  -> Output class averages and selfiles for each projection direction "<<endl;
    }

    cerr << " ================================================================="<<endl;
  }
} 

// Usage ===================================================================
void Prog_projection_matching_prm::usage() {
  cerr << "Usage:  projection_matching [options] "<<endl;
  cerr << "   -i <selfile>                : Selfile with input images \n"
       << "   -vol <volume>               : Reference volume \n"
       << " [ -o <rootname=\"out\"> ]       : Output rootname \n"
       << " [ -sam <float=10> ]           : Sampling rate for rot, tilt & psi (degrees) \n"
       << " [ -max_shift <float=5> ]      : Maximum change in origin offset (+/- pixels) \n"
       << " [ -show_all_options ]         : Show all program options\n";
}

// Extended usage ===================================================================
void Prog_projection_matching_prm::extended_usage() {
  cerr << "Additional options: \n"
       << " [ -ang_search <float=-1> ]    : Maximum change in rot & tilt  (+/- degrees) \n"
       << " [ -tilt0 <float=0.> ]         : Lower-value for restricted tilt angle search \n"
       << " [ -tiltF <float=180.> ]       : Higher-value for restricted tilt angle search \n"
       << " [ -Ri <float=0> ]             : Inner radius to limit rotational search \n"
       << " [ -Ro <float=dim/2> ]         : Outer radius to limit rotational search \n"
       << " [ -sym <symfile> ]            : Limit angular search to asymmetric part \n"
       << " [ -output_refs ]              : Output reference projections, sel and docfile \n"
       << " [ -output_classes ]           : Output averages and selfiles for all projection directions\n"
       << " [ -ref <selfile>  ]           : Selfile with reference projections (instead of volume)\n"
       << " [ -ang <docfile> ]            : Angles for projection library (instead of  \n"
       << " [ -dont_modify_header ]       : Do not store alignment parameters in the image headers \n";
  exit(1);
}

// Side info stuff ===================================================================
void Prog_projection_matching_prm::produce_Side_info() {

    VolumeXmipp      vol;
    ImageXmipp       img,empty;
    Projection       proj;
    DocFile          DF,DFi,DF2;
    SelFile          SFr,emptySF;
    SymList          SL;
    FileName         fn_tmp, fn_refs;
    double           mean_ref,stddev_ref,dummy,psi=0.;
    matrix1D<double> dataline(3);
    int              nl;

    // Set nr_psi
    nr_psi=CEIL(360./sampling);
 
    // Create rotational-search mask
    rotmask.resize(dim,dim);
    rotmask.set_Xmipp_origin();
    if (Ri<0.) Ri=0.;
    if (Ro<0.) Ro=(double)dim/2;
    BinaryCrownMask(rotmask,Ri,Ro,INNER_MASK);
    nr_pixels_rotmask=(int)rotmask.sum();

    // Initialize empty image
    if (output_classes) 
    {
	empty().resize(dim,dim);
	empty().set_Xmipp_origin();
	empty.clear_header();
    }
    // Read symmetry file into memory
    if (fn_sym!="") 
    { 
	SL.read_sym_file(fn_sym);
    }

    if (fn_ref!="") 
    // Read projections from selfile
    {
	SFr.read(fn_ref);
	nl=SFr.ImgNo();
	ref_img.clear();
	ref_rot=(double*)malloc(nl*sizeof(double));
	ref_tilt=(double*)malloc(nl*sizeof(double));
	ref_mean=(double*)malloc(nl*sizeof(double));
	ref_stddev=(double*)malloc(nl*sizeof(double));
	SFr.go_beginning();
	nr_dir=0;
	while (!SFr.eof()) 
	{
	    proj.read(SFr.NextImg());
	    proj().set_Xmipp_origin();
	    ref_rot[nr_dir]=proj.rot();
	    ref_tilt[nr_dir]=proj.tilt();
	    compute_stats_within_binary_mask(rotmask,proj(),dummy,dummy,mean_ref,stddev_ref);
	    proj()-=mean_ref;
	    apply_binary_mask(rotmask,proj(),proj(),0.);
	    ref_img.push_back(proj());
	    ref_stddev[nr_dir]=stddev_ref;
	    ref_mean[nr_dir]=mean_ref;
	    nr_dir++;
	    if (output_classes)
	    {
		empty.rot()=proj.rot();
		empty.tilt()=proj.tilt();
		class_avgs.push_back(empty);
		class_selfiles.push_back(emptySF);
	    }
	}
    } 
    else 
    // Generate reference projections from sampling
    {

	if (ang_search>=0) 
	{
	    // ignore -sym or -ref option and use -sam to generate all projections on the Ewald sphere
	    // then select only those that are within the search_ranges of all experimental projections
	    // THIS CODE IS AT LEAST SLOPPY...
      
	    // 1. Create even distribution over the entire Ewald sphere
	    if (verb>0) cerr << "--> Making even distribution on entire Ewald sphere "<<endl;
	    make_even_distribution(DF,sampling,SL,true);
	    if (tilt_range0>0. || tilt_rangeF<180.) 
		limit_tilt_range(DF,tilt_range0,tilt_rangeF);
	    // 2. Get all angles from all experimental images
	    double act_rot_range, ref_rot,ref_tilt,img_rot,img_tilt;
	    int nn,c;
	    SF.go_beginning();
	    DFi.clear();
	    DF2.clear();
	    while (!SF.eof()) 
	    {
		proj.read(SF.NextImg());
		dataline(0)=proj.rot();
		dataline(1)=proj.tilt();
		dataline(2)=proj.psi();
		DFi.append_data_line(dataline);
	    }
	    // 3. Check which angles of DF to use
	    if (verb>0) cerr << "--> Selecting relevant library projection directions ..."<<endl;
	    if (verb>0) 
	    {
		nn=DF.dataLineNo();
		init_progress_bar(nn);
		c=MAX(1,nn/60);
	    }
	    DF.go_first_data_line();
	    int ii=0;
	    while (!DF.eof()) 
	    {
		ref_rot=DF(0);
		ref_tilt=DF(1);
		// act_rot_range is tilt-angle dependent!
		if (ref_tilt>0 && ref_tilt<180) 
		    act_rot_range=ang_search/sin(DEG2RAD(ref_tilt)); 
		else 
		    act_rot_range=361.;
		bool search=false;
		DFi.go_first_data_line();
		while (!DFi.eof()) 
		{
		    img_rot=DFi(0);
		    img_tilt=DFi(1);
		    if ( ABS(realWRAP(img_rot-ref_rot,-180.,180.)) <= act_rot_range && 
			 ABS(realWRAP(img_tilt-ref_tilt,-180.,180.)) <= ang_search ) 
		    {
			search=true;
			break;
		    }
		    DFi.next_data_line();
		}
		if (search) 
		{
		    dataline(0)=ref_rot;
		    dataline(1)=ref_tilt;
		    dataline(2)=0.;
		    DF2.append_data_line(dataline);
		}
		ii++;
		if (verb>0) if (ii%c==0) progress_bar(ii);
		DF.next_data_line();
	    }
	    if (verb>0) progress_bar(nn);
	    DF=DF2;
	    DF2.clear();
	} 
	else if (fn_ang!="") 
	// Generate reference projections from docfile
	{
	    DF.read(fn_ang);
	} 
	else 
	// Generate reference projections from even distribution
	{
	    // Create evenly-distributed reference projection angles
	    if (verb>0) cerr << "--> Making even angular distribution ..."<<endl;
	    make_even_distribution(DF,sampling,SL,true);
	    if (tilt_range0>0. || tilt_rangeF<180.) 
		limit_tilt_range(DF,tilt_range0,tilt_rangeF);
	}

	// At this point we have a docfile with all projection directions
	// Now create reference projection images
	vol.read(fn_vol);
	vol().set_Xmipp_origin();
	nl=DF.dataLineNo();
	ref_img.clear();
	ref_rot=(double*)malloc(nl*sizeof(double));
	ref_tilt=(double*)malloc(nl*sizeof(double));
	ref_mean=(double*)malloc(nl*sizeof(double));
	ref_stddev=(double*)malloc(nl*sizeof(double));
	SFr.reserve(nl);
	SFr.go_beginning();
	DF.go_beginning();

	if (verb>0) cerr << "--> Projecting the reference volume ..."<<endl;
	if (verb>0) init_progress_bar(nl);

	fn_refs=fn_root+"_lib";
	DF.adjust_to_data_line();
	nr_dir=0;
	while (!DF.eof()) 
	{
	    ref_rot[nr_dir]=DF(0);
	    ref_tilt[nr_dir]=DF(1);
	    project_Volume(vol(),proj,dim,dim,ref_rot[nr_dir],ref_tilt[nr_dir],psi);
	    if (output_refs) 
	    {
		fn_tmp.compose(fn_refs,nr_dir+1,"proj");
		proj.write(fn_tmp);
		SFr.insert(fn_tmp);
	    }
	    if (output_classes)
	    {
		empty.rot()=DF(0);
		empty.tilt()=DF(1);
		class_avgs.push_back(empty);
		class_selfiles.push_back(emptySF);
	    }
	    compute_stats_within_binary_mask(rotmask,proj(),dummy,dummy,mean_ref,stddev_ref);
	    proj()-=mean_ref;
	    ref_img.push_back(proj());
	    ref_stddev[nr_dir]=stddev_ref;
	    ref_mean[nr_dir]=mean_ref;
	    DF.next_data_line();
	    nr_dir++;
	    if (verb>0 && (nr_dir%MAX(1,nl/60)==0)) progress_bar(nr_dir);
	}
	if (output_refs) 
	{
	    fn_tmp=fn_refs+".doc";
	    DF.write(fn_tmp);
	    fn_tmp=fn_refs+".sel";
	    SFr.write(fn_tmp);
	}
	if (verb>0) progress_bar(nl);
	if (verb>0) cerr << " ================================================================="<<endl;

    }


}


void Prog_projection_matching_prm::PM_process_one_image(matrix2D<double> &Mexp,
							float &img_rot, float &img_tilt, float &img_psi, 
							int &opt_dirno, double &opt_psi,
							double &opt_xoff, double &opt_yoff, 
							double &maxCC, double &Zscore) {


  // Rotational search ====================================================
  matrix2D<double> Mimg,Mref,Maux,Mcorr;
  double act_rot_range,psi,thisCC,oldCC,aveCC=0.,varCC=0.;
  double stddev_img,mean_img,dummy,xmax,ymax;
  int c=0,ioptpsi=0,ioptflip=0;
  bool search;
  vector<matrix2D<double> >::iterator ipp;

  maxCC=-99.e99; 
  Mimg.resize(dim,dim);
  Mimg.set_Xmipp_origin();
  Mref.resize(dim,dim);
  Mref.set_Xmipp_origin();
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();

  // Calculate mean_img,stddev_img and apply rotmask
  Maux=Mexp;
  compute_stats_within_binary_mask(rotmask,Mexp,dummy,dummy,mean_img,stddev_img);
  Maux-=mean_img;
  apply_binary_mask(rotmask,Maux,Maux,0.);

  // Calculate correlation coefficients for all angles
  FOR_ALL_ROTATIONS() {
    psi=(double)(ipsi*360./nr_psi);
    Mimg=Maux.rotate(psi,DONT_WRAP);
    ipp=ref_img.begin();
    FOR_ALL_DIRECTIONS() {
      search=true;
      // For some strange reason I need to access the vector via its pointer
      // otherwise it goes 50x slower on jumilla (Alpha-Unix)
      Mref=*(ipp);
      if (ang_search >0) {
	// act_rot_range is tilt-angle dependent!
	if (ref_tilt[dirno]>0 && ref_tilt[dirno]<180) act_rot_range=ang_search/sin(DEG2RAD(ref_tilt[dirno])); 
	else act_rot_range=361.;
	if ( ABS(realWRAP(img_rot-ref_rot[dirno],-180.,180.)) > act_rot_range || 
	     ABS(realWRAP(img_tilt-ref_tilt[dirno],-180.,180.)) > ang_search ) search=false;
      }
      if (search) {
	thisCC=0.;
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mimg) {
	  thisCC+=dMij(Mref,i,j)*dMij(Mimg,i,j);
	}
	thisCC/=ref_stddev[dirno]*stddev_img*nr_pixels_rotmask;
	c++;
        oldCC=aveCC;
        aveCC+=(thisCC-oldCC)/(c+1);
        if (c>1) varCC=(1.-1./(double)c)*varCC+(c+1.)*(oldCC-aveCC)*(oldCC-aveCC);
	//cout << "rot= "<<ref_rot[dirno]<<" tilt= "<<ref_tilt[dirno]<<" psi= "<<psi<<" CC= "<<thisCC<<" maxCC= "<<maxCC<<" "<<dirno<<" "<<opt_dirno<<endl;
	if (thisCC>maxCC) {
	  maxCC=thisCC;
	  opt_psi=psi;
	  opt_dirno=dirno;
	}
      }
      ipp++;
    }
  }
                                                                                
  // Calculate Z-score on rotational permutations
  Zscore=(maxCC-aveCC)/(sqrt(varCC));

  // Interpolated translational search for optimal angles ===================================
  Mimg=Mexp.rotate(opt_psi,DONT_WRAP);
  Mref=ref_img[opt_dirno];
  if (max_shift>0) best_shift(Mimg,Mref,xmax,ymax);
  else xmax=ymax=0.;
  if (xmax*xmax+ymax*ymax>max_shift*max_shift) {
    xmax=0.; 
    ymax=0.;
  } 
  opt_xoff=-xmax*COSD(opt_psi)-ymax*SIND(opt_psi);
  opt_yoff=xmax*SIND(opt_psi)-ymax*COSD(opt_psi);

  // Calculate optimal correlation coefficient
  Mimg=Mimg.translate(vector_R2(-xmax,-ymax)); 
  Mimg-=mean_img;
  maxCC=0.;
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mimg) {
    maxCC+=dMij(Mref,i,j)*dMij(Mimg,i,j);
  }
  maxCC/=ref_stddev[opt_dirno]*stddev_img*nr_pixels_rotmask;

}

void Prog_projection_matching_prm::PM_loop_over_all_images(SelFile &SF, DocFile &DFo, 
							   double &sumCC) {


  ImageXmipp img;
  FileName fn_img;
  matrix1D<double> dataline(8);
  matrix2D<double> A(3,3);
  double opt_psi,opt_xoff,opt_yoff,maxCC,Zscore;
  int c,nn,imgno,opt_dirno;

  if (verb>0) cerr << "--> Projection matching ... "<<endl;

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
    img.read(fn_img,false,false,false,true);
    img().set_Xmipp_origin();

    // Perform the projection matching for each image separately
    PM_process_one_image(img(),img.Phi(),img.Theta(),img.Psi(),
			 opt_dirno,opt_psi,opt_xoff,opt_yoff,maxCC,Zscore);

    opt_xoff+=img.Xoff();
    opt_yoff+=img.Yoff();

    sumCC+=maxCC;
    dataline(0)=ref_rot[opt_dirno];      // rot
    dataline(1)=ref_tilt[opt_dirno];     // tilt
    dataline(2)=opt_psi;                 // psi
    dataline(3)=opt_xoff;                // Xoff
    dataline(4)=opt_yoff;                // Yoff
    dataline(5)=opt_dirno+1;             // optimal direction number
    dataline(6)=maxCC;                   // maximum CC
    dataline(7)=Zscore;                   // maximum CC
    DFo.append_comment(img.name());
    DFo.append_data_line(dataline);

    

    if (modify_header) 
    {
	// Re-read image to get the untransformed image matrix again
	img.read(fn_img);
	img.set_eulerAngles(ref_rot[opt_dirno],ref_tilt[opt_dirno],opt_psi);
	img.set_originOffsets(opt_xoff,opt_yoff);
	img.write(fn_img);
    }
    if (output_classes)
    {
	// Re-read image to get the untransformed image matrix again
	img.read(fn_img);
	img.set_eulerAngles(ref_rot[opt_dirno],ref_tilt[opt_dirno],opt_psi);
	img.set_originOffsets(opt_xoff,opt_yoff);
	img().self_apply_geom_Bspline(img.get_transformation_matrix(),3,IS_INV,WRAP);
	class_avgs[opt_dirno]()+=img();
	class_avgs[opt_dirno].weight()+=1.;
	class_selfiles[opt_dirno].insert(img.name());
    }


    if (verb>0) if (imgno%c==0) progress_bar(imgno);
    imgno++;
  }

  if (verb>0) progress_bar(nn);
  if (verb>0) cerr << " ================================================================="<<endl;

  // free memory
  free(ref_mean);
  free(ref_rot);
  free(ref_tilt);
  ref_img.clear();
}
void Prog_projection_matching_prm::write_classes() 
{

    FileName fn_base,fn_img,fn_sel;
    SelFile SF,SF2;

    fn_base=fn_root+"_class";
    SF.clear();
    SF2.clear();
    FOR_ALL_DIRECTIONS() {
	fn_img.compose(fn_base,dirno+1,"xmp");
	SF.insert(fn_img);
	fn_sel.compose(fn_base,dirno+1,"sel");
	class_avgs[dirno]()/=class_avgs[dirno].weight();
	class_avgs[dirno].write(fn_img);
	class_selfiles[dirno].write(fn_sel);
    }
    fn_base+="es.sel";
    SF.write(fn_base);
}
