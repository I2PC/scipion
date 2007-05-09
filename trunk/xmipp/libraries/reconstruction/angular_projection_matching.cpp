/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
 *             Roberto Marabini
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

// Read arguments ==========================================================
void Prog_projection_matching_prm::read(int argc, char **argv)  {

  // Read command line
  if (check_param(argc,argv,"-show_all_options")) { usage(); extended_usage();}
  fn_vol=get_param(argc,argv,"-vol");
  SF.read(get_param(argc,argv,"-i"));
  SF.ImgSize(dim,dim);
  fn_root=get_param(argc,argv,"-o","out");
  sampling=AtoF(get_param(argc,argv,"-sam","10"));
  max_shift=AtoF(get_param(argc,argv,"-max_shift","5"));

  // Additional commands
  ang_search=AtoF(get_param(argc,argv,"-ang_search","360."));
  Ri=AtoF(get_param(argc,argv,"-Ri","-1"));
  Ro=AtoF(get_param(argc,argv,"-Ro","-1"));
  symmetry=get_param(argc,argv,"-symmetry","cn");
  sym_order=AtoI(get_param(argc,argv,"-sym_order","1"));
  output_refs=check_param(argc,argv,"-output_refs");
  modify_header=!check_param(argc,argv,"-dont_modify_header");
  output_classes=check_param(argc,argv,"-output_classes");
  tilt_range0=AtoF(get_param(argc,argv,"-tilt0","0."));
  tilt_rangeF=AtoF(get_param(argc,argv,"-tiltF","180."));

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
    cerr    << "symmetry group:     " << symmetry << endl;
    cerr    << "symmetry order:     " << sym_order << endl;
    if (ang_search>0) {
      cerr << "  -> Limit search of rot and tilt angle to  +/- "<<ang_search<<" degrees"<<endl;
    }
    if (tilt_range0>0. || tilt_rangeF<180.)
    {
	cerr << "  -> Limited tilt range       : "<<tilt_range0<<"  "<<tilt_rangeF<<endl;
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
       << " [ -symmetry <cn> ]            :  One of the 17 possible symmetries in\n" 
       << "                                  single particle electronmicroscopy\n"
       << "                                  i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv, dnh, t, td, th, o, oh, i, ih\n"
       << " [ -sym_order <1> ]            : For infinite groups symmetry order\n"
       << " [ -output_refs ]              : Output reference projections, sel and docfile \n"
       << " [ -output_classes ]           : Output averages and selfiles for all projection directions\n"
       << " [ -dont_modify_header ]       : Do not store alignment parameters in the image headers \n";
  exit(1);
}

// Side info stuff ===================================================================
void Prog_projection_matching_prm::produce_Side_info() {

    VolumeXmipp      vol;
    ImageXmipp       img,empty;
    Projection       proj;
    DocFile          DF;
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

    // Set up angular sampling	
    mysampling.SetSampling(sampling); // STILL ADD OPTION TO USE A SELFILE
    mysampling.SetNeighborhoodRadius(ang_search);
    mysampling.Compute_sampling_points(false); // STILL ADD MIN&MAX TILT!
    mysampling.create_sym_file(symmetry,sym_order);
    mysampling.remove_redundant_points(symmetry,sym_order);
    mysampling.compute_neighbors();

    // Generate reference projections from sampling
    vol.read(fn_vol);
    vol().set_Xmipp_origin();
    nl=mysampling.no_redundant_sampling_points_vector.size();
    ref_img.clear();
    ref_mean=(double*)malloc(nl*sizeof(double));
    ref_stddev=(double*)malloc(nl*sizeof(double));
    SFr.reserve(nl);
    SFr.go_beginning();

    if (verb>0) cerr << "--> Projecting the reference volume ..."<<endl;
    if (verb>0) init_progress_bar(nl);

    fn_refs=fn_root+"_lib";
    for (int i=0; i<mysampling.no_redundant_sampling_points_vector.size(); i++)
    {
	double rot=XX(mysampling.no_redundant_sampling_points_angles[i]);
	double tilt=YY(mysampling.no_redundant_sampling_points_angles[i]);

	project_Volume(vol(),proj,dim,dim,rot,tilt,psi);
	if (output_refs)
	{
	    fn_tmp.compose(fn_refs,i+1,"proj");
	    proj.write(fn_tmp);
	    SFr.insert(fn_tmp);
	    DF.insert_data_line(mysampling.no_redundant_sampling_points_angles[i]);
	}
	if (output_classes)
	{
	    empty.rot()=rot;
	    empty.tilt()=tilt;
	    class_avgs.push_back(empty);
	    class_selfiles.push_back(emptySF);
	}
	compute_stats_within_binary_mask(rotmask,proj(),dummy,dummy,mean_ref,stddev_ref);
	proj()-=mean_ref;
	ref_img.push_back(proj());
	ref_stddev[i]=stddev_ref;
	ref_mean[i]=mean_ref;
	if (verb>0 && (i%MAX(1,nl/60)==0)) progress_bar(i);
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


void Prog_projection_matching_prm::PM_process_one_image(matrix2D<double> &Mexp,
							float &img_rot, float &img_tilt, float &img_psi,
							int &opt_dirno, double &opt_psi,
							double &opt_xoff, double &opt_yoff,
							double &maxCC, double &Zscore) {


  // Rotational search ====================================================
  matrix2D<double> Mimg,Mref,Maux,Mcorr;
  double psi,psi_min,psi_max,psi_ref,thisCC,oldCC,aveCC=0.,varCC=0.;
  double stddev_img,mean_img,dummy,xmax,ymax;
  int dirno,c=0,ioptpsi=0,ioptflip=0;

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

  // INSERT THE FOLLOWING FUNCTION:
  //int my_i=mysampling.get_my_nearest_sampling_point(img_rot,img_tilt,img_psi)
  int my_i=0;

  // Calculate correlation coefficients for all angles
  psi_min=img_psi-MIN(180.,ang_search);
  psi_max=img_psi+MIN(180.,ang_search);
  for (psi=psi_min; psi<psi_max; psi+=sampling)
  {
    Mimg=Maux.rotate(psi,DONT_WRAP);
    for (int j = 0; j < mysampling.my_neighbors[my_i].size();j++)
    {
      dirno=mysampling.my_neighbors[my_i][j];
      Mref=ref_img[dirno];
      psi_ref=mysampling.my_neighbors_psi[my_i][j];
      if (psi_ref!=0.) 
      { 
          // or -psi_ref, or psi_ref + pi divido por e
          Mimg=Mimg.rotate(psi_ref,DONT_WRAP);
      }
      thisCC=0.;
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mimg) {
	thisCC+=dMij(Mref,i,j)*dMij(Mimg,i,j);
      }
      thisCC/=ref_stddev[dirno]*stddev_img*nr_pixels_rotmask;
      c++;
      oldCC=aveCC;
      aveCC+=(thisCC-oldCC)/(c+1);
      if (c>1) varCC=(1.-1./(double)c)*varCC+(c+1.)*(oldCC-aveCC)*(oldCC-aveCC);
      if (thisCC>maxCC) {
	maxCC=thisCC;
	opt_psi=psi;
	opt_dirno=dirno;
      }
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
  double opt_rot,opt_tilt,opt_psi,opt_xoff,opt_yoff,maxCC,Zscore;
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
    opt_rot=XX(mysampling.no_redundant_sampling_points_angles[opt_dirno]);
    opt_tilt=YY(mysampling.no_redundant_sampling_points_angles[opt_dirno]);

    sumCC+=maxCC;
    dataline(0)=opt_rot;                 // rot
    dataline(1)=opt_tilt;                // tilt
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
	img.set_eulerAngles(opt_rot,opt_tilt,opt_psi);
	img.set_originOffsets(opt_xoff,opt_yoff);
	img.write(fn_img);
    }
    if (output_classes)
    {
	// Re-read image to get the untransformed image matrix again
	img.read(fn_img);
	img.set_eulerAngles(opt_rot,opt_tilt,opt_psi);
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
  ref_img.clear();
}
void Prog_projection_matching_prm::write_classes()
{

    FileName fn_base,fn_img,fn_sel;
    SelFile SF,SF2;

    fn_base=fn_root+"_class";
    SF.clear();
    SF2.clear();
    for (int i=0; i<mysampling.no_redundant_sampling_points_vector.size(); i++)
    {
	fn_img.compose(fn_base,i+1,"xmp");
	SF.insert(fn_img);
	fn_sel.compose(fn_base,i+1,"sel");
	class_avgs[i]()/=class_avgs[i].weight();
	class_avgs[i].write(fn_img);
	class_selfiles[i].write(fn_sel);
    }
    fn_base+="es.sel";
    SF.write(fn_base);
}
