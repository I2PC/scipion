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
void Prog_projection_matching_prm::read(int argc, char **argv) _THROW  {

  // Read command line
  if (check_param(argc,argv,"-show_all_options")) { usage(); extended_usage();}
  fn_vol=get_param(argc,argv,"-vol");
  SF.read(get_param(argc,argv,"-i"));
  SF.ImgSize(dim,dim);
  fn_out=get_param(argc,argv,"-doc","out.doc");
  sampling=AtoF(get_param(argc,argv,"-sam","10"));
  max_shift=AtoF(get_param(argc,argv,"-max_shift","5"));

  // Additional commands
  rot_range=AtoF(get_param(argc,argv,"-rot_search","-1"));
  tilt_range=AtoF(get_param(argc,argv,"-tilt_search","-1"));
  psi_range=AtoF(get_param(argc,argv,"-psi_search","-1"));
  fn_sym=get_param(argc,argv,"-sym","");
  output_refs=check_param(argc,argv,"-output_refs");
  fn_refs=get_param(argc,argv,"-ref_name","ref");
  modify_header=!check_param(argc,argv,"-dont_modify_header");

  // Hidden stuff
  verb=AtoI(get_param(argc,argv,"-verb","1"));

  if ( (fn_sym!="") && (rot_range>=0||tilt_range>=0)) 
    REPORT_ERROR(1," Limited rot & tilt_search not implemented for symmetry!");

}

// Show ====================================================================
void Prog_projection_matching_prm::show() {

  if (verb>0) {
    cerr << "  Input images            : "<< SF.name()<<" ("<<SF.ImgNo()<<")"<<endl;
    cerr << "  Reference volume        : "<< fn_vol<<endl;
    cerr << "  Output docfile          : "<< fn_out<<endl;
    cerr << "  Angular sampling rate   : "<< sampling <<endl;
    if (max_shift>0) {
      cerr << "  -> Limit search of origin offsets to  +/- "<<max_shift<<" pixels"<<endl;
    }
    if (rot_range>0) {
      cerr << "  -> Limit search of rot-angle to  +/- "<<rot_range<<" degrees"<<endl;
    }
    if (tilt_range>0) {
      cerr << "  -> Limit search of tilt-angle to +/- "<<tilt_range<<" degrees"<<endl;
    }
    if (rot_range>0) {
      cerr << "  -> Limit search of psi-angle to  +/- "<<psi_range<<" degrees"<<endl;
    }
    if (fn_sym!="") {
      cerr << "  -> Limit angular search to asymmetric part, as defined by: "<<fn_sym<<endl;
    }
    if (!modify_header) {
      cerr << "  -> Do not modify the image headers (only output docfile)"<<endl;
    }
    if (output_refs) {
      cerr << "  -> Output library projections, sel and docfile with rootname: "<<fn_refs<<endl;
    }

    cerr << " ================================================================="<<endl;
  }

} 

// Usage ===================================================================
void Prog_projection_matching_prm::usage() {
  cerr << "Usage:  projection_matching [options] "<<endl;
  cerr << "   -i <selfile>                : Selfile with input images \n"
       << "   -vol <volume>               : Reference volume \n"
       << " [ -doc <docfile> ]            : Output document file (default= out.doc) \n"
       << " [ -sam <float=10> ]           : Sampling rate for rot, tilt & psi (degrees) \n"
       << " [ -max_shift <float=5> ]      : Maximum change in origin offset (+/- pixels) \n"
       << " [ -show_all_options ]         : Show all program options\n";
}

// Extended usage ===================================================================
void Prog_projection_matching_prm::extended_usage() {
  cerr << "Additional options: \n"
       << " [ -rot_search <float=-1> ]    : Maximum change in rot  (+/- degrees) \n"
       << " [ -tilt_search <float=-1> ]   : Maximum change in tilt (+/- degrees) \n"
       << " [ -psi_search <float=-1> ]    : Maximum change in psi  (+/- degrees) \n"
       << " [ -sym <symfile> ]            : Limit angular search to asymmetric part \n"
       << " [ -output_refs  ]             : Output reference projections, sel and docfile \n"
       << " [ -ref_name <root=\"ref\"> ]        => corresponding rootname \n"
       << " [ -dont_modify_header ]       : Do not store alignment parameters in the image headers \n";
  exit(1);
}

/* Check whether projection directions are unique ---------------------- */
bool Prog_projection_matching_prm::directions_are_unique(double rot,  double tilt,
							 double rot2, double tilt2, 
							 double rot_limit, double tilt_limit,
							 SymList &SL) {
  
  bool are_unique=true;
  double rot2p,tilt2p,psi2p, psi2=0.;
  double diff_rot,diff_tilt;
  matrix2D<double>  L(4,4), R(4,4);
  
  for (int isym=0; isym<=SL.SymsNo(); isym++) {
    
    if (isym==0) {rot2p=rot2; tilt2p=tilt2; psi2p=psi2;}
    else {
      SL.get_matrices(isym-1,L,R);
      L.resize(3,3); // Erase last row and column
      R.resize(3,3); // as only the relative orientation
      // is useful and not the translation
      Euler_apply_transf(L,R,rot2,tilt2,psi2,rot2p,tilt2p,psi2p);
    }

    diff_rot=rot-rot2p;
    diff_tilt=tilt-tilt2p;
    diff_rot=ABS(realWRAP(diff_rot,-180,180));
    diff_tilt=ABS(realWRAP(diff_tilt,-180,180));
    if ((rot_limit-diff_rot)>1e-3 && (tilt_limit-diff_tilt)>1e-3) are_unique=false;
    Euler_another_set(rot2p,tilt2p,psi2p,rot2p,tilt2p,psi2p);
    diff_rot=rot-rot2p;
    diff_tilt=tilt-tilt2p;
    diff_rot=ABS(realWRAP(diff_rot,-180,180));
    diff_tilt=ABS(realWRAP(diff_tilt,-180,180));
    if ((rot_limit-diff_rot)>1e-3 && (tilt_limit-diff_tilt)>1e-3) are_unique=false;
  }

  return are_unique;

}

/* Fill DF with evenly distributed rot & tilt  ----------------------------- */
void Prog_projection_matching_prm::make_even_distribution(DocFile &DF, double &sampling, 
							  SymList &SL, bool exclude_mirror) {

  int rot_nstep,tilt_nstep=ROUND(180./sampling)+1;
  double rotp,tiltp,psip,rot_sam,tilt,rot,tilt_sam,psi=0.;
  bool append;
  matrix1D<double> dataline(3);
  tilt_sam=(180./tilt_nstep);

  DF.clear();
  // Create evenly distributed angles
  for (int tilt_step=0; tilt_step<tilt_nstep; tilt_step++) {
    tilt=((double)tilt_step/(tilt_nstep-1))*180.;
    if (tilt>0) rot_nstep=CEIL(360.*sin(DEG2RAD(tilt))/sampling);
    else rot_nstep=1;
    rot_sam=360./(double)rot_nstep;
    for (double rot=0.; rot<360.; rot+=rot_sam) {
      // Check whether by symmetry or mirror the angle has been included already
      append=true;
      DF.go_first_data_line();
      while (!DF.eof()) {
	if (!directions_are_unique(rot,tilt,DF(0),DF(1),rot_sam,tilt_sam,SL)) 
	  append=false;
	if (exclude_mirror) {
	  Euler_up_down(rot,tilt,psi,rotp,tiltp,psip);
	  if (!directions_are_unique(rotp,tiltp,DF(0),DF(1),rot_sam,tilt_sam,SL)) 
	  append=false;
	}
	DF.next_data_line();
      }
      if (append) {
	dataline(0)=rot;
	dataline(1)=tilt;
	dataline(2)=0.;
	DF.append_data_line(dataline);
      }
    }
  }
  DF.go_beginning();

}

// Side info stuff ===================================================================
void Prog_projection_matching_prm::produce_Side_info() _THROW {

  VolumeXmipp     vol;
  Projection      proj;
  matrix2D<double> A(3,3);
  DocFile         DF;
  SelFile         SF;
  SymList         SL;
  FileName        fn_tmp;

  // Set nr_psi
  nr_psi=CEIL(360./sampling);
  sampling=360./nr_psi;
 
  // Create max_shift mask
  shiftmask.resize(dim,dim);
  shiftmask.set_Xmipp_origin();
  if (max_shift<0.) max_shift=(double)dim/2.;
  BinaryCircularMask(shiftmask,max_shift,INNER_MASK);

  // Create evenly-distributed reference projection angles
  DF.clear();
  if (fn_sym!="") SL.read_sym_file(fn_sym);
  make_even_distribution(DF,sampling,SL,false);

  // Create reference projection images
  double          mean_ref,stddev_ref,dummy,psi=0.;
  int             nl;

  ref_img.clear();
  ref_rot.clear();
  ref_tilt.clear();
  ref_mean.clear();
  ref_stddev.clear();
  vol.read(fn_vol);
  vol().set_Xmipp_origin();

  nl=DF.dataLineNo();
  SF.reserve(nl);
  SF.go_beginning();
  DF.go_beginning();
  if (verb>0) cerr << "--> Projecting the reference volume ..."<<endl;
  if (verb>0) init_progress_bar(nl);

  DF.adjust_to_data_line();
  nr_dir=0;
  while (!DF.eof()) {
    ref_rot.push_back(DF(0));
    ref_tilt.push_back(DF(1));
    project_Volume(vol(),proj,dim,dim,ref_rot[nr_dir],ref_tilt[nr_dir],psi);
    if (output_refs) {
      fn_tmp.compose(fn_refs,nr_dir+1,"xmp");
      proj.write(fn_tmp);
      SF.insert(fn_tmp);
    }
    proj().compute_stats(mean_ref,stddev_ref,dummy,dummy);
    proj()-=mean_ref;
    ref_img.push_back(proj());
    ref_stddev.push_back(stddev_ref);
    ref_mean.push_back(mean_ref);
    DF.next_data_line();
    nr_dir++;
    if (verb>0 && (nr_dir%MAX(1,nl/60)==0)) progress_bar(nr_dir);
  }
  if (output_refs) {
    fn_tmp=fn_refs+".doc";
    DF.write(fn_tmp);
    fn_tmp=fn_refs+".sel";
    SF.write(fn_tmp);
  }

  if (verb>0) progress_bar(nl);
  if (verb>0) cerr << " ================================================================="<<endl;

}


void Prog_projection_matching_prm::PM_process_one_image(matrix2D<double> &Mexp,
							float &img_rot, float &img_tilt, float &img_psi, 
							int &opt_dirno, double &opt_psi,
							double &opt_xoff, double &opt_yoff, 
							double &maxCC, double &Z) _THROW {


  // Rotational search ====================================================
  matrix2D<double> Mimg,Mref,Mcorr;
  double act_rot_range,psi,thisCC,oldCC,aveCC=0.,varCC=0.;
  double mean_ref,stddev_ref,stddev_img,mean_img,dummy;
  int c=0,ioptpsi=0,ioptflip=0;
  bool search;

  maxCC=-99.e99; 
  Mimg.resize(dim,dim);
  Mimg.set_Xmipp_origin();
  Mref.resize(dim,dim);
  Mref.set_Xmipp_origin();

  // Calculate correlations for all angles
  FOR_ALL_ROTATIONS() {
    psi=(double)(ipsi*360./nr_psi);
    Mimg=Mexp.rotate(psi,DONT_WRAP);
    Mimg.compute_stats(mean_ref,stddev_img,dummy,dummy);
    Mimg-=mean_img;
    FOR_ALL_DIRECTIONS() {
      search=true;
      Mref=ref_img[dirno];
      if (rot_range>0) {
	// Rot_range is tilt-angle dependent!
	if (ref_tilt[dirno]>0 && ref_tilt[dirno]<180) act_rot_range=rot_range/sin(DEG2RAD(ref_tilt[dirno])); 
	else act_rot_range=361.;
	if (ABS(realWRAP(img_rot-ref_rot[dirno],-180.,180.)) > act_rot_range) search=false;
      }
      if (search && tilt_range>0)
	if (ABS(realWRAP(img_tilt-ref_tilt[dirno],-180.,180.)) > tilt_range) search=false;
      if (search && psi_range>0) 
	if (ABS(realWRAP(img_psi-psi,-180.,180.)) > psi_range) search=false;
      if (search) {
	thisCC=0.;
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mimg) {
	  thisCC+=dMij(Mref,i,j)*dMij(Mimg,i,j);
	}
	thisCC/=ref_stddev[dirno]*stddev_img*dim*dim;
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
  }

  // Calculate Z-score on rotational permutations
  Z=(maxCC-aveCC)/(sqrt(varCC));

  // Interpolated translational search =======================================================
  int              imax,jmax,i_actual,j_actual;
  double           max,xmax,ymax,sumcorr,avecorr,stdcorr;
  float            xshift,yshift,shift;
  int              n_max=-1;
  bool             neighbourhood=TRUE;

  // Calculate cross-correlation matrix to search shifts
  Mimg=Mexp.rotate(opt_psi,DONT_WRAP);
  Mref=ref_img[opt_dirno];
  Mref+=ref_mean[opt_dirno];
  correlation_matrix(Mimg,Mref,Mcorr);
  // Adjust statistics within shiftmask to average 0 and stddev 1
  compute_stats_within_binary_mask(shiftmask,Mcorr,dummy,dummy,avecorr,stdcorr);
  FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr) {
    if (MAT_ELEM(shiftmask,i,j)) 
      MAT_ELEM(Mcorr,i,j)=(MAT_ELEM(Mcorr,i,j)-avecorr)/stdcorr;
    else MAT_ELEM(Mcorr,i,j)=0.;
  }
  Mcorr.max_index(imax,jmax);
  max=MAT_ELEM(Mcorr,imax,jmax);

  while (neighbourhood) {
    n_max ++;
    for (int i=-n_max; i <= n_max; i++)
      for (int j=-n_max; j <= n_max; j++) {
        i_actual = i+imax;
        j_actual = j+jmax;
        if (i_actual < Mcorr.startingY()  || j_actual < Mcorr.startingX() &&
            i_actual > Mcorr.finishingY() || j_actual > Mcorr.finishingX() )
          neighbourhood=FALSE;
        else if (max/1.414 > MAT_ELEM(Mcorr,i_actual,j_actual))
            neighbourhood=FALSE;
      }
  }

  // We have the neighbourhood => looking for the gravity centre
  xmax = ymax = sumcorr = 0.;
  for (int i=-n_max; i <= n_max; i++)
    for (int j=-n_max; j <= n_max; j++) {
      i_actual = i+imax;
      j_actual = j+jmax;
      if (i_actual >= Mcorr.startingY()  && j_actual >= Mcorr.startingX() &&
          i_actual <= Mcorr.finishingY() && j_actual <= Mcorr.finishingX() ) {
        ymax += i_actual*MAT_ELEM(Mcorr,i_actual,j_actual);
        xmax += j_actual*MAT_ELEM(Mcorr,i_actual,j_actual);
        sumcorr += MAT_ELEM(Mcorr,i_actual,j_actual);
      }
    }
  xmax/= sumcorr; ymax/= sumcorr;

  opt_xoff=-xmax*COSD(opt_psi)-ymax*SIND(opt_psi);
  opt_yoff=xmax*SIND(opt_psi)-ymax*COSD(opt_psi);

  // Calculate correlation coefficient
  Mimg=Mimg.translate(vector_R2(-xmax,-ymax)); 
  maxCC=correlation_index(Mimg,Mref);

}

void Prog_projection_matching_prm::PM_loop_over_all_images(SelFile &SF, DocFile &DFo, double &sumCC, double &sumZ) _THROW {


  ImageXmipp img;
  FileName fn_img;
  matrix1D<double> dataline(7);  
  double opt_psi,opt_xoff,opt_yoff,maxCC,Z;
  int c,nn,imgno,opt_dirno;

  if (verb>0) cerr << "--> Projection matching ... "<<endl;

  // Initialize
  nn=SF.ImgNo();
  if (verb>0) init_progress_bar(nn);
  c=MAX(1,nn/60);
  
  // Loop over all images
  sumCC=0.;
  sumZ=0.;
  imgno=0;
  SF.go_beginning();
  while ((!SF.eof())) {
    fn_img=SF.NextImg();
    img.read(fn_img,FALSE,FALSE,FALSE,TRUE);
    img().set_Xmipp_origin();

    // Perform the projection matching for each image separately
    PM_process_one_image(img(),img.Phi(),img.Theta(),img.Psi(),opt_dirno,opt_psi,opt_xoff,opt_yoff,maxCC,Z);

    opt_xoff+=img.Xoff();
    opt_yoff+=img.Yoff();

    sumCC+=maxCC;
    sumZ+=Z;
    dataline(0)=ref_rot[opt_dirno];      // rot
    dataline(1)=ref_tilt[opt_dirno];     // tilt
    dataline(2)=opt_psi;                 // psi
    dataline(3)=opt_xoff;                // Xoff
    dataline(4)=opt_yoff;                // Yoff
    dataline(5)=maxCC;                   // maximum CC
    dataline(6)=Z;                       // Z-score
    DFo.append_comment(img.name());
    DFo.append_data_line(dataline);

    if (modify_header) {
      // Re-read image to get the untransformed image matrix again
      img.read(fn_img);
      img.set_eulerAngles(ref_rot[opt_dirno],ref_tilt[opt_dirno],opt_psi);
      img.set_originOffsets(opt_xoff,opt_yoff);
      img.write(fn_img);
    }

    if (verb>0) if (imgno%c==0) progress_bar(imgno);
    imgno++;
  }

  if (verb>0) progress_bar(nn);
  if (verb>0) cerr << " ================================================================="<<endl;

}
