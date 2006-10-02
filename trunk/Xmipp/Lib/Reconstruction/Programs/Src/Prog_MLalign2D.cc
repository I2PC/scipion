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
#include "../Prog_MLalign2D.hh"

// Read arguments ==========================================================
void Prog_MLalign2D_prm::read(int argc, char ** argv, bool ML3D)  {


  // Generate new command line for restart procedure
  cline="";
  if (check_param(argc,argv,"-restart")) {
    string comment;
    FileName fn_sel;
    DocFile DFi;
    DFi.read(get_param(argc,argv,"-restart"));
    DFi.go_beginning();
    comment=(DFi.get_current_line()).get_text();
    if (strstr(comment.c_str(),"MLalign2D-logfile")==NULL) {
      cerr << "Error!! Docfile is not of MLalign2D-logfile type. "<<endl;
      exit(1);
    } else {
      char *copy;
      int n=0;
      int nmax=DFi.dataLineNo();
      SFr.reserve(nmax);
      copy=NULL; 
      DFi.next();
      comment=" -frac "+DFi.name();
      if (!ML3D) {
	fn_sel=DFi.name();
	fn_sel=fn_sel.without_extension()+"_restart.sel";
	comment+=" -ref "+fn_sel;
      } 
      comment+=(DFi.get_current_line()).get_text();
      DFi.next();
      cline=(DFi.get_current_line()).get_text();
      comment=comment+cline;
      generate_command_line(comment,argc,argv,copy);
      if (!ML3D) {
	// Read images names from restart file
	DFi.next();
	while (n<nmax) {
	  n++;
	  DFi.next();
	  if (DFi.get_current_line().Is_comment()) fn_sel=((DFi.get_current_line()).get_text()).erase(0,3);
	  SFr.insert(fn_sel,SelLine::ACTIVE);
	  DFi.adjust_to_data_line();
	}
	fn_sel=DFi.name();
	fn_sel=fn_sel.without_extension()+"_restart.sel";
	SFr.write(fn_sel);
	SFr.clear();
      }
    }
  } else {
    for (int i=1; i<argc; i++) {
      cline=cline+(string)argv[i]+" ";
    }
  }
    
  // Read command line
  if (check_param(argc,argv,"-more_options")) { usage(); extended_usage();}
  n_ref=AtoI(get_param(argc,argv,"-nref","0"));
  fn_ref=get_param(argc,argv,"-ref","");
  SF.read(get_param(argc,argv,"-i"));
  fn_root=get_param(argc,argv,"-o","MLalign2D");
  psi_step=AtoF(get_param(argc,argv,"-psi_step","5"));
  Niter=AtoI(get_param(argc,argv,"-iter","100"));
  istart=AtoI(get_param(argc,argv,"-istart","1"));
  sigma_noise=AtoF(get_param(argc,argv,"-noise","1"));
  sigma_offset=AtoF(get_param(argc,argv,"-offset","3"));
  do_mirror=check_param(argc,argv,"-mirror");
  eps=AtoF(get_param(argc,argv,"-eps","5e-5"));
  fn_frac=get_param(argc,argv,"-frac","");
  write_docfile=!check_param(argc,argv,"-dont_output_docfile");
  write_selfiles=!check_param(argc,argv,"-dont_output_selfiles");
  write_intermediate=!check_param(argc,argv,"-dont_output_intermediate");
  fix_fractions=check_param(argc,argv,"-fix_fractions");
  fix_sigma_offset=check_param(argc,argv,"-fix_sigma_offset");
  fix_sigma_noise=check_param(argc,argv,"-fix_sigma_noise");
  verb=AtoI(get_param(argc,argv,"-verb","1"));
  maxCC_rather_than_ML=check_param(argc,argv,"-maxCC");
  fast_mode=check_param(argc,argv,"-fast");
  C_fast=AtoF(get_param(argc,argv,"-C","1e-12"));
  max_shift=AtoF(get_param(argc,argv,"-max_shift","-1"));
  save_mem1=check_param(argc,argv,"-save_memA");
  save_mem2=check_param(argc,argv,"-save_memB");
  save_mem3=check_param(argc,argv,"-save_memC");
  search_shift=AtoF(get_param(argc,argv,"-search_shift","999."));
  fn_doc=get_param(argc,argv,"-doc","");

  // Hidden arguments
  do_esthetics=check_param(argc,argv,"-esthetics");

  //only for interaction with MLrefine3D:
  search_rot=AtoF(get_param(argc,argv,"-search_rot","999."));
  zero_offsets=check_param(argc,argv,"-zero_offsets");

}

// Show ====================================================================
void Prog_MLalign2D_prm::show(bool ML3D) {

  if (verb>0) {
    // To screen
    if (!ML3D) {
    cerr << " -----------------------------------------------------------------"<<endl;
    cerr << " | Read more about this program in the following publications:   |"<<endl;
    cerr << " |  Scheres ea. (2005) J.Mol.Biol. 348(1), 139-49                |"<<endl;
    cerr << " |  Scheres ea. (2005) Bioinform. 21(suppl.2), ii243-4   (-fast) |"<<endl;
    cerr << " |                                                               |"<<endl;
    cerr << " |  *** Please cite them if this program is of use to you! ***   |"<<endl;
    cerr << " -----------------------------------------------------------------"<<endl;
    }
    cerr << "--> Maximum-likelihood multi-reference refinement "<<endl; 
    cerr << "  Input images            : "<< SF.name()<<" ("<<SF.ImgNo()<<")"<<endl;
    if (fn_ref!="")
      cerr << "  Reference image(s)      : "<< fn_ref<<endl;
    else
      cerr << "  Number of references:   : "<< n_ref<<endl;
    cerr << "  Output rootname         : "<< fn_root<<endl;
    cerr << "  Stopping criterium      : "<< eps<<endl;
    cerr << "  initial sigma noise     : "<< sigma_noise<<endl;
    cerr << "  initial sigma offset    : "<< sigma_offset<<endl;
    cerr << "  Psi sampling interval   : "<< psi_step<<endl;
    if (do_mirror) 
    cerr << "  Check mirrors           : true"<<endl;
    else 
    cerr << "  Check mirrors           : false"<<endl;
    if (fn_frac!="") 
    cerr << "  Initial model fractions : "<< fn_frac<<endl;
    if (maxCC_rather_than_ML) {
      cerr << "  -> Use a maxCC instead of a maximum likelihood target."<<endl;
    }
    if (fast_mode) {
      cerr << "  -> Use fast, reduced search-space approach with C = "<<C_fast<<endl;
      if (zero_offsets)
	cerr << "    + Start from all-zero translations"<<endl;
    }
    if (search_shift<999.)
      cerr << "    + Limit translational search to +/- "<<search_shift<<" pixels"<<endl;
    if (search_rot<180.)
      cerr << "    + Limit orientational search to +/- "<<search_rot<<" degrees"<<endl;
    if (save_mem1) 
      cerr << "  -> Save_memory A: recalculate real-space rotations in -fast"<<endl;
    if (save_mem2) 
      cerr << "  -> Save_memory B: limit translations to 3 sigma_offset "<<endl;   
    if (save_mem3)
      cerr << "  -> Save_memory C: do not store rotated references; rotate experimental image instead "<<endl;


    // Hidden stuff
    if (!write_intermediate) {
      cerr << "  -> Do not write out images after each iteration."<<endl;
    }
    if (fix_fractions && !maxCC_rather_than_ML) {
      cerr << "  -> Do not update estimates of model fractions."<<endl;
    }
    if (fix_sigma_offset && !maxCC_rather_than_ML) {
      cerr << "  -> Do not update sigma-estimate of origin offsets."<<endl;
    }
    if (fix_sigma_noise && !maxCC_rather_than_ML) {
      cerr << "  -> Do not update sigma-estimate of noise."<<endl;
    }
    if (do_esthetics) {
      cerr << "  -> Perform esthetics on (0,0)-pixel artifacts"<<endl;
    }
    cerr << " -----------------------------------------------------------------"<<endl;

  }

} 

// Usage ===================================================================
void Prog_MLalign2D_prm::usage() {
  cerr << "Usage:  MLalign2D [options] "<<endl;
  cerr << "   -i <selfile>                : Selfile with input images \n";
  cerr << "   -ref <selfile/image>        : Selfile with initial references/single reference image \n";
  cerr << "      OR -nref <int>               OR number of references to generate automatically (bias-free)\n";
  cerr << " [ -o <rootname> ]             : Output rootname (default = \"MLalign2D\")\n";
  cerr << " [ -mirror ]                   : Also check mirror image of each reference \n";
  cerr << " [ -fast ]                     : Use pre-centered images to pre-calculate significant orientations\n";
  cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_MLalign2D_prm::extended_usage(bool ML3D) {
  cerr << "Additional MLalign2D-options: "<<endl;
  cerr << " [ -eps <float=5e-5> ]         : Stopping criterium \n";
  cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
  cerr << " [ -psi_step <float=5> ]       : In-plane rotation sampling interval [deg]\n";
  cerr << " [ -noise <float=1> ]          : Expected standard deviation for pixel noise \n";
  cerr << " [ -offset <float=3> ]         : Expected standard deviation for origin offset [pix]\n";
  cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
  cerr << " [ -C <double=1e-12> ]         : Significance criterion for fast approach \n";
  if (!ML3D) cerr << " [ -restart <logfile> ]        : restart a run with all parameters as in the logfile \n";
  if (!ML3D) cerr << " [ -istart <int> ]             : number of initial iteration \n";
  cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
  cerr << " [ -fix_sigma_offset]          : Do not re-estimate the standard deviation in the origin offsets \n";
  cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
  cerr << " [ -dont_output_docfile ]      : Do not write out docfile with most likely angles & translations \n";
  cerr << " [ -dont_output_selfiles ]     : Do not write out selfiles with most likely class assignments \n";
  cerr << " [ -maxCC ]                    : Use maximum cross-correlation instead of maximum likelihood target \n";
  cerr << " [ -search_shift <float=999>]  : Limit of translational search [pix] (does NOT use FFT) \n";
  cerr << " [ -max_shift <float=dim/4>]   : Dont trust shifts larger than max_shift \n";
  cerr << " [ -doc <docfile=\"\"> ]         : Read initial angles and offsets from docfile \n";
  cerr << endl;
  exit(1);
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void Prog_MLalign2D_prm::produce_Side_info() {

  FileName                    fn_img,fn_tmp,fn_base,fn_tmp2;
  ImageXmipp                  img;
  FourierImageXmipp           fourimg;
  SelLine                     SL;
  SelFile                     SFtmp;
  matrix1D<double>            offsets(2);
  matrix2D<double>            A(3,3),Maux,Maux2;
  vector<int>                 tmppointp,tmppointp_nolow,tmppointi,tmppointj;
  bool                        uniqname,nomoredirs;
  float                       xx,yy;
  double                      av,psi;
  int                         im,jm;

  // Get image size
  SF.ImgSize(dim,dim);
  hdim=dim/2;
  dim2=dim*dim;

  // Get total number of images
  nr_exp_images=SF.ImgNo();

  // Get number of references
  if (fn_ref!="") {
    if (Is_ImageXmipp(fn_ref)) n_ref=1;
    else  {
      SFr.read(fn_ref);
      n_ref=SFr.ImgNo();
    }
  }

  // Check the uniqueness of all filenames (for names of temporary offsets files)
  uniqname=false;
  nomoredirs=false;
  offsets_keepdir=0;
  while ( (uniqname==false) && (nomoredirs==false) ) {
    SF.go_beginning();
    SFtmp.clear();
    while (!SF.eof()) {
      fn_tmp=SF.NextImg();
      fn_tmp2=fn_tmp.remove_directories(offsets_keepdir);
      if (fn_tmp==fn_tmp2) nomoredirs=true;
      SFtmp.insert(fn_tmp2);
    }
    SFtmp.sort_by_filenames();
    SFtmp.go_beginning();
    uniqname=true;
    while (!SFtmp.eof()) {
      fn_tmp=SFtmp.NextImg();
      fn_tmp2=SFtmp.NextImg();
      if (fn_tmp==fn_tmp2) {
	uniqname=false;
	offsets_keepdir++;
	break;
      }
    }
  }
  SFtmp.clear();
  if (!uniqname) 
    REPORT_ERROR(1,"Prog_MLalign2D_prm: Provide a selfile with unique image names (preferably all in one directory)");

  // Set nr_psi & nr_flip and construct flipping matrices
  if (save_mem3) {
    nr_psi=1;
    nr_flip=nr_nomirror_flips=CEIL(360./psi_step);
    psi_step=360./nr_psi;
    // store all rotation (and mirror) matrices
    FOR_ALL_FLIPS() {
      double ang =(double)(iflip*360./nr_flip)+SMALLANGLE;
      A=rot2D_matrix(ang);
      F.push_back(A);
    }
    if (do_mirror) {
      FOR_ALL_FLIPS() {
        double ang =(double)(iflip*360./nr_flip);
        A=rot2D_matrix(ang);
        A(0,0)*=-1.;
        A(0,1)*=-1.;
        F.push_back(A);
      }
      nr_flip*=2;
    }
  } else {
    psi_max=90.;
    nr_psi=CEIL(psi_max/psi_step);
    psi_step=psi_max/nr_psi; 
    nr_flip=nr_nomirror_flips=4;
    // 0, 90, 180 & 270 degree flipping, as well as mirror 
    A.init_identity();
    F.push_back(A);
    A(0,0)=0.; A(1,1)=0.; A(1,0)=1.; A(0,1)=-1;
    F.push_back(A);
    A(0,0)=-1.; A(1,1)=-1.; A(1,0)=0.; A(0,1)=0;
    F.push_back(A);
    A(0,0)=0.; A(1,1)=0.; A(1,0)=-1.; A(0,1)=1;
    F.push_back(A);
    if (do_mirror) {
      nr_flip=8;
      A.init_identity(); A(0,0)=-1;
      F.push_back(A);
      A(0,0)=0.; A(1,1)=0.; A(1,0)=1.; A(0,1)=1;
      F.push_back(A);
      A(0,0)=1.; A(1,1)=-1.; A(1,0)=0.; A(0,1)=0;
      F.push_back(A);
      A(0,0)=0.; A(1,1)=0.; A(1,0)=-1.; A(0,1)=-1;
      F.push_back(A);
    }
  }

  // Set some stuff for maxCC-mode
  if (maxCC_rather_than_ML) {
    fix_sigma_noise=true;
    fix_sigma_offset=true;
    fix_fractions=true;
  }

  // Set & check max_shift
  if (max_shift<0) max_shift=dim/4.;

  // Set limit_rot & limit_trans
  if (search_rot<180.) limit_rot=true;
  else limit_rot=false;
  if (search_shift<999) limit_trans=true;
  else limit_trans=false;

  // Fill limited translation shift-vectors
  if (limit_trans) {
    nr_trans=0;
    Maux.resize(dim,dim);
    Maux.set_Xmipp_origin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Maux) {
      double rr=(double)i*i+(double)j*j;
      if (rr<=search_shift*search_shift) {
	offsets(0)=(double)j;
	offsets(1)=(double)i;
	Vtrans.push_back(offsets);
	if (i==0 && j==0) zero_trans=nr_trans;
	nr_trans++;
      }
    }
  }


}

// Generate initial references =============================================
void Prog_MLalign2D_prm::generate_initial_references()  {

  SelFile SFtmp, SFout;
  ImageXmipp Iave,Itmp;
  double dummy;
  FileName fn_tmp;

  // Make random subsets and calculate average images
  if (verb>0) {
    cerr << "  Generating initial references by averaging over random subsets" <<endl;
    init_progress_bar(n_ref);
  }
  SFtmp=SF.randomize();
  int Nsub=ROUND((double)SFtmp.ImgNo()/n_ref);
  for (int refno=0; refno<n_ref; refno++) {
    SFout.clear();
    SFout.reserve(Nsub);
    SFtmp.go_beginning();
    SFtmp.jump_lines(Nsub*refno);
    if (refno==n_ref-1) Nsub=SFtmp.ImgNo()-refno*Nsub;
    for (int nn=0; nn<Nsub; nn++) {
      SFout.insert(SFtmp.current());
      SFtmp.NextImg();
    }
    SFout.get_statistics(Iave,Itmp,dummy,dummy);
    fn_tmp=fn_root+"_it";
    fn_tmp.compose(fn_tmp,0,"");
    fn_tmp=fn_tmp+"_ref";
    fn_tmp.compose(fn_tmp,refno+1,"");
    fn_tmp=fn_tmp+".xmp";
    Iave.write(fn_tmp);
    SFr.insert(fn_tmp,SelLine::ACTIVE);
    if (verb>0) progress_bar(refno);
  }
  if (verb>0) progress_bar(n_ref);
  fn_ref=fn_root+"_it";
  fn_ref.compose(fn_ref,0,"sel");
  SFr.write(fn_ref);

}

// Read reference images to memory and initialize offset vectors
// This side info is NOT general, i.e. in parallel mode it is NOT the
// same for all processors! (in contrast to produce_Side_info)
void Prog_MLalign2D_prm::produce_Side_info2() {

  int                       c,refno=0;
  DocFile                   DF,DF2;
  DocLine                   DL;
  double                    offx,offy,sumfrac=0.;
  FileName                  fn_tmp;
  ImageXmipp                img;
  matrix1D<double>          offsets(2);

  // Read in all reference images in memory
  if (Is_ImageXmipp(fn_ref)) {
    SFr.reserve(1);
    SFr.insert(fn_ref);
  } else {
    SFr.read(fn_ref);
  }
  n_ref=0;
  SFr.go_beginning();
  while ((!SFr.eof())) {
    img.read(SFr.NextImg(),false,false,true,false);
    img().set_Xmipp_origin();
    Iref.push_back(img);
    if (!save_mem3) Iold.push_back(img);
    // Default start is all equal model fractions
    alpha_k.push_back((double)1/SFr.ImgNo());
    Iref[refno].weight()=alpha_k[refno]*(double)nr_exp_images;
    // Default start is half-half mirrored images
    if (do_mirror) mirror_fraction.push_back(0.5);
    else mirror_fraction.push_back(0.);
    n_ref++;
    refno++;
  }

  // Read optimal origin offsets from fn_doc
  if (fn_doc!="") {
    DF.read(fn_doc);
    SF.go_beginning();
    while (!SF.eof()) {
      fn_tmp=SF.NextImg();
      if (DF.search_comment(fn_tmp)) {
	imgs_oldxoff.push_back(DF(3));
	imgs_oldyoff.push_back(DF(4));
      } else {
	REPORT_ERROR(1,(string)"Prog_MLalign2D_prm: Cannot find "+fn_tmp+" in docfile "+fn_doc);
      }
    }
    DF.go_beginning();
  }

  // For limited orientational search: fill imgs_oldphi & imgs_oldtheta 
  // (either read from fn_doc or initialize to -999.)
  if (limit_rot) {
    imgs_oldphi.clear();
    imgs_oldtheta.clear();
    SF.go_beginning();
    while (!SF.eof()) {
      fn_tmp=SF.NextImg();
      if (fn_doc!="") {
	if (DF.search_comment(fn_tmp)) {
	  imgs_oldphi.push_back(DF(0));
	  imgs_oldtheta.push_back(DF(1));
	} else {
	  REPORT_ERROR(1,(string)"Prog_MLalign2D_prm: Cannot find "+fn_tmp+" in docfile "+fn_doc);
	}
      } else {
	imgs_oldphi.push_back(-999.);
	imgs_oldtheta.push_back(-999.);
      }
    }
  }
  DF.clear();

  // read in model fractions if given on command line
  if (fn_frac!="") {
    DF.read(fn_frac);
    DF.go_first_data_line();
    for (refno=0; refno<n_ref; refno++) {
      DL=DF.get_current_line();
      alpha_k[refno]=DL[0];
      if (do_mirror) {
	if (DL[1]>1.||DL[1]<0.) 
	  REPORT_ERROR(1,"Prog_MLalign2D_prm: Mirror fraction (2nd column) should be [0,1]!");
	mirror_fraction[refno]=DL[1];
      }
      sumfrac+=alpha_k[refno];
      DF.next_data_line();
    }
    if (ABS(sumfrac-1.)>1e-3) 
      if (verb>0) cerr << " ->WARNING: Sum of all expected model fractions ("<<sumfrac<<") is not one!"<<endl;
    for (refno=0; refno<n_ref; refno++) { alpha_k[refno]/=sumfrac; }
  }

}

void Prog_MLalign2D_prm::write_offsets(FileName fn, vector<double> &data) {

  ofstream fh;
  int itot;

  fh.open((fn).c_str(),ios::out);
  if (!fh) {
    fh.clear();
    // Create the directory if it does not exist yet, and try again
    string dirname;
    int last_slash=((string)fn).rfind("/");
    dirname=((string)fn).erase(last_slash);
    if (!exists(dirname)) system(((string)"mkdir -p "+dirname).c_str());
    fh.open((fn).c_str(),ios::out);
    if (!fh) 
      REPORT_ERROR(3008,(string)"Prog_MLalign2D_prm: Cannot write file: "+fn);
  }
  
  itot=data.size();
  fh<<itot<<"\n";
  for (int i=0; i<itot; i+=2) {
    fh<<data[i]<<" "<<data[i+1]<<"\n";
  }
  fh.close();
  data.clear();

}

bool Prog_MLalign2D_prm::read_offsets(FileName fn, vector<double> &data) {

  ifstream fh;
  int ii,itot,nr_off,itoth,nr_offh;
  double remain;
  vector<double> data1;

  if (!exists(fn)) return false;
  else {
    fh.open((fn).c_str(),ios::in);
    if (!fh) return false;
    else {
      fh>>itot;
      if (do_mirror) nr_off=n_ref*4;
      else nr_off=n_ref*2;
      if (itot!=nr_off) {
	fh.close();
	return false;
      } else {
        data.clear();
        data.resize(itot);
        for (int i=0; i<itot; i+=2) {
          fh>>data[i];
          fh>>data[i+1];
        }
	fh.close(); 
	return true;
      }
    }
  }

}

// Calculate probability density function of all in-plane transformations phi
void Prog_MLalign2D_prm::calculate_pdf_phi() {

  double r2,pdfpix,sum;
  P_phi.resize(dim,dim);
  P_phi.set_Xmipp_origin();
  Mr2.resize(dim,dim);
  Mr2.set_Xmipp_origin();

  FOR_ALL_ELEMENTS_IN_MATRIX2D(P_phi) {
    r2=(double)(j*j + i*i);
    if (sigma_offset>0.) {
      pdfpix=exp(-r2/(2*sigma_offset*sigma_offset));
      pdfpix/=2*PI*sigma_offset*sigma_offset*nr_psi*nr_nomirror_flips;
    } else {
      if (j==0 && i==0) pdfpix=1.;
      else pdfpix=0.;
    }
    MAT_ELEM(P_phi,i,j)=pdfpix;
    MAT_ELEM(Mr2,i,j)=(float)r2;
  }

}

// Rotate reference for all models and rotations and fill Fref vectors =============
void Prog_MLalign2D_prm::rotate_reference(vector<ImageXmipp> &Iref, 
						     bool fill_real_space, 
						     bool fill_fourier_space, 
						     vector <vector< matrix2D<double> > > &Mref,
						     vector <vector< matrix2D<complex<double> > > > &Fref) {

  double AA,stdAA,psi,dum,avg,mean_ref,stddev_ref,dummy;
  matrix2D<double> Maux;
  matrix2D<complex<double> > Faux;
  vector<matrix2D<complex <double> > > dumF;
  vector<matrix2D<double> > dumM;
  matrix2D<int> mask,omask;
  matrix2D<double> cmask;

  Maux.init_zeros(dim,dim);
  Maux.set_Xmipp_origin();
  Fref.clear();
  Mref.clear();
  A2.clear();

  // prepare masks
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  BinaryCircularMask(mask,hdim,INNER_MASK);
  omask.resize(dim,dim);
  omask.set_Xmipp_origin();
  BinaryCircularMask(omask,hdim,OUTSIDE_MASK);

  FOR_ALL_MODELS() {
    Mref.push_back(dumM);
    Fref.push_back(dumF);
    compute_stats_within_binary_mask(omask,Iref[refno](),dum,dum,avg,dum);
    FOR_ALL_ROTATIONS() {
      // Add arbitrary number (small_angle) to avoid 0-degree rotation (lacking interpolation)
      psi=(double)(ipsi*psi_max/nr_psi)+SMALLANGLE;
      Maux=Iref[refno]().rotate_Bspline(3,psi,WRAP);
      apply_binary_mask(mask,Maux,Maux,avg);
      // Normalize the magnitude of the rotated references to 1st rot of that ref
      // This is necessary because interpolation due to rotation can lead to lower overall Fref 
      // This would result in lower probabilities for those rotations
      AA=Maux.sum2();      
      if (ipsi==0) {
	stdAA=AA;
	A2.push_back(AA);
      }
      // Subtract mean_ref from image prior to FFT for maxCC
      if (maxCC_rather_than_ML) {
	Maux.compute_stats(mean_ref,stddev_ref,dummy,dummy);
	Maux-=mean_ref;
	if (ipsi==0) A2[refno]=stddev_ref;
      }
      if (fill_real_space) {
	Mref[refno].push_back(Maux);
	if (AA>0) Mref[refno][ipsi]*=sqrt(stdAA/AA);
      }
      if (fill_fourier_space) {
	FourierTransformHalf(Maux,Faux);
	Faux*=dim*dim;
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux) {
	  dMij(Faux,i,j)=conj(dMij(Faux,i,j));
	}
	Fref[refno].push_back(Faux);
	if (AA>0) Fref[refno][ipsi]*=sqrt(stdAA/AA);
      }
    }
    // If we dont use save_mem1 Iref[refno] is useless from here on
    if (!save_mem1) Iref[refno]().resize(0,0);
  }

}

// Collect all rotations and sum to update Iref() for all models ==========
void Prog_MLalign2D_prm::reverse_rotate_reference(
	     vector <vector< matrix2D<complex<double> > > > &Fref, 
	     vector <vector< matrix2D<double> > > &Mref, bool real_space, 
	     vector<matrix2D<double > > &Mnew) {

  double psi,dum,avg,ang;
  matrix2D<double> Maux,Maux2;
  matrix2D<complex<double> > Faux;
  matrix2D<int> mask, omask;
  Maux.resize(dim,dim);
  Maux2.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Maux2.set_Xmipp_origin();

  Mnew.clear();
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  BinaryCircularMask(mask,hdim,INNER_MASK);
  omask.resize(dim,dim);
  omask.set_Xmipp_origin();
  BinaryCircularMask(omask,hdim,OUTSIDE_MASK);

  FOR_ALL_MODELS() {
    Maux.init_zeros();
    Mnew.push_back(Maux);
    FOR_ALL_ROTATIONS() {
      // Add arbitrary number to avoid 0-degree rotation without interpolation effects
      psi=(double)(ipsi*psi_max/nr_psi)+SMALLANGLE;
      if (real_space) {
	Maux=Mref[refno][ipsi];
      } else {
	InverseFourierTransformHalf(Fref[refno][ipsi],Maux,dim);
	Maux/=dim*dim;
	CenterFFT(Maux,true);
      }
      compute_stats_within_binary_mask(omask,Maux,dum,dum,avg,dum);
      Maux2=Maux.rotate_Bspline(3,-psi,WRAP);
      apply_binary_mask(mask,Maux2,Maux2,avg);
      Mnew[refno]+=Maux2;
    }  

    // perform correction of origin-pixel artifacts
    if (do_esthetics)
      MAT_ELEM(Mnew[refno],0,0)=(MAT_ELEM(Mnew[refno],1,0)+MAT_ELEM(Mnew[refno],0,1)+
				 MAT_ELEM(Mnew[refno],-1,0)+MAT_ELEM(Mnew[refno],0,-1))/4;

  }

}

void Prog_MLalign2D_prm::preselect_directions(float &phi, float &theta,
							 vector<double> &pdf_directions) {

  float phi_ref,theta_ref,angle,angle2;
  matrix1D<double> u,v;

  pdf_directions.clear();
  pdf_directions.resize(n_ref);
  FOR_ALL_MODELS() {
    if (!limit_rot || (phi==-999. && theta==-999.) ) pdf_directions[refno]=1.;
    else {
      phi_ref=Iref[refno].Phi();
      theta_ref=Iref[refno].Theta();
      Euler_direction(phi,theta,0.,u);
      Euler_direction(phi_ref,theta_ref,0.,v);
      u.normalize();
      v.normalize();
      angle=RAD2DEG(acos(dot_product(u,v)));
      angle=fabs(realWRAP(angle,-180,180));
      // also check mirror
      angle2=180.+angle;
      angle2=fabs(realWRAP(angle2,-180,180));
      angle=MIN(angle,angle2);
      if (fabs(angle)>search_rot) pdf_directions[refno]=0.;
      else pdf_directions[refno]=1.;
    }
  }

}

// Pre-selection of significant refno and ipsi, based on current optimal translation =======
void Prog_MLalign2D_prm::preselect_significant_model_phi(
			 matrix2D<double> &Mimg, vector<double > &offsets, 
			 vector <vector< matrix2D<double > > > &Mref, 
			 matrix2D<int> &Msignificant,
			 vector<double> &pdf_directions) {


  matrix2D<double> Maux,Maux2,Mrefl,Mdsig(n_ref,nr_psi*nr_flip);
  double ropt,sigma_noise2,aux,fracpdf;
  double Xi2,A2_plus_Xi2,CC;
  double mindiff=99.e99;
  double maxweight=-99.e99;
  int irot,irefmir;
  vector<double> maxw_ref(2*n_ref);
  matrix1D<double> trans(2);
  vector<matrix2D<double> > Mrot;

  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Maux2.resize(dim,dim);
  Maux2.set_Xmipp_origin();
  sigma_noise2=sigma_noise*sigma_noise;
  Xi2=Mimg.sum2();
  Msignificant.init_zeros();

  // Flip images and calculate correlations and maximum correlation
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      A2_plus_Xi2=0.5*(A2[refno]+Xi2);
      if (save_mem1) {
	Mrot.clear();
	FOR_ALL_ROTATIONS() {
	  double psi=(double)(ipsi*psi_max/nr_psi)+SMALLANGLE; 
	  Maux=Iref[refno]().rotate_Bspline(3,psi,WRAP);
	  Mrot.push_back(Maux);
	}
      }
      FOR_ALL_FLIPS() {
	irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	// Do not trust optimal offsets if they are larger than 3*sigma_offset: 
	ropt=sqrt(offsets[2*irefmir]*offsets[2*irefmir]+offsets[2*irefmir+1]*offsets[2*irefmir+1]);
	if (ropt>3*sigma_offset)  {
	  FOR_ALL_ROTATIONS() {
	    irot=iflip*nr_psi+ipsi;
	    dMij(Msignificant,refno,irot)=1;
	  }
	} else {
	  trans(0)=offsets[2*irefmir];
	  trans(1)=offsets[2*irefmir+1];
	  Maux=Mimg.translate(trans,true);
	  apply_geom(Maux2,F[iflip],Maux,IS_INV,WRAP);
	  FOR_ALL_ROTATIONS() {
	    irot=iflip*nr_psi+ipsi;
	    dMij(Msignificant,refno,irot)=0;
	    CC=A2_plus_Xi2;
	    if (save_mem1) Mrefl=Mrot[ipsi];
	    else Mrefl=Mref[refno][ipsi];
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Maux2) {
	      CC-=dMij(Maux2,i,j)*dMij(Mrefl,i,j);
	    }
	    dMij(Mdsig,refno,irot)=CC;
	    if (CC<mindiff) mindiff=CC;
	  }
	}
      }
    }
  }

  // Now that we have mindiff calculate the weighting matrices and maxweight
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      FOR_ALL_FLIPS() {
	irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  if (!dMij(Msignificant,refno,irot)) {
	    if (iflip<nr_nomirror_flips) fracpdf=alpha_k[refno]*(1.-mirror_fraction[refno]);
	    else fracpdf=alpha_k[refno]*mirror_fraction[refno];
	    aux=(dMij(Mdsig,refno,irot)-mindiff)/sigma_noise2;
	    // next line because of numerical precision of exp-function
	    if (aux>1000.) aux=0.;
	    else aux=exp(-aux)*fracpdf*MAT_ELEM(P_phi,(int)offsets[2*irefmir+1],(int)offsets[2*irefmir]);
	    dMij(Mdsig,refno,irot)=aux;
	    if (aux>maxw_ref[irefmir]) maxw_ref[irefmir]=aux;
	  }
	}
      }
    }
  }

  // Now that we have maxweight calculate which weighting matrices are significant
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      FOR_ALL_FLIPS() {
	irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  if (!dMij(Msignificant,refno,irot)) {
	    if (dMij(Mdsig,refno,irot)>=C_fast*maxw_ref[irefmir]) dMij(Msignificant,refno,irot)=1;
	    else dMij(Msignificant,refno,irot)=0;
	  }
	}
      }
    }
  }

}

// Calculate translated matrices for all limited translations
// for each of the flipped variants
void Prog_MLalign2D_prm::calculate_realspace_offsets(
					      matrix2D<double> &Mimg, vector<double > &offsets,
					      vector<double > &pdf_directions,
					      vector<vector<matrix2D<double> > > &Mimg_trans,
					      matrix2D<int> &Moffsets, matrix2D<int> &Moffsets_mirror) {

  int irefmir,ix,iy,opt_ix,opt_iy,iflip,irot,opt_iflip,nr_mir,iflip_start,iflip_stop,count;
  double ropt2,maxCC,CC,dxx,dyy;
  vector<matrix2D<double> > Mflip, dum;
  matrix1D<double> trans(2);
  matrix2D<double> Maux2,Maux;
  vector<matrix2D<double> > Finv;

  if (do_mirror) nr_mir=2;
  else nr_mir=1;

  Moffsets.resize(dim,dim);
  Moffsets.set_Xmipp_origin();
  Moffsets.init_constant(-1);
  Moffsets_mirror.resize(dim,dim);
  Moffsets_mirror.set_Xmipp_origin();
  Moffsets_mirror.init_constant(-1);
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Mimg_trans.clear();

  // Pre-calculate flipped image matrices
  Mflip.clear();
  FOR_ALL_FLIPS() {
    Maux.set_Xmipp_origin();
    apply_geom(Maux,F[iflip],Mimg,IS_INV,WRAP);
    Mflip.push_back(Maux);
  }

  // Calculate inverted flipping matrices
  // This is to have the offsets consistent with those in fast_mode
  FOR_ALL_FLIPS() {
    Maux=F[iflip].inv();
    Finv.push_back(Maux);
  }

  // If offsets > max_shift: reset to zero...
  count=0;
  FOR_ALL_MODELS() {
    for (int imir=0; imir<nr_mir; imir++) {
      irefmir=imir*n_ref+refno;
      iflip_start=imir*nr_nomirror_flips;
      iflip_stop=imir*nr_nomirror_flips+nr_nomirror_flips;
      ropt2=offsets[2*irefmir]*offsets[2*irefmir]+offsets[2*irefmir+1]*offsets[2*irefmir+1];
      if (ropt2>max_shift*max_shift) {
	offsets[2*irefmir]=0.;
	offsets[2*irefmir+1]=0.;
      }
      if (!limit_rot || pdf_directions[refno]>0.) {
	FOR_ALL_LIMITED_TRANSLATIONS() {
	  ix=ROUND(offsets[2*irefmir]+Vtrans[itrans](0));
	  iy=ROUND(offsets[2*irefmir+1]+Vtrans[itrans](1));
	  dxx=(double)intWRAP(ix,Moffsets.startingX(),Moffsets.finishingX());
	  dyy=(double)intWRAP(iy,Moffsets.startingY(),Moffsets.finishingY());
	  // For non-mirrors
	  if (imir==0 && MAT_ELEM(Moffsets,ROUND(dyy),ROUND(dxx))<0) {
	    Mimg_trans.push_back(dum);
	    for (int iflip=iflip_start; iflip<iflip_stop; iflip++) {
	      trans(0)=dxx*DIRECT_MAT_ELEM(Finv[iflip],0,0)+dyy*DIRECT_MAT_ELEM(Finv[iflip],0,1);
	      trans(1)=dxx*DIRECT_MAT_ELEM(Finv[iflip],1,0)+dyy*DIRECT_MAT_ELEM(Finv[iflip],1,1);
	      Maux=Mflip[iflip].translate(trans,WRAP);
	      Mimg_trans[count].push_back(Maux);
	    }
	    MAT_ELEM(Moffsets,ROUND(dyy),ROUND(dxx))=count;
	    count++;
	  }
	  // For mirrors use a separate offset-matrix
	  else if (imir==1 && MAT_ELEM(Moffsets_mirror,ROUND(dyy),ROUND(dxx))<0) {
	    Mimg_trans.push_back(dum);
	    for (int iflip=iflip_start; iflip<iflip_stop; iflip++) {
	      trans(0)=dxx*DIRECT_MAT_ELEM(Finv[iflip],0,0)+dyy*DIRECT_MAT_ELEM(Finv[iflip],0,1);
	      trans(1)=dxx*DIRECT_MAT_ELEM(Finv[iflip],1,0)+dyy*DIRECT_MAT_ELEM(Finv[iflip],1,1);
	      Maux=Mflip[iflip].translate(trans,WRAP);
	      Mimg_trans[count].push_back(Maux);
	    }
	    MAT_ELEM(Moffsets_mirror,ROUND(dyy),ROUND(dxx))=count;
	    count++;
	  }
	}
      }
    }
  }

}

void Prog_MLalign2D_prm::ML_integrate_locally(
          matrix2D<double> &Mimg, vector <vector< matrix2D<double> > > &Mref, 
	  vector <vector< matrix2D<double> > > &Mwsum_imgs, 
	  double &wsum_sigma_noise, double &wsum_sigma_offset, 
	  vector<double> &sumw, vector<double> &sumw_mirror, 
	  double &LL, double &fracweight, int &opt_refno, double &opt_psi, 
	  matrix1D<double> &opt_offsets, vector<double> &opt_offsets_ref,
	  vector<double> &pdf_directions) {

  matrix3D<double> Mweight;
  matrix2D<double> Maux,Mdzero,Mtrans;
  vector<matrix2D<double> > Mflip;
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Mtrans.resize(Maux);
  matrix2D<int> Moffsets,Moffsets_mirror;
  vector<vector<matrix2D<double> > > Mimg_trans;
  vector<double> refw(n_ref), refw_mirror(n_ref);
  double sigma_noise2,Xi2,aux,fracpdf,A2_plus_Xi2,weight,maxw,pdf,diff,mindiff2=99.e99;
  double sum,wsum_corr=0., sum_refw=0., wsum_A2=0., maxweight=-99.e99;
  int point_trans,nr_mir,irefmir,iflip_start,iflip_stop,ixx,iyy,irot;
  int opt_ipsi=0,opt_iflip=0,opt_irefmir=0,opt_itrans=0;
  int ii,imax=n_ref*nr_flip/nr_nomirror_flips;
  vector<double> Pmax_refmir(imax);
  for (int i=0; i<imax; i++) Pmax_refmir[i]=-99.e99;
  sigma_noise2=sigma_noise*sigma_noise;

  // Calculate all flipped and translated versions of Mimg
  calculate_realspace_offsets(Mimg,opt_offsets_ref,pdf_directions, 
			      Mimg_trans,Moffsets,Moffsets_mirror);

  // Calculate all squared differences & mindiff2 (for optimal offsets only)
  Mweight.init_zeros(nr_trans,n_ref,nr_flip*nr_psi);
  Xi2=Mimg.sum2();
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      A2_plus_Xi2=0.5*(A2[refno]+Xi2);
      FOR_ALL_FLIPS() {
	irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	ixx=ROUND(opt_offsets_ref[2*irefmir]);
	iyy=ROUND(opt_offsets_ref[2*irefmir+1]);
	if (iflip<nr_nomirror_flips) point_trans=MAT_ELEM(Moffsets,iyy,ixx);
	else point_trans=MAT_ELEM(Moffsets_mirror,iyy,ixx);
	Mtrans=Mimg_trans[point_trans][iflip%nr_nomirror_flips];
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  Maux=Mref[refno][ipsi];
	  diff=A2_plus_Xi2;
	  for (int ii=0; ii<dim2; ii++) {
	    diff-=(Maux).__m[ii]*(Mtrans).__m[ii];
	  }
	  dVkij(Mweight,zero_trans,refno,irot)=diff;
	  if (diff<mindiff2) mindiff2=diff;
	}
      }
    }
  }
  
  // Now that we have mindiff2 calculate the weighting matrices and maxweight
  FOR_ALL_MODELS() {
    refw[refno]=0.;
    refw_mirror[refno]=0.;
    if (!limit_rot || pdf_directions[refno]>0.) {
      FOR_ALL_FLIPS() {
	irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	ixx=ROUND(opt_offsets_ref[2*irefmir]);
	iyy=ROUND(opt_offsets_ref[2*irefmir+1]);
	if (iflip<nr_nomirror_flips) fracpdf=(1.-mirror_fraction[refno]);
	else fracpdf=mirror_fraction[refno];
	fracpdf*=alpha_k[refno]*pdf_directions[refno]*MAT_ELEM(P_phi,iyy,ixx);
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  diff=dVkij(Mweight,zero_trans,refno,irot);
	  aux=(diff-mindiff2)/sigma_noise2;
	  // next line because of numerical precision of exp-function
	  if (aux>1000.) weight=0.;
	  else weight=exp(-aux)*fracpdf;
	  wsum_corr+=weight*diff;
	  dVkij(Mweight,zero_trans,refno,irot)=weight;
	  if (weight>Pmax_refmir[irefmir]) Pmax_refmir[irefmir]=weight;
	  if (weight>maxweight) {
	    maxweight=weight;
	    opt_refno=refno;
	    opt_ipsi=ipsi;
	    opt_iflip=iflip;
	    opt_itrans=zero_trans;
	    opt_irefmir=irefmir;
	  }
	  // Accumulate sum weights
	  if (iflip<nr_nomirror_flips) refw[refno]+=weight;
	  else refw_mirror[refno]+=weight;  
	  sum_refw+=weight;
	}
      }
    }
  }
  
  // Now for all irefmir, check significant rotations
  // and calculate their limited_translations probabilities
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      A2_plus_Xi2=0.5*(A2[refno]+Xi2);
      FOR_ALL_FLIPS() {
	irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	if (iflip<nr_nomirror_flips) fracpdf=(1.-mirror_fraction[refno]);
	else fracpdf=mirror_fraction[refno];
	fracpdf*=alpha_k[refno]*pdf_directions[refno];
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  if (dVkij(Mweight,zero_trans,refno,irot)>C_fast*Pmax_refmir[irefmir]) {
	    Maux=Mref[refno][ipsi];
	    // expand for all limited translations
	    FOR_ALL_LIMITED_TRANSLATIONS() {
	      if (itrans!=zero_trans) {
		ixx=ROUND(opt_offsets_ref[2*irefmir]+Vtrans[itrans](0));
		iyy=ROUND(opt_offsets_ref[2*irefmir+1]+Vtrans[itrans](1));
		if (iflip<nr_nomirror_flips) point_trans=MAT_ELEM(Moffsets,iyy,ixx);
		else point_trans=MAT_ELEM(Moffsets_mirror,iyy,ixx);
		Mtrans=Mimg_trans[point_trans][iflip%nr_nomirror_flips];
		diff=A2_plus_Xi2;
		for (int ii=0; ii<dim2; ii++) {
		  diff-=(Maux).__m[ii]*(Mtrans).__m[ii];
		}
		aux=(diff-mindiff2)/sigma_noise2;
		// next line because of numerical precision of exp-function
		if (aux>1000.) weight=0.;
		else weight=exp(-aux)*fracpdf*MAT_ELEM(P_phi,iyy,ixx);
		wsum_corr+=weight*diff;
		dVkij(Mweight,itrans,refno,irot)=weight;
		if (weight>maxweight) {
		  maxweight=weight;
		  opt_refno=refno;
		  opt_ipsi=ipsi;
		  opt_iflip=iflip;
		  opt_itrans=itrans;
		  opt_irefmir=irefmir;
		}
		// Accumulate sum weights
		if (iflip<nr_nomirror_flips) refw[refno]+=weight;
		else refw_mirror[refno]+=weight; 
		sum_refw+=weight;
	      }
	    }
	  }
	}
      }
    }
  }

  // Normalize all weighted sums by sum_refw such that sum over all weights is one!
  // And accumulate weighted sums
  wsum_sigma_noise+=(2*wsum_corr/sum_refw);
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      sumw[refno]+=(refw[refno]+refw_mirror[refno])/sum_refw;
      sumw_mirror[refno]+=refw_mirror[refno]/sum_refw;
      FOR_ALL_FLIPS() {
	irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  FOR_ALL_LIMITED_TRANSLATIONS() {
	    weight=dVkij(Mweight,itrans,refno,irot);
	    if (weight>SIGNIFICANT_WEIGHT_LOW*maxweight) {
	      weight/=sum_refw;
	      ixx=ROUND(opt_offsets_ref[2*irefmir]+Vtrans[itrans](0));
	      iyy=ROUND(opt_offsets_ref[2*irefmir+1]+Vtrans[itrans](1));
	      if (iflip<nr_nomirror_flips) point_trans=MAT_ELEM(Moffsets,iyy,ixx);
	      else point_trans=MAT_ELEM(Moffsets_mirror,iyy,ixx);
	      // weighted sums
	      wsum_sigma_offset+=weight*(double)(ixx*ixx+iyy*iyy);
	      Mwsum_imgs[refno][ipsi]+=weight*Mimg_trans[point_trans][iflip%nr_nomirror_flips];
	    }
	  }
	}
      }
    }
  }

  // Update the optimal origin in-plane transformations
  opt_offsets(0)=opt_offsets_ref[2*opt_irefmir]+Vtrans[opt_itrans](0);
  opt_offsets(1)=opt_offsets_ref[2*opt_irefmir+1]+Vtrans[opt_itrans](1);
  opt_psi=-psi_step*(opt_iflip*nr_psi+opt_ipsi)-SMALLANGLE;
  fracweight=maxweight/sum_refw; 

  if (do_mirror) nr_mir=2;
  else nr_mir=1;
  for (int i=0; i<imax; i++) Pmax_refmir[i]=-99.e99;
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      for (int imir=0; imir<nr_mir; imir++) {
	irefmir=imir*n_ref+refno;
	iflip_start=imir*nr_nomirror_flips;
	iflip_stop=imir*nr_nomirror_flips+nr_nomirror_flips;
	opt_itrans=zero_trans;
	for (int iflip=iflip_start; iflip<iflip_stop; iflip++) {
	  FOR_ALL_ROTATIONS() {
	    irot=iflip*nr_psi+ipsi;
	    FOR_ALL_LIMITED_TRANSLATIONS() {
	      weight=dVkij(Mweight,itrans,refno,irot);
	      if (weight>Pmax_refmir[irefmir]) {
		Pmax_refmir[irefmir]=weight;
		opt_itrans=itrans;
	      }
	    }
	  }
	}
	opt_offsets_ref[2*irefmir]+=Vtrans[opt_itrans](0);
	opt_offsets_ref[2*irefmir+1]+=Vtrans[opt_itrans](1);
      }
    }
  }

  // Compute Log Likelihood
  // 1st term: log(refw_i)
  // 2nd term: for subtracting mindiff2
  // 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
  LL+= log(sum_refw) - mindiff2/sigma_noise2 - dim*dim*log(2.50663*sigma_noise);

}


// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void Prog_MLalign2D_prm::ML_integrate_complete(
          matrix2D<double> &Mimg, vector <vector< matrix2D<complex<double> > > > &Fref, 
          matrix2D<int> &Msignificant,
	  vector <vector< matrix2D<complex<double> > > > &Fwsum_imgs, 
	  double &wsum_sigma_noise, double &wsum_sigma_offset, 
	  vector<double> &sumw, vector<double> &sumw_mirror, 
	  double &LL, double &fracweight, int &opt_refno, double &opt_psi, 
	  matrix1D<double> &opt_offsets, vector<double> &opt_offsets_ref,
	  vector<double> &pdf_directions) {

  matrix2D<double> Maux,Mdzero;
  matrix2D<complex<double> > Fimg,Faux,Faux2;
  vector<matrix2D<complex<double> > > Fimg_flip;
  vector<matrix2D<double> > dumM;
  vector<bool> dumb;
  vector <vector< matrix2D<double> > > Mweight;
  vector<double> refw(n_ref), refw_mirror(n_ref);
  double sigma_noise2,XiA,Xi2,aux,fracpdf,A2_plus_Xi2,maxw,mind,mindiff2=99.e99;
  double sum,wsum_corr=0., sum_refw=0., wsum_A2=0., maxweight=-99.e99;
  int irot,irefmir,sigdim,xmax,ymax;
  int ioptx=0,iopty=0,ioptpsi=0,ioptflip=0,imax=0;
  if (fast_mode) imax=n_ref*nr_flip/nr_nomirror_flips;
  vector<int> ioptx_ref(imax),iopty_ref(imax),ioptflip_ref(imax);
  vector<double> maxw_ref(imax);


  /* Not to store all 360-degrees rotations of the references (and pdf, Fwsum_imgs etc.) in memory,
     the experimental image is rotated over 0, 90, 180 & 270 degrees (called FLIPS), and only
     rotations over 90 degrees of the references are stored.
     This save a factor of 4 in memory requirements, and these FLIPS do not require interpolation
     of the experiemental image and thereby deterioration of the process.
     If do_mirror there are 8 flips, now also including all mirrored versions.
     Rotation is done in real space, and then the Fourier Transform is calculated four/eight times.
     In principle, rotation could be done in Fourier space, but then there is a problem with
     even-sized images, where the origin is not exactly in the center, and 1 pixel wrapping is required.
     Anyway, the total number of (I)FFT's is determined in much greater extent by n_ref and n_rot!
  */
  

  // Only translations smaller than save_mem2 (default=6) sigma_offset are considered!
  // This saves a lot of memory and CPU! (typically a factor 2, depending on sigma_offset vs. dim)
  if (save_mem2) sigdim=2*CEIL(sigma_offset*3);
  else sigdim=2*CEIL(sigma_offset*6);
  sigdim++; // (to get uneven number)
  sigdim=MIN(dim,sigdim);

  if (fast_mode) {
    for (int i=0; i<imax; i++) {
      maxw_ref[i]=-99.e99;
      ioptx_ref[i]=0;
      iopty_ref[i]=0;
      ioptflip_ref[i]=0;
    }
  }
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Fimg_flip.clear();
  Mweight.clear();
  sigma_noise2=sigma_noise*sigma_noise;  
  Xi2=Mimg.sum2();

  // Flip images and calculate correlation matrices and maximum correlation
  FOR_ALL_FLIPS() {
    Maux.set_Xmipp_origin();
    apply_geom(Maux,F[iflip],Mimg,IS_INV,WRAP);
    FourierTransformHalf(Maux,Fimg);
    Fimg*=dim*dim;
    Fimg_flip.push_back(Fimg);  
    FOR_ALL_MODELS() {
      Mweight.push_back(dumM); 
      if (!limit_rot || pdf_directions[refno]>0.) {
	A2_plus_Xi2=0.5*(A2[refno]+Xi2);
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  if (dMij(Msignificant,refno,irot)) {
	    mul_elements(Fimg,Fref[refno][ipsi],Faux);
	    InverseFourierTransformHalf(Faux,Maux,dim);
	    Maux/=dim*dim;
	    CenterFFT(Maux,true);
	    Mweight[refno].push_back(Mdzero);
	    Mweight[refno][irot].resize(sigdim,sigdim);
	    Mweight[refno][irot].set_Xmipp_origin();
	    FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot]) {
	      MAT_ELEM(Mweight[refno][irot],i,j)=A2_plus_Xi2-MAT_ELEM(Maux,i,j);
	    }
	    mind=Mweight[refno][irot].compute_min();
	    if (mind<mindiff2) mindiff2=mind;
	  } else {
	    Mweight[refno].push_back(Mdzero);
	  }
	}
      }
    }
  }

  // Now that we have mindiff2 calculate the weighting matrices and maxweight
  FOR_ALL_MODELS() {
    refw[refno]=0.;
    refw_mirror[refno]=0.;
    if (!limit_rot || pdf_directions[refno]>0.) {
      FOR_ALL_ROTATIONS() {
	FOR_ALL_FLIPS() {
	  irot=iflip*nr_psi+ipsi;
	  irefmir=FLOOR(iflip/nr_nomirror_flips)*n_ref+refno;
	  if (dMij(Msignificant,refno,irot)) {
	    if (iflip<nr_nomirror_flips) fracpdf=alpha_k[refno]*(1.-mirror_fraction[refno]);
	    else fracpdf=alpha_k[refno]*mirror_fraction[refno];
	    sum=0.;
	    FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot]) {
	      aux=(MAT_ELEM(Mweight[refno][irot],i,j)-mindiff2)/sigma_noise2;
	      // next line because of numerical precision of exp-function
	      if (aux>1000.) aux=0.;
	      else aux=exp(-aux)*fracpdf*MAT_ELEM(P_phi,i,j);
	      wsum_corr+=aux*MAT_ELEM(Mweight[refno][irot],i,j);
	      MAT_ELEM(Mweight[refno][irot],i,j)=aux;
	      sum+=aux;
	    }
	    if (iflip<nr_nomirror_flips) refw[refno]+=sum;
	    else refw_mirror[refno]+=sum;
	    Mweight[refno][irot].max_index(ymax,xmax);
	    maxw=MAT_ELEM(Mweight[refno][irot],ymax,xmax);
	    if (maxw>maxweight) {
	      maxweight=maxw;
	      iopty=ymax;
	      ioptx=xmax;
	      ioptpsi=ipsi;
	      ioptflip=iflip;
	      opt_refno=refno;
	    }
	    if (fast_mode && maxw>maxw_ref[irefmir]) {
	      maxw_ref[irefmir]=maxw;
	      iopty_ref[irefmir]=ymax;
	      ioptx_ref[irefmir]=xmax;
	      ioptflip_ref[irefmir]=iflip;
	    } 
	  }
	}
      }
    }
    sum_refw+=refw[refno]+refw_mirror[refno];
  }

  // Normalize all weighted sums by sum_refw such that sum over all weights is one!
  // And accumulate the FT of the weighted, shifted images.
  wsum_sigma_noise+=(2*wsum_corr/sum_refw);
  FOR_ALL_MODELS() {
    if (!limit_rot || pdf_directions[refno]>0.) {
      sumw[refno]+=(refw[refno]+refw_mirror[refno])/sum_refw;
      sumw_mirror[refno]+=refw_mirror[refno]/sum_refw;
      FOR_ALL_ROTATIONS() {
	FOR_ALL_FLIPS() {
	  irot=iflip*nr_psi+ipsi;
	  if (dMij(Msignificant,refno,irot)) {
	    if (Mweight[refno][irot].compute_max()>SIGNIFICANT_WEIGHT_LOW*maxweight) {
	      Mweight[refno][irot]/=sum_refw;
	      // Use Maux, because Mweight is smaller than dim x dim!
	      Maux.init_zeros();
	      FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot]) {
		MAT_ELEM(Maux,i,j)=MAT_ELEM(Mweight[refno][irot],i,j);
		wsum_sigma_offset+=MAT_ELEM(Maux,i,j)*MAT_ELEM(Mr2,i,j);
	      }
	      FourierTransformHalf(Maux,Faux);
	      Faux*=dim*dim;
	      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux) {
		dMij(Fwsum_imgs[refno][ipsi],i,j)+=conj(dMij(Faux,i,j))*dMij(Fimg_flip[iflip],i,j);
	      }
	    }
	  }
	}
      }
    }
  }

  // Calculate optimal transformation parameters
  if (fast_mode) {
    for (int i=0; i<imax;i++) {
      opt_offsets_ref[2*i]=-(double)ioptx_ref[i]*DIRECT_MAT_ELEM(F[ioptflip_ref[i]],0,0)-
	                     (double)iopty_ref[i]*DIRECT_MAT_ELEM(F[ioptflip_ref[i]],0,1);
      opt_offsets_ref[2*i+1]=-(double)ioptx_ref[i]*DIRECT_MAT_ELEM(F[ioptflip_ref[i]],1,0)-
	                     (double)iopty_ref[i]*DIRECT_MAT_ELEM(F[ioptflip_ref[i]],1,1);
    }
  }
  opt_offsets(0)=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],0,0)-
                  (double)iopty*DIRECT_MAT_ELEM(F[ioptflip],0,1);
  opt_offsets(1)=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],1,0)-
                  (double)iopty*DIRECT_MAT_ELEM(F[ioptflip],1,1);
  opt_psi=-psi_step*(ioptflip*nr_psi+ioptpsi)-SMALLANGLE;
  fracweight=maxweight/sum_refw; 

  // Compute Log Likelihood
  // 1st term: log(refw_i)
  // 2nd term: for subtracting mindiff2
  // 3rd term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
  LL+= log(sum_refw) - mindiff2/sigma_noise2 - dim*dim*log(2.50663*sigma_noise);

}

void Prog_MLalign2D_prm::maxCC_search_complete(matrix2D<double> &Mimg, 
                vector <vector< matrix2D<complex<double> > > > &Fref, 
                vector <vector< matrix2D<double> > > &Mref, 
		double &search_shift, vector <vector< matrix2D<double> > > &Msum_imgs, 
		vector<double> &sumw, vector<double> &sumw_mirror, 
		double &maxCC, int &opt_refno, double &opt_psi, matrix1D<double> &opt_offsets,
		vector<double> &pdf_directions) {

  matrix2D<double> Maux,Maux2;
  matrix2D<complex<double> > Fimg, Faux;
  double sigma_noise2,aux,avg,std,CC;
  int irot,sigdim,xmax=0,ymax=0;
  int ioptx=0,iopty=0,ioptpsi=0,ioptflip=0,imax=0;
  double stddev_img,mean_img,dummy;

  maxCC=-99.e99; 
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Maux2.resize(dim,dim);
  Maux2.set_Xmipp_origin();
  Faux.resize(hdim+1,dim);
  sigma_noise2=sigma_noise*sigma_noise;
  matrix2D<int> shiftmask;

  if (search_shift>0.) {
    shiftmask.resize(dim,dim);
    shiftmask.set_Xmipp_origin();
    BinaryCircularMask(shiftmask,search_shift,INNER_MASK);
  }

  Mimg.compute_stats(mean_img,stddev_img,dummy,dummy);
  Maux2=Mimg;
  Maux2-=mean_img;

  // Flip images and calculate correlation matrices and maximum correlation
  FOR_ALL_FLIPS() {
    apply_geom(Maux,F[iflip],Maux2,IS_INV,WRAP);
    if (search_shift>0.) {
      FourierTransformHalf(Maux,Fimg);
      Fimg*=dim*dim;
    }
    FOR_ALL_MODELS() {
      if (!limit_rot || pdf_directions[refno]>0.) {
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  if (search_shift>0.) {
	    mul_elements(Fimg,Fref[refno][ipsi],Faux);
	    InverseFourierTransformHalf(Faux,Maux,dim);
	    Maux/=dim*dim;
	    CenterFFT(Maux,true);
	    apply_binary_mask(shiftmask,Maux,Maux,0.);
	    Maux.max_index(ymax,xmax);
	    CC=MAT_ELEM(Maux,ymax,xmax);
	  } else {
	    CC=0.;
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Maux) {
	      CC+=dMij(Maux,i,j)*dMij(Mref[refno][ipsi],i,j);
	    }
	  }
	  CC/=A2[refno]*stddev_img; // For maxCC-mode, A2[refno] holds stddev_ref!
	  if (CC>maxCC) {
	    maxCC=CC;
	    iopty=ymax;
	    ioptx=xmax;
	    ioptpsi=ipsi;
	    ioptflip=iflip;
	    opt_refno=refno;
	  }
	}
      }
    }
  }
  maxCC/=dim*dim;

  // Calculate optimal transformation parameters
  opt_offsets(0)=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],0,0)
                 -(double)iopty*DIRECT_MAT_ELEM(F[ioptflip],0,1);
  opt_offsets(1)=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],1,0)
                 -(double)iopty*DIRECT_MAT_ELEM(F[ioptflip],1,1);
  if (save_mem3)
    opt_psi=-ioptflip*360./nr_flip-SMALLANGLE;
  else
    opt_psi=-psi_step*(ioptflip*nr_psi+ioptpsi)-SMALLANGLE;

  // Store sums of the aligned images
  Mimg.translate(opt_offsets,Maux,true);
  apply_geom(Maux2,F[ioptflip],Maux,IS_INV,WRAP);
  Msum_imgs[opt_refno][ioptpsi]+=Maux2;
  sumw[opt_refno]+=1.;
  if (ioptflip>3) sumw_mirror[opt_refno]+=1.;

}


void Prog_MLalign2D_prm::ML_sum_over_all_images(SelFile &SF, vector<ImageXmipp> &Iref, 
			  double &LL, double &sumcorr, DocFile &DFo, 
			  vector<matrix2D<double> > &wsum_Mref,
			  double &wsum_sigma_noise, double &wsum_sigma_offset, 
			  vector<double> &sumw, vector<double> &sumw_mirror) {


  ImageXmipp img;
  SelLine line;
  FileName fn_img,fn_trans;
  vector <vector< matrix2D<double> > > Mref,Msum_imgs;
  vector <vector< matrix2D<complex<double> > > > Fref,Fwsum_imgs;
  vector<matrix2D<complex <double> > > dum;
  vector<matrix2D<double> > dum2;
  vector<double> allref_offsets,pdf_directions(n_ref);
  matrix2D<complex <double> > Fdzero;
  matrix2D<double>  Mdzero;
  matrix2D<int> Msignificant;
  Msignificant.resize(n_ref,nr_psi*nr_flip);
  matrix1D<double> dataline(8),opt_offsets(2),trans(2);  

  float old_phi=-999., old_theta=-999.;
  double opt_psi,opt_flip,maxcorr;
  double opt_xoff, opt_yoff;
  int c,nn,imgno,opt_refno;
  bool fill_real_space,fill_fourier_space;

  // Generate (FT of) each rotated version of all references 
  if ( limit_trans || (maxCC_rather_than_ML && !(search_shift>0.)) ) {
    fill_real_space=true;
    fill_fourier_space=false;
  }
  else if (fast_mode) {
    fill_fourier_space=true;
    if (save_mem1) fill_real_space=false;
    else fill_real_space=true;
  } else {
    fill_real_space=false;
    fill_fourier_space=true;
  }
  rotate_reference(Iref,fill_real_space,fill_fourier_space,Mref,Fref);

  // Initialize
  nn=SF.ImgNo();
  if (verb>0) init_progress_bar(nn);
  c=MAX(1,nn/60);
  Fwsum_imgs.clear();
  Msum_imgs.clear();
  sumw.clear();
  sumw_mirror.clear();
  Fdzero.resize(hdim+1,dim);
  Mdzero.resize(dim,dim);
  Mdzero.set_Xmipp_origin();
  LL=0.;
  wsum_sigma_noise=0.;
  wsum_sigma_offset=0.;
  sumcorr=0.;
  trans.init_zeros();
  FOR_ALL_MODELS() {
    sumw.push_back(0.);
    sumw_mirror.push_back(0.);
    if (maxCC_rather_than_ML || limit_trans ) {
      Msum_imgs.push_back(dum2);
      FOR_ALL_ROTATIONS() {Msum_imgs[refno].push_back(Mdzero);}
    } else {
      Fwsum_imgs.push_back(dum);
      FOR_ALL_ROTATIONS() {Fwsum_imgs[refno].push_back(Fdzero);}
    }
  }

  // Loop over all images
  imgno=0;
  SF.go_beginning();
  while ((!SF.eof())) {

    fn_img=SF.NextImg();
    fn_trans=fn_img.remove_directories(offsets_keepdir);
    fn_trans=fn_root+"_offsets/"+fn_trans+".off";

    img.read(fn_img,false,false,false,false);
    img().set_Xmipp_origin();
    if (fn_doc!="") {
      trans(0)=(double)ROUND(imgs_oldxoff[imgno]);
      trans(1)=(double)ROUND(imgs_oldyoff[imgno]);
      img().self_translate(trans,true);
    } 

    // Read optimal offsets for all references from disc
    if (fast_mode || limit_trans) {
      if (!read_offsets(fn_trans,allref_offsets)) {
        int itot=n_ref*2;
        if (do_mirror) itot*=2;
        allref_offsets.clear();
        allref_offsets.resize(itot);
        if (zero_offsets) for (int i=0; i<itot; i++) allref_offsets[i]=0.;
        else for (int i=0; i<itot; i++) allref_offsets[i]=-999.;
      }
    }

    // Read optimal orientations from memory
    if (limit_rot) {
      old_phi=imgs_oldphi[imgno];
      old_theta=imgs_oldtheta[imgno];
    }

    // For limited orientational search: preselect relevant directions
    preselect_directions(old_phi,old_theta,pdf_directions);

    if (maxCC_rather_than_ML) {
      // A. Use a maximum cross-correlation target function
      
      maxCC_search_complete(img(),Fref,Mref,search_shift,Msum_imgs,sumw,sumw_mirror,
			    maxcorr,opt_refno,opt_psi,opt_offsets,pdf_directions);
      
    } else if (limit_trans) {
      // B. Use a maximum-likelihood target function in real space
      //    with limited translational searches

      ML_integrate_locally(img(),Mref,Msum_imgs,wsum_sigma_noise,wsum_sigma_offset,
			   sumw,sumw_mirror,LL,maxcorr,opt_refno,opt_psi,
			   opt_offsets,allref_offsets,pdf_directions);
      
    } else {
      // C. Use a maximum-likelihood target function in real space
      //    with complete or reduced-space translational searches (-fast)
      
      if (fast_mode) preselect_significant_model_phi(img(),allref_offsets,Mref,
						     Msignificant,pdf_directions);
      else Msignificant.init_constant(1);
      ML_integrate_complete(img(),Fref,Msignificant,
			    Fwsum_imgs,wsum_sigma_noise,wsum_sigma_offset,sumw,sumw_mirror, 
			    LL,maxcorr,opt_refno,opt_psi,opt_offsets,allref_offsets,pdf_directions);

    }
    
    // Write optimal offsets for all references to disc
    if (fast_mode || limit_trans) {
      write_offsets(fn_trans,allref_offsets);
    }

    // Store optimal phi and theta in memory
    if (limit_rot) {
      imgs_oldphi[imgno]=Iref[opt_refno].Phi();
      imgs_oldtheta[imgno]=Iref[opt_refno].Theta();
    }

    // Output docfile 
    sumcorr+=maxcorr;
    if (write_docfile) {
      opt_flip=0.;
      if (-opt_psi>360.) { 
	opt_psi+=360.;
	opt_flip=1.;
      }
      dataline(0)=Iref[opt_refno].Phi();   // rot
      dataline(1)=Iref[opt_refno].Theta(); // tilt
      dataline(2)=opt_psi+360.;            // psi
      dataline(3)=trans(0)+opt_offsets(0); // Xoff
      dataline(4)=trans(1)+opt_offsets(1); // Yoff
      dataline(5)=(double)(opt_refno+1);   // Ref
      dataline(6)=opt_flip;                // Mirror 
      dataline(7)=maxcorr;                 // P_max/P_tot or Corr      
      DFo.append_comment(img.name());
      DFo.append_data_line(dataline);

    }

    if (verb>0) if (imgno%c==0) progress_bar(imgno);
    imgno++;
  }
  if (verb>0) progress_bar(nn);

  // reverse rotation of the weighted sums 
  if (maxCC_rather_than_ML || limit_trans) 
    reverse_rotate_reference(Fwsum_imgs,Msum_imgs,true,wsum_Mref);
  else 
    reverse_rotate_reference(Fwsum_imgs,Msum_imgs,false,wsum_Mref);

}

// Update all model parameters
void Prog_MLalign2D_prm::update_parameters(vector<matrix2D<double> > &wsum_Mref,
					   double &wsum_sigma_noise, double &wsum_sigma_offset, 
					   vector<double> &sumw, vector<double> &sumw_mirror, 
					   double &sumcorr, double &sumw_allrefs) {

  matrix1D<int> center(2),radial_count;
  matrix2D<complex<double> > Faux, Faux2;
  matrix2D<double> Maux;
  vector<matrix1D<double> > ssnr;
  FileName fn_tmp;
  double rr,thresh,aux;
  int c;

  // Pre-calculate sumw_allrefs & average Pmax/sumP or cross-correlation
  sumw_allrefs=0.;
  FOR_ALL_MODELS() { sumw_allrefs+=sumw[refno]; }
  sumcorr/=sumw_allrefs;

  // Update the reference images
  FOR_ALL_MODELS() {
    if (sumw[refno]>0.) {
      Iref[refno]()=wsum_Mref[refno];
      Iref[refno]()/=sumw[refno];
      Iref[refno].weight()=sumw[refno];
    } else {
      Iref[refno].weight()=0.;
      Iref[refno]().init_zeros(dim,dim);
    }
  }

  // Update the model fractions
  if (!fix_fractions) {
    FOR_ALL_MODELS() {
      if (sumw[refno]>0.) {
	alpha_k[refno]=sumw[refno]/sumw_allrefs;
	mirror_fraction[refno]=sumw_mirror[refno]/sumw[refno];
      } else {
	alpha_k[refno]=0.;
	mirror_fraction[refno]=0.;
      }
    }
  }

  // Update sigma of the origin offsets
  if (!fix_sigma_offset) sigma_offset=sqrt(wsum_sigma_offset/(2*sumw_allrefs));

  // Update the noise parameters
  if (!fix_sigma_noise) sigma_noise=sqrt(wsum_sigma_noise/(sumw_allrefs*dim*dim));

}

// Check convergence 
bool Prog_MLalign2D_prm::check_convergence(vector<double> &conv) {

  bool converged=true;
  double convv;
  matrix2D<double> Maux;

  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  
  conv.clear();
  FOR_ALL_MODELS() {
    if (Iref[refno].weight()>0.) {
      Maux=mul_elements(Iold[refno](),Iold[refno]());
      convv=1./(Maux.compute_avg());
      Maux=Iold[refno]()-Iref[refno]();
      Maux=mul_elements(Maux,Maux);
      convv*=Maux.compute_avg();
      conv.push_back(convv);
      if (convv>eps) converged=false;
    } else {
      conv.push_back(-1.);
    }
  }

  return converged;
} 

// Output to screen
void Prog_MLalign2D_prm::output_to_screen(int &iter, double &sumcorr, double &LL) {
  if (verb>0) { 
    if (maxCC_rather_than_ML) cout <<"  iter "<<iter<<" <CC>= "+FtoA(sumcorr,10,5);
    else {
      cout <<"  iter "<<iter<<" noise= "<<FtoA(sigma_noise,10,7)<<" offset= "<<FtoA(sigma_offset,10,7);
      cout <<"  LL= "<<LL<<" <Pmax/sumP>= "<<sumcorr<<endl;
      cout <<"  Model  fraction  mirror-fraction "<<endl;
      FOR_ALL_MODELS() {
	cout <<"  "<<ItoA(refno+1,5)<<" "<<FtoA(alpha_k[refno],10,7)<<" "<<FtoA(mirror_fraction[refno],10,7)<<endl;
      }
    }
  }

}

void Prog_MLalign2D_prm::write_output_files(const int iter, DocFile &DFo, 
					    double &sumw_allrefs, double &LL, double &avecorr, 
					    vector<double> &conv) {

  FileName          fn_tmp,fn_base,fn_tmp2;
  matrix1D<double>  fracline(3);
  SelFile           SFo,SFc;
  DocFile           DFl;
  string            comment;
 
  DFl.clear();
  SFo.clear();
  SFc.clear();

  fn_base=fn_root;
  if (iter>=0) {
    fn_base+="_it";
    fn_base.compose(fn_base,iter,"");
  }

  // Write out current reference images and fill sel & log-file
  FOR_ALL_MODELS() {
    fn_tmp=fn_base+"_ref";
    fn_tmp.compose(fn_tmp,refno+1,"");
    fn_tmp=fn_tmp+".xmp";
    Iref[refno].write(fn_tmp);
    SFo.insert(fn_tmp,SelLine::ACTIVE);
    fracline(0)=alpha_k[refno];
    fracline(1)=mirror_fraction[refno];
    fracline(2)=1000*conv[refno]; // Output 1000x the change for precision
    DFl.insert_comment(fn_tmp);
    DFl.insert_data_line(fracline);
  }

  // Write out sel & log-file
  fn_tmp=fn_base+".sel";
  SFo.write(fn_tmp);

  DFl.go_beginning();
  comment="MLalign2D-logfile: Number of images= "+FtoA(sumw_allrefs);
  if (maxCC_rather_than_ML) comment+=" <CC>= "+FtoA(avecorr,10,5);
  else { 
    comment+=" LL= "+FtoA(LL,10,5)+" <Pmax/sumP>= "+FtoA(avecorr,10,5);
    DFl.insert_comment(comment);
    comment="-noise "+FtoA(sigma_noise,10,7)+" -offset "+FtoA(sigma_offset,10,7)+" -istart "+ItoA(iter+1);
  }
  DFl.insert_comment(comment);
  DFl.insert_comment(cline);
  DFl.insert_comment("columns: model fraction (1); mirror fraction (2); 1000x signal change (3)");
  fn_tmp=fn_base+".log";
  DFl.write(fn_tmp);

  if (write_docfile) {
    // Write out docfile with optimal transformation & references
    fn_tmp=fn_base+".doc";
    DFo.write(fn_tmp);
  }

  if (write_selfiles) {
    // Also write out selfiles of all experimental images, 
    // classified according to optimal reference image
    for (int refno=0;refno<n_ref; refno++) { 
      DFo.go_beginning();
      SFo.clear();
      for (int n=0; n<DFo.dataLineNo(); n++ ) {
	DFo.next();
	fn_tmp=((DFo.get_current_line()).get_text()).erase(0,3);
	DFo.adjust_to_data_line();
	if ((refno+1)==(int)DFo(5)) SFo.insert(fn_tmp,SelLine::ACTIVE);
      }
      fn_tmp=fn_root+"_ref";
      fn_tmp.compose(fn_tmp,refno+1,"sel");
      SFo.write(fn_tmp);
    }
  }


}

