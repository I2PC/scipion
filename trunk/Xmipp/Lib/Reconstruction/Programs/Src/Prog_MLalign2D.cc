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
void Prog_MLalign2D_prm::read(int argc, char **argv) _THROW  {


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
      fn_sel=DFi.name();
      fn_sel=fn_sel.without_extension()+"_restart.sel";
      comment=" -frac "+DFi.name()+" -ref "+fn_sel;
      comment+=(DFi.get_current_line()).get_text();
      DFi.next();
      cline=(DFi.get_current_line()).get_text();
      comment=comment+cline;
      generate_command_line(comment,argc,argv,copy);
      // Read images names from restart file
      DFi.next(); //
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
  } else {
    for (int i=1; i<argc; i++) {
      cline=cline+(string)argv[i]+" ";
    }
  }
    
  // Read command line
  n_ref=AtoI(get_param(argc,argv,"-nref","0"));
  fn_ref=get_param(argc,argv,"-ref","");
  SF.read(get_param(argc,argv,"-i"));
  fn_root=get_param(argc,argv,"-o","MLalign2D");
  LSQ_rather_than_ML=check_param(argc,argv,"-LSQ");
  psi_step=AtoF(get_param(argc,argv,"-psi_step","5"));
  Niter=AtoI(get_param(argc,argv,"-iter","100"));
  sigma_noise=AtoF(get_param(argc,argv,"-noise","1"));
  sigma_offset=AtoF(get_param(argc,argv,"-offset","3"));
  do_mirror=check_param(argc,argv,"-mirror");
  eps=AtoF(get_param(argc,argv,"-eps","5e-5"));
  fn_frac=get_param(argc,argv,"-frac","");
  write_docfile=check_param(argc,argv,"-output_docfile");
  write_selfiles=check_param(argc,argv,"-output_selfiles");
  write_intermediate=!check_param(argc,argv,"-dont_output_intermediate");
  apply_shifts=check_param(argc,argv,"-apply_shifts");
  fix_fractions=check_param(argc,argv,"-fix_fractions");
  fix_sigma_offset=check_param(argc,argv,"-fix_sigma_offset");
  fix_sigma_noise=check_param(argc,argv,"-fix_sigma_noise");
  verb=AtoI(get_param(argc,argv,"-verb","1"));
  // Hidden arguments
  istart=AtoI(get_param(argc,argv,"-istart","1"));
  do_esthetics=check_param(argc,argv,"-esthetics");
  fast_mode=check_param(argc,argv,"-fast");
  do_precenter=check_param(argc,argv,"-precenter");
  max_shift=AtoF(get_param(argc,argv,"-max_shift","-1"));

  // Fill all memory vectors with appropriate variables
  read_all_input_in_memory();

}

// Read all input selfiles in memory 
void Prog_MLalign2D_prm::read_all_input_in_memory() _THROW {

  FileName fn_img;
  ImageXmipp img;
  matrix1D<double> offsets(2);
  matrix2D<double> A(3,3);
  int refno=0;

  // Set nr_psi
  nr_psi=ROUND(90./psi_step);
  psi_step=90./nr_psi;

  // Create references from random subset averages, or read them from selfile
  if (fn_ref=="") {
    if (n_ref!=0) {
      generate_initial_references();
    } else {
      REPORT_ERROR(1,"Please provide -ref or -nref");
    }
  }
 
  // Read image- and reference- selfiles
  if (Is_ImageXmipp(fn_ref)) {
    SFr.reserve(1);
    SFr.insert(fn_ref);
  } else {
    SFr.read(fn_ref);
  }

  // Construct flipping (0, 90, 180 & 270 degree rotation) matrices
  if (do_mirror) nr_flip=8;
  else nr_flip=4;
  A.init_identity();
  F.push_back(A);
  A(0,0)=0.; A(1,1)=0.; A(1,0)=1.; A(0,1)=-1;
  F.push_back(A);
  A(0,0)=-1.; A(1,1)=-1.; A(1,0)=0.; A(0,1)=0;
  F.push_back(A);
  A(0,0)=0.; A(1,1)=0.; A(1,0)=-1.; A(0,1)=1;
  F.push_back(A);
  if (do_mirror) {
    A.init_identity(); A(0,0)=-1;
    F.push_back(A);
    A(0,0)=0.; A(1,1)=0.; A(1,0)=1.; A(0,1)=1;
    F.push_back(A);
    A(0,0)=1.; A(1,1)=-1.; A(1,0)=0.; A(0,1)=0;
    F.push_back(A);
    A(0,0)=0.; A(1,1)=0.; A(1,0)=-1.; A(0,1)=-1;
    F.push_back(A);
  }

  SF.ImgSize(dim,dim);

  // Read in all reference images in memory
  n_ref=0;
  SFr.go_beginning();
  while ((!SFr.eof())) {
    img.read(SFr.NextImg());
    img().set_Xmipp_origin();
    if (apply_shifts) {
      offsets(0)=ROUND(img.Xoff());
      offsets(1)=ROUND(img.Yoff());
      img().self_translate(offsets,WRAP);
    }
    Iref.push_back(img);
    Iold.push_back(img);
    // Default start is all equal model fractions
    alpha_k.push_back((double)1/SFr.ImgNo());
    // Default start is half-half mirrored images
    if (do_mirror) mirror_fraction.push_back(0.5);
    else mirror_fraction.push_back(0.);
    n_ref++;
    refno++;
  }

  if (fn_frac!="") {
    DocFile  DF;
    DocLine DL;
    DF.read(fn_frac);
    DF.go_first_data_line();
    double sumfrac=0.;
    for (refno=0; refno<n_ref; refno++) {
      DL=DF.get_current_line();
      alpha_k[refno]=DL[0];
      if (do_mirror) {
	if (DL[1]>1.||DL[1]<0.) REPORT_ERROR(1,"Mirror fraction (2nd column) should be [0,1]!");
	mirror_fraction[refno]=DL[1];
      }
      sumfrac+=alpha_k[refno];
      DF.next_data_line();
    }
    if (ABS(sumfrac-1.)>1e-4) 
      cerr << " ->WARNING: Sum of all expected model fractions ("<<FtoA(sumfrac)<<") is not one!"<<endl;
    for (refno=0; refno<n_ref; refno++) { alpha_k[refno]/=sumfrac; }
  }

  if (max_shift>0) {
    if (!LSQ_rather_than_ML) REPORT_ERROR(1,"-max_shift is only for -LSQ!!!");
    sigma_noise=1;
    fix_sigma_noise=true;
    fix_sigma_offset=true;
    fix_fractions=true;
  }

  if (fast_mode) {
    while (!SF.eof()) {
      img.read(SF.NextImg());
      offset_x.push_back(0.);
      offset_y.push_back(0.);
    }
  }
}

// Generate initial references =============================================
void Prog_MLalign2D_prm::generate_initial_references() _THROW  {

  SelFile SFtmp, SFout;
  ImageXmipp Iave,Itmp;
  double dummy;
  FileName fn_tmp;

  // Make random subsets and calculate average images
  cerr << " Generating initial references by averaging over random subsets" <<endl;
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
  }
  fn_ref=fn_root+"_it";
  fn_ref.compose(fn_ref,0,"sel");
  SFr.write(fn_ref);

}

// Show ====================================================================
void Prog_MLalign2D_prm::show() {

  if (verb>0) {
    // To screen
    cerr << " Input images            : "<< SF.name()<<" ("<<SF.ImgNo()<<")"<<endl;
    cerr << " Reference images        : "<< fn_ref<<" ("<<SFr.ImgNo()<<")"<<endl;
    cerr << " Output rootname         : "<< fn_root<<endl;
    cerr << " Number of iterations    : "<< Niter<<endl;
    cerr << " Stopping criterium      : "<< eps<<endl;
    cerr << " initial sigma noise     : "<< sigma_noise<<endl;
    cerr << " initial sigma offset    : "<< sigma_offset<<endl;
    cerr << " Psi sampling interval   : "<< psi_step<<endl;
    if (fn_frac!="") {
      cerr << " -> Read initial model fractions from "<< fn_frac<<endl;
    }
    if (do_mirror) {
      cerr << " -> Check mirror image of each reference as well."<<endl;
    }
    if (LSQ_rather_than_ML) {
      cerr << " -> Use least-squares instead of maximum likelihood target."<<endl;
    }
    if (write_docfile) {
      cerr << " -> Write out docfile with most likely angles & translations. "<<endl;
    }
    if (write_selfiles) {
      cerr << " -> Write out selfiles with most likely reference assignments. "<<endl;
    }

    // Hidden stuff
    if (!write_intermediate) {
      cerr << " -> Do not write out images after each iteration."<<endl;
    }
    if (apply_shifts) {
      cerr << " -> Apply shifts stored in header of 2D-images."<<endl;
    }
    if (fix_fractions) {
      cerr << " -> Do not update estimates of model fractions."<<endl;
    }
    if (fix_sigma_offset) {
      cerr << " -> Do not update sigma-estimate of origin offsets."<<endl;
    }
    if (fix_sigma_noise) {
      cerr << " -> Do not update sigma-estimate of noise."<<endl;
    }
    if (do_esthetics) {
      cerr << " -> Perform esthetics on (0,0)-pixel artifacts"<<endl;
    }
    if (max_shift>0) {
      cerr << " -> Use maximum shift criterion instead of pdf; max_shift = "<<max_shift<<endl;
    }
    cerr << " -----------------------------------------------------------------"<<endl;
  }

} 

// Usage ===================================================================
void Prog_MLalign2D_prm::usage() {
  cerr << "Usage:  MLalign2D [options] "<<endl;
  cerr << "   -i <selfile>                : Selfile with input images \n"
       << "   -ref <selfile/image>        : Selfile with initial reference images/single reference image \n"
       << "      OR -nref <int>               OR number of bias-free references to generate automatically\n"
       << " [ -o <rootname=\"MLalign2D\"> ] : Output rootname \n"
       << " [ -noise <float=1> ]          : Expected standard deviation for pixel noise \n"
       << " [ -offset <float=3> ]         : Expected standard deviation for origin offset [pix]\n"
       << " [ -eps <float=5e-5> ]         : Stopping criterium \n"
       << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n"
       << " [ -psi_step <float=5> ]       : In-plane rotation sampling interval [deg]\n"
       << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n"
       << " [ -mirror ]                   : Also check mirror image of each reference \n"
       << " [ -LSQ ]                      : Use least-squares instead of maximum likelihood target \n"
       << " [ -output_docfile ]           : Write out docfile with most likeliy angles & translations \n"
       << " [ -output_selfiles ]          : Write out selfiles with most likely reference assignments \n"
       << endl;
}

// Calculate probability density function of all in-plane transformations phi
void Prog_MLalign2D_prm::calculate_pdf_phi(double &sigma_offset, matrix2D<double> &P_phi, matrix2D<double> &Mr2) _THROW {

  double r2,pdfpix,sum;
  P_phi.resize(dim,dim);
  P_phi.set_Xmipp_origin();
  Mr2.resize(dim,dim);
  Mr2.set_Xmipp_origin();

  FOR_ALL_ELEMENTS_IN_MATRIX2D(P_phi) {
    r2=(double)(j*j + i*i);
    if (sigma_offset>0.) {
      pdfpix=exp(-r2/(2*sigma_offset*sigma_offset));
      pdfpix/=2*PI*sigma_offset*sigma_offset*nr_psi*4;
    } else {
      if (j==0 && i==0) pdfpix=1.;
      else pdfpix=0.;
    }
    MAT_ELEM(P_phi,i,j)=pdfpix;
    MAT_ELEM(Mr2,i,j)=(float)r2;
  }

  if (max_shift>0.) {
    matrix2D<int> shiftmask;
    shiftmask.resize(dim,dim);
    shiftmask.set_Xmipp_origin();
    BinaryCircularMask(shiftmask,max_shift,INNER_MASK);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(P_phi) {
      MAT_ELEM(P_phi,i,j)=(double)MAT_ELEM(shiftmask,i,j);
    }  
  }

}

// Rotate reference for all models and rotations and fill Fref vectors =============
void Prog_MLalign2D_prm::rotate_reference(vector<ImageXmipp> &Iref, bool &real_space,
					       vector <vector< matrix2D<double> > > &Mref,
					       vector <vector< matrix2D<complex<double> > > > &Fref) _THROW {

  double psi,dum,avg,sum;
  matrix2D<double> Maux;
  matrix2D<complex<double> > Faux;
  vector<matrix2D<complex <double> > > dumF;
  vector<matrix2D<double> > dumM;
  matrix2D<int> mask, omask;

  Maux.init_zeros(dim,dim);
  Maux.set_Xmipp_origin();
  Faux.set_Xmipp_origin();
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  omask.resize(dim,dim);
  omask.set_Xmipp_origin();
  BinaryCircularMask(mask,dim/2,INNER_MASK);
  BinaryCircularMask(omask,dim/2,OUTSIDE_MASK);

  Fref.clear();
  Mref.clear();

  FOR_ALL_MODELS() {
    if (real_space) Mref.push_back(dumM);
    else Fref.push_back(dumF);
    compute_stats_within_binary_mask(omask,Iref[refno](),dum,dum,avg,dum);
    FOR_ALL_ROTATIONS() {

      // Add arbitrary number to avoid 0-degree rotation (lacking interpolation effects)
      psi=(double)(ipsi*90./nr_psi)+1.75; 
      Maux=Iref[refno]().rotate(psi,DONT_WRAP);
      apply_binary_mask(mask,Maux,Maux,avg);

      // Normalize the magnitude of the rotated references to 1st rot of 1st ref
      // This is necessary because interpolation due to rotation can lead to lower overall Fref 
      // This would result in lower probabilities for those rotations
      if (ipsi==0 && refno==0) A2=Maux.sum2();
      sum=Maux.sum2();

      if (real_space) {
	Mref[refno].push_back(Maux);
	Mref[refno][ipsi]*=sqrt(A2/sum);
      } else {
	FourierTransform(Maux,Faux);
	Faux*=dim*dim;
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux) {
	  dMij(Faux,i,j)=conj(dMij(Faux,i,j));
	}
	Fref[refno].push_back(Faux);
	Fref[refno][ipsi]*=sqrt(A2/sum);
      }
    }
  }

}


// Collect all rotations and sum to update Iref() for all models ==========
void Prog_MLalign2D_prm::reverse_rotate_reference(vector <vector< matrix2D<double> > > &Mnew, 
						       vector <vector< matrix2D<complex<double> > > > &Fnew, 
						       bool &real_space, 
						       vector<matrix2D<double> > &Mref) _THROW {

  double psi,dum,avg,ang;
  int s;
  matrix2D<double> Maux,Maux2;
  matrix2D<int> mask, omask;
  Maux.resize(dim,dim);
  Maux2.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Maux2.set_Xmipp_origin();

  Mref.clear();
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  omask.resize(dim,dim);
  omask.set_Xmipp_origin();
  BinaryCircularMask(mask,dim/2,INNER_MASK);
  BinaryCircularMask(omask,dim/2,OUTSIDE_MASK);

  FOR_ALL_MODELS() {
    Maux.init_zeros();
    Mref.push_back(Maux);
    FOR_ALL_ROTATIONS() {
      // Add arbitrary number to avoid 0-degree rotation without interpolation effects
      psi=(double)(ipsi*90./nr_psi)+1.75;
      if (real_space) {
	compute_stats_within_binary_mask(omask,Mnew[refno][ipsi],dum,dum,avg,dum);
	Maux=Mnew[refno][ipsi].rotate(-psi,DONT_WRAP);
      } else {
	InverseFourierTransform(Fnew[refno][ipsi],Maux2);
	Maux2/=dim*dim;
	CenterFFT(Maux2,true);
	compute_stats_within_binary_mask(omask,Maux2,dum,dum,avg,dum);
	Maux=Maux2.rotate(-psi,DONT_WRAP);
      }
      apply_binary_mask(mask,Maux,Maux,avg);
      Mref[refno]+=Maux;
    }
  }

}

// For fast_mode: pre-center images before entering ML-refinement
void Prog_MLalign2D_prm::precenter_images() _THROW{

  int ioptx,iopty,imgno=0;
  ImageXmipp img;
  matrix2D<double> Maux,Maveref;

  Maveref.resize(dim,dim);
  Maveref.set_Xmipp_origin();
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();

  FOR_ALL_MODELS() { Maveref+=Iref[refno](); }
  Maveref/=n_ref;

  cerr << "  Pre-centering images against average of all reference images .. "<<endl;
  SF.go_beginning();
  while (!SF.eof()) {
    img.read(SF.NextImg());
    img().set_Xmipp_origin();
    correlation_matrix(img(),Maveref,Maux);
    Maux.max_index(iopty,ioptx);
    offset_x[imgno]=-(double)ioptx;
    offset_y[imgno]=-(double)iopty;
    imgno++;
  }

}



// Maximum Likelihood calculation for one image ============================================
// Integration over all models and in-plane rotations
void Prog_MLalign2D_prm::ML_integrate_model_phi(matrix2D<double> &Mimg, vector <vector< matrix2D<double > > > &Mref, 
						vector<double> &P_model, vector<vector<matrix2D<double> > > &Mwsum_imgs,
						double &wsum_sigma_noise, vector<double> &sumw, vector<double> &sumw_mirror, 
						double &LL, double &maxcorr, 
						int &opt_refno, double &opt_psi, double &opt_xoff, double &opt_yoff) _THROW {

  matrix2D<double> Maux;
  vector<vector<double> > Vweight;
  vector<matrix2D<double> > Mimg_flip;
  vector<double >  dum2;
  vector<double> refw, refw_mirror;
  double sigma_noise2,XiA,Xi2,aux,fracpdf;
  double wsum_corr=0., sum_refw=0.;
  double CC,maxCC=-99.e99,maxweight=-99.e99;
  int irot,opt_ipsi,opt_iflip,ioptx,iopty;

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

  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  sigma_noise2=sigma_noise*sigma_noise;
  Xi2=Mimg.sum2();
  Vweight.clear();
  Mimg_flip.clear();
  refw.clear();
  refw_mirror.clear();

  // Flip images and calculate correlations and maximum correlation
  FOR_ALL_FLIPS() {
    apply_geom(Maux,F[iflip],Mimg,IS_INV,WRAP);
    Mimg_flip.push_back(Maux);
    FOR_ALL_MODELS() {
      Vweight.push_back(dum2);
      refw.push_back(0.);
      refw_mirror.push_back(0.);
      FOR_ALL_ROTATIONS() {
	irot=iflip*nr_psi+ipsi;
	CC=0.;
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Maux) {
	  CC+=dMij(Maux,i,j)*dMij(Mref[refno][ipsi],i,j);
	}
	Vweight[refno].push_back(CC);
	if (CC>maxCC) maxCC=CC;
      }
    }
  }

  // Now that we have maxCC calculate the weighting matrices and maxweight
  FOR_ALL_MODELS() {
    FOR_ALL_FLIPS() {
      FOR_ALL_ROTATIONS() {
	irot=iflip*nr_psi+ipsi;
	if (iflip<4) fracpdf=P_model[refno]*(1.-mirror_fraction[refno]);
	else fracpdf=P_model[refno]*mirror_fraction[refno];
	aux=Vweight[refno][irot]-maxCC;
	aux=exp(aux/sigma_noise2)*fracpdf;
	wsum_corr+=aux*Vweight[refno][irot];
	if (aux>maxweight) {
	  maxweight=aux;
	  opt_refno=refno;
	  opt_ipsi=ipsi;
	  opt_iflip=iflip;
	  maxcorr=Vweight[refno][irot];
	}
	Vweight[refno][irot]=aux;
      }
    }
  }

  // Now that we have maxweight calculate which weighting matrices are significant
  // For LSQ, convert weights in a delta-fucntion at the maximum
  FOR_ALL_MODELS() {
    FOR_ALL_ROTATIONS() {
      FOR_ALL_FLIPS() {
	irot=iflip*nr_psi+ipsi;
	if (LSQ_rather_than_ML) {
	  if (refno==opt_refno && opt_ipsi==ipsi && opt_iflip==iflip) {
	    Vweight[refno][irot]=1.;
	    if (iflip<4) refw[refno]+=1;
	    else refw_mirror[refno]+=1;
	  } else Vweight[refno][irot]=0.;
	} else {
	  if (iflip<4) refw[refno]+=Vweight[refno][irot];
	  else refw_mirror[refno]+=Vweight[refno][irot];
	}
      }
    }
    sum_refw+=refw[refno]+refw_mirror[refno];
  }

  // In ML-case: write out maxweight/sumweight instead of maxcorr
  if (!LSQ_rather_than_ML) maxcorr=maxweight/sum_refw; 

  // Normalize all weighted sums by sum_refw such that sum over all weights is one!
  // And accumulate the FT of the weighted, shifted images.
  if (LSQ_rather_than_ML) wsum_sigma_noise+=A2+Xi2-2*maxcorr;
  else   wsum_sigma_noise+=A2+Xi2-2*wsum_corr/sum_refw;

  FOR_ALL_MODELS() {
    sumw[refno]+=(refw[refno]+refw_mirror[refno])/sum_refw;
    sumw_mirror[refno]+=refw_mirror[refno]/sum_refw;
    FOR_ALL_ROTATIONS() {
      FOR_ALL_FLIPS() {
	irot=iflip*nr_psi+ipsi;
	Vweight[refno][irot]/=sum_refw;
	Mwsum_imgs[refno][ipsi]+=(Vweight[refno][irot]*Mimg_flip[iflip]);
      }
    }
  }

  // Compute Log Likelihood
  // 1st term: log(refw_i)
  // 2nd term: for subtracting maxc
  // 3rd term: for only considering Xi*A instead of (A-Xi)^2
  // 4th term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
  LL+= log(sum_refw) + maxCC/sigma_noise2 - (A2+Xi2)/(2*sigma_noise2) - dim*dim*log(2.50663*sigma_noise);


  // Now for optimal model and phi get the origin offsets
  correlation_matrix(Mimg_flip[opt_iflip],Mref[opt_refno][opt_ipsi],Maux);
  Maux.max_index(iopty,ioptx);
  opt_xoff=-(double)ioptx*DIRECT_MAT_ELEM(F[opt_iflip],0,0)-(double)iopty*DIRECT_MAT_ELEM(F[opt_iflip],0,1);
  opt_yoff=-(double)ioptx*DIRECT_MAT_ELEM(F[opt_iflip],1,0)-(double)iopty*DIRECT_MAT_ELEM(F[opt_iflip],1,1);
  opt_psi=-psi_step*(opt_iflip*nr_psi+opt_ipsi);


}




// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void Prog_MLalign2D_prm::ML_integrate_model_phi_trans(matrix2D<double> &Mimg, vector <vector< matrix2D<complex<double> > > > &Fref, 
			vector<double> &P_model, matrix2D<double> &P_phi, matrix2D<double> &Mr2,
                        vector <vector< matrix2D<complex<double> > > > &Fwsum_imgs, 
			double &wsum_sigma_noise, double &wsum_sigma_offset, vector<double> &sumw, vector<double> &sumw_mirror, 
			double &LL, double &maxcorr, 
                        int &opt_refno, double &opt_psi, double &opt_xoff, double &opt_yoff) _THROW {

  matrix2D<double> Maux,Mdzero;
  matrix2D<complex<double> > Fimg, Faux;
  vector<matrix2D<complex<double> > > Fimg_flip;
  vector<matrix2D<double> > dumM;
  vector<bool> dumb;
  vector <vector< matrix2D<double> > > Mweight;
  vector<vector<bool> >  significant_weight;
  vector<double> refw, refw_mirror;
  double sigma_noise2,XiA,Xi2,aux,fracpdf,sigw,add,max,maxc=-99.e99;
  double wsum_corr=0., sum_refw=0.;
  double testweight,maxweight=-99.e99;
  int irot,sigdim,ioptx,iopty,ioptpsi,ioptflip;
  int itmp1,itmp2,itmp3,itmp4;

  opt_refno=-1;
  ioptpsi=0;
  ioptx=0;
  iopty=0;
  ioptflip=0;
    

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

  // Only translations smaller than 6 sigma_offset are considered!
  // This saves a lot of memory and CPU! (typically a factor 2, depending on sigma_offset vs. dim)
  sigdim=2*CEIL(sigma_offset*6);
  sigdim++; // (to get uneven number)
  sigdim=MIN(dim,sigdim);

  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Faux.resize(dim,dim);
  Faux.set_Xmipp_origin();
  refw.clear();
  refw_mirror.clear();
  Fimg_flip.clear();
  Mweight.clear();
  sigma_noise2=sigma_noise*sigma_noise;
  Xi2=Mimg.sum2();

  // Flip images and calculate correlation matrices and maximum correlation
  FOR_ALL_FLIPS() {
    apply_geom(Maux,F[iflip],Mimg,IS_INV,WRAP);
    FourierTransform(Maux,Fimg);
    Fimg*=dim*dim;
    Fimg_flip.push_back(Fimg);
    FOR_ALL_MODELS() {
      Mweight.push_back(dumM); 
      significant_weight.push_back(dumb);
      refw.push_back(0.);
      refw_mirror.push_back(0.);
      FOR_ALL_ROTATIONS() {
	irot=iflip*nr_psi+ipsi;
	mul_elements(Fimg,Fref[refno][ipsi],Faux);
	InverseFourierTransform(Faux,Maux);
	Maux/=dim*dim;
	CenterFFT(Maux,true);
	significant_weight[refno].push_back(true);
	Mweight[refno].push_back(Mdzero);
	Mweight[refno][irot].resize(sigdim,sigdim);
	Mweight[refno][irot].set_Xmipp_origin();
	FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot]) {
	  MAT_ELEM(Mweight[refno][irot],i,j)=MAT_ELEM(Maux,i,j);
	}
	max=Mweight[refno][irot].compute_max();
	if (max>maxc) maxc=max;
      }
    }
  }

  // Now that we have maxc calculate the weighting matrices and maxweight
  FOR_ALL_MODELS() {
    FOR_ALL_ROTATIONS() {
      FOR_ALL_FLIPS() {
	irot=iflip*nr_psi+ipsi;
	if (iflip<4) fracpdf=P_model[refno]*(1.-mirror_fraction[refno]);
	else fracpdf=P_model[refno]*mirror_fraction[refno];
	FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot]) {
	  aux=MAT_ELEM(Mweight[refno][irot],i,j)-maxc;
	  aux=exp(aux/sigma_noise2)*fracpdf*MAT_ELEM(P_phi,i,j);
	  wsum_corr+=aux*MAT_ELEM(Mweight[refno][irot],i,j);
	  if (aux>maxweight) {
	    maxweight=aux;
	    iopty=i; ioptx=j; 
	    ioptpsi=ipsi;
	    ioptflip=iflip;
	    opt_refno=refno;
	    maxcorr=MAT_ELEM(Mweight[refno][irot],i,j);
	  }
	  MAT_ELEM(Mweight[refno][irot],i,j)=aux;
	}
      }
    }
  }

  // Calculate optimal transformation parameters
  opt_xoff=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],0,0)-(double)iopty*DIRECT_MAT_ELEM(F[ioptflip],0,1);
  opt_yoff=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],1,0)-(double)iopty*DIRECT_MAT_ELEM(F[ioptflip],1,1);
  opt_psi=-psi_step*(ioptflip*nr_psi+ioptpsi);

  // Now that we have maxweight calculate which weighting matrices are significant,
  // For LSQ, convert weights in a delta-fucntion at the maximum
  FOR_ALL_MODELS() {
    FOR_ALL_ROTATIONS() {
      FOR_ALL_FLIPS() {
	irot=iflip*nr_psi+ipsi;
	if (LSQ_rather_than_ML) {
	  Mweight[refno][irot].init_zeros();
	  significant_weight[refno][irot]=false;
	  if (refno==opt_refno && ioptpsi==ipsi && ioptflip==iflip) { 
	    MAT_ELEM(Mweight[refno][irot],iopty,ioptx)=1;
	    significant_weight[refno][irot]=true;
	    if (iflip<4) refw[refno]+=1;
	    else refw_mirror[refno]+=1;
	  }
	} else {
	  if (Mweight[refno][irot].compute_max()>SIGNIFICANT_WEIGHT_LOW*maxweight) {
	    significant_weight[refno][irot]=true;
	    if (iflip<4) refw[refno]+=Mweight[refno][irot].sum();
	    else refw_mirror[refno]+=Mweight[refno][irot].sum();
	  } else significant_weight[refno][irot]=false;
	}
      }
    }
    sum_refw+=refw[refno]+refw_mirror[refno];
  }

  // In ML-case: write out maxweight/sumweight instead of maxcorr
  if (!LSQ_rather_than_ML) maxcorr=maxweight/sum_refw; 

  // Normalize all weighted sums by sum_refw such that sum over all weights is one!
  // And accumulate the FT of the weighted, shifted images.
  if (LSQ_rather_than_ML) wsum_sigma_noise+=A2+Xi2-2*maxcorr;
  else   wsum_sigma_noise+=A2+Xi2-2*wsum_corr/sum_refw;

  FOR_ALL_MODELS() {
    sumw[refno]+=(refw[refno]+refw_mirror[refno])/sum_refw;
    sumw_mirror[refno]+=refw_mirror[refno]/sum_refw;
    FOR_ALL_ROTATIONS() {
      FOR_ALL_FLIPS() {
	irot=iflip*nr_psi+ipsi;
	if (significant_weight[refno][irot]) {
	  Mweight[refno][irot]/=sum_refw;
	  // Use Maux, because Mweight is smaller than dim x dim!
	  Maux.init_zeros();
	  FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot]) {
	    MAT_ELEM(Maux,i,j)=MAT_ELEM(Mweight[refno][irot],i,j);
	    wsum_sigma_offset+=MAT_ELEM(Maux,i,j)*MAT_ELEM(Mr2,i,j);
	  }
	  FourierTransform(Maux,Faux);
	  Faux*=dim*dim;
	  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux) {
	    dMij(Fwsum_imgs[refno][ipsi],i,j)+=conj(dMij(Faux,i,j))*dMij(Fimg_flip[iflip],i,j);
	  }
	}
      }
    }
  }

  // Compute Log Likelihood
  // 1st term: log(refw_i)
  // 2nd term: for subtracting maxc
  // 3rd term: for only considering Xi*A instead of (A-Xi)^2
  // 4th term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
  LL+= log(sum_refw) + maxc/sigma_noise2 - (A2+Xi2)/(2*sigma_noise2) - dim*dim*log(2.50663*sigma_noise);


}


void Prog_MLalign2D_prm::ML_sum_over_all_images(SelFile &SF, vector<ImageXmipp> &Iref, 
			  matrix2D<double> &P_phi, matrix2D<double> &Mr2, vector<double> &P_model, vector<double> &mirror_fraction, 
			  double &LL, double &avecorr, DocFile &DFo, 
			  vector<matrix2D<double> > &wsum_Mref,
			  double &wsum_sigma_noise, double &wsum_sigma_offset, vector<double> &sumw, vector<double> &sumw_mirror) _THROW {


  ImageXmipp img;
  FileName fn_img;
  vector <vector< matrix2D<complex<double> > > > Fref, Fwsum_imgs;
  vector <vector< matrix2D<double> > > Mref,Mwsum_imgs;
  matrix2D<complex <double> > Fdzero;
  matrix2D<double> Mdzero;
  vector<matrix2D<complex <double> > > dum;
  vector<matrix2D<double > > dum2;
  matrix1D<double> offsets(2), dataline(8);  
  double opt_psi,opt_flip,opt_xoff,opt_yoff,maxcorr;
  int c,nn,imgno,opt_refno;

  double avex=0., avey=0.;

  // Generate (FT of) each rotated version of all references 
  rotate_reference(Iref,fast_mode,Mref,Fref);

  nn=SF.ImgNo();
  if (verb>0) init_progress_bar(nn);
  c=MAX(1,nn/60);

  Fwsum_imgs.clear();
  Mwsum_imgs.clear();
  sumw.clear();
  sumw_mirror.clear();

  Fdzero.resize(dim,dim);
  Fdzero.set_Xmipp_origin();
  Mdzero.resize(dim,dim);
  Mdzero.set_Xmipp_origin();
  LL=0.;
  wsum_sigma_noise=0.;
  wsum_sigma_offset=0.;
  avecorr=0.;
  FOR_ALL_MODELS() {
    if (fast_mode) Mwsum_imgs.push_back(dum2);
    else Fwsum_imgs.push_back(dum);
    sumw.push_back(0.);
    sumw_mirror.push_back(0.);
    FOR_ALL_ROTATIONS() {
      if (fast_mode) Mwsum_imgs[refno].push_back(Mdzero); 
      else Fwsum_imgs[refno].push_back(Fdzero); 
    }
  }
  imgno=0;
  SF.go_beginning();
  while ((!SF.eof())) {
    fn_img=SF.NextImg();
    img.read(fn_img);
    img().set_Xmipp_origin();
    offsets.init_zeros();
    if (apply_shifts) {
      offsets(0)=ROUND(img.Xoff());
      offsets(1)=ROUND(img.Yoff());
    }
    if (fast_mode) {
      offsets(0)+=offset_x[imgno];
      offsets(1)+=offset_y[imgno];      
    }
    if (apply_shifts || fast_mode) img().self_translate(offsets,WRAP);

    // Perform the integration over all in-plane transformations for a single image (against all refs)
    if (fast_mode) {
      ML_integrate_model_phi(img(),Mref,P_model,
				   Mwsum_imgs,wsum_sigma_noise, sumw,sumw_mirror, 
				   LL,maxcorr,opt_refno,opt_psi,opt_xoff,opt_yoff);
      avex+=ABS(opt_xoff);
      avey+=ABS(opt_yoff);

      offset_x[imgno]+=opt_xoff;
      offset_y[imgno]+=opt_yoff;
      opt_xoff=offset_x[imgno];
      opt_yoff=offset_y[imgno];

    } else {
      ML_integrate_model_phi_trans(img(),Fref,P_model,P_phi,Mr2,
				   Fwsum_imgs,wsum_sigma_noise, wsum_sigma_offset,sumw,sumw_mirror, 
				   LL,maxcorr,opt_refno,opt_psi,opt_xoff,opt_yoff);
    }
    avecorr+=maxcorr;
    if (write_docfile) {
      opt_flip=0.;
      if (apply_shifts) {
	opt_xoff+=ROUND(img.Xoff());
	opt_yoff+=ROUND(img.Yoff());
      }
      if (-opt_psi>360.) { 
	opt_psi+=360.;
	opt_flip=1.;
	opt_xoff=-opt_xoff;
      }
      dataline(0)=Iref[opt_refno].Phi();
      dataline(1)=Iref[opt_refno].Theta();
      dataline(2)=opt_psi+360.;
      dataline(3)=opt_xoff;
      dataline(4)=opt_yoff;
      dataline(5)=(double)(opt_refno+1);
      dataline(6)=opt_flip;
      dataline(7)=maxcorr;
      DFo.append_comment(img.name());
      DFo.append_data_line(dataline);
    }


    if (verb>0) if (imgno%c==0) progress_bar(imgno);
    //    progress_bar(imgno);
    imgno++;
  }
  if (verb>0) progress_bar(nn);
  avecorr/=imgno;

  if (fast_mode) {
    avex/=imgno; avey/=imgno;
    cerr <<"  Average residual origin offsets: X= "<<avex<<" Y="<<avey<<endl;
  }

  reverse_rotate_reference(Mwsum_imgs,Fwsum_imgs,fast_mode,wsum_Mref);

}

void Prog_MLalign2D_prm::write_output_files(const int iter, SelFile &SF, DocFile &DF, 
			 double &sumw_allrefs, double &LL, double &avecorr, vector<double> &conv ) _THROW {

  FileName fn_img,fn_tmp;
  matrix1D<double> fracline(3);
  string comment;
 
  DF.clear();
  if (iter<0) SF.clear();

  // Write out current reference images and fill sel & log-file
  FOR_ALL_MODELS() {
    fn_img=fn_root;
    if (iter>=0) {
      fn_img+="_it";
      fn_img.compose(fn_img,iter,"");
    }
    fn_img+="_ref";
    fn_img.compose(fn_img,refno+1,"");
    fn_img=fn_img+".xmp";
    Iref[refno].write(fn_img);
    SF.insert(fn_img,SelLine::ACTIVE);
    fracline(0)=alpha_k[refno];
    fracline(1)=mirror_fraction[refno];
    fracline(2)=1000*conv[refno]; // Output 1000x the change for precision
    DF.insert_comment(fn_img);
    DF.insert_data_line(fracline);
  }
  // Write out sel & log-file
  fn_tmp=fn_root;
  if (iter>=0) {
    fn_tmp+="_alliter";
  }
  fn_tmp+=".sel";
  SF.write(fn_tmp);

  DF.go_beginning();
  comment="MLalign2D-logfile: Number of images= "+FtoA(sumw_allrefs);
  if (!LSQ_rather_than_ML) comment+=" LL= "+FtoA(LL,10,5)+" R= "+FtoA(avecorr,10,5);
  DF.insert_comment(comment);
  comment="-noise "+FtoA(sigma_noise,10,7)+" -offset "+FtoA(sigma_offset,10,7)+" -istart "+ItoA(iter+1);
  DF.insert_comment(comment);
  DF.insert_comment(cline);
  DF.insert_comment("columns: model fraction (1); mirror fraction (2); 1000x signal change (3)");
  fn_tmp=fn_root;
  if (iter>=0) {
    fn_tmp+="_it";
    fn_tmp.compose(fn_tmp,iter,"");
  }
  fn_tmp+=".log";
  DF.write(fn_tmp);

}

// Write out results  ========================================================
void Prog_MLalign2D_prm::MLalign2D() _THROW {

  int c,nn,imgno,opt_refno,iter ;
  double LL,sumw_allrefs,convv,avecorr;
  bool converged;
  vector<double> conv;
  double wsum_sigma_noise, wsum_sigma_offset;
  vector<matrix2D<double > > wsum_Mref;
  vector<ImageXmipp> Ireg;
  vector<double> sumw,sumw_mirror;
  matrix2D<double> P_phi,Mr2,Maux;
  FileName fn_img,fn_tmp;
  matrix1D<double> oneline(0);
  DocFile DFo,DFf;
  SelFile SFo,SFa;
  
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  DFo.reserve(2*SF.ImgNo()+1);
  DFf.reserve(2*SFr.ImgNo()+4);
  SFa.reserve(Niter*n_ref);
  SFa.clear();

  if (fast_mode && do_precenter) precenter_images();

  // Loop over all iterations
  for (iter=istart; iter<=Niter; iter++) {

    if (verb>0) cerr << "  multi-reference refinement:  iteration " << iter <<" of "<< Niter<<endl;

    FOR_ALL_MODELS() { Iold[refno]()=Iref[refno](); };

    conv.clear();
    DFo.clear();
    DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), Corr (8)");

    // Pre-calculate pdfs
    calculate_pdf_phi(sigma_offset,P_phi,Mr2);

    // Integrate over all images
    ML_sum_over_all_images(SF,Iref, 
			   P_phi,Mr2,alpha_k,mirror_fraction, 
			   LL,avecorr,DFo, 
			   wsum_Mref,
			   wsum_sigma_noise,wsum_sigma_offset,sumw,sumw_mirror); 

    // Update model parameters
    sumw_allrefs=0.;
    FOR_ALL_MODELS() {
      if (sumw[refno]>0.) {
	Iref[refno]()=wsum_Mref[refno];
	Iref[refno]()/=sumw[refno];
	Iref[refno].weight()=sumw[refno];
	sumw_allrefs+=sumw[refno];
	if (do_esthetics) MAT_ELEM(Iref[refno](),0,0)=
			    (MAT_ELEM(Iref[refno](),1,0)+MAT_ELEM(Iref[refno](),0,1)+
			     MAT_ELEM(Iref[refno](),-1,0)+MAT_ELEM(Iref[refno](),0,-1))/4;
	if (!fix_fractions) alpha_k[refno]=sumw[refno]/SF.ImgNo();
	if (!fix_fractions) mirror_fraction[refno]=sumw_mirror[refno]/sumw[refno];
      } else {
	Iref[refno].weight()=0.;
	Iref[refno]().init_zeros();
	alpha_k[refno]=0.;
	mirror_fraction[refno]=0.;
      }
    }
    if (!fix_sigma_offset) sigma_offset=sqrt(wsum_sigma_offset/(2*SF.ImgNo()));
    if (!fix_sigma_noise)  sigma_noise=sqrt(wsum_sigma_noise/(SF.ImgNo()*dim*dim));

   // Check convergence 
    converged=true;
    FOR_ALL_MODELS() {
      if (Iref[refno].weight()>0.) {
	Maux=mul_elements(Iold[refno](),Iold[refno]());
	convv=1/(Maux.compute_avg());
	Maux=Iold[refno]()-Iref[refno]();
	Maux=mul_elements(Maux,Maux);
	convv*=Maux.compute_avg();
	conv.push_back(convv);
	if (convv>eps) converged=false;
      } else {
	conv.push_back(-1.);
      }
    }

    if (write_intermediate) {
      write_output_files(iter,SFa,DFf,sumw_allrefs,LL,avecorr,conv);
    } else {
      // Output current parameter values to screen 
      if (verb>0) { 
	cout <<"  iter "<<iter<<" noise= "<<FtoA(sigma_noise,10,7)<<" offset= "<<FtoA(sigma_offset,10,7);
	if (!LSQ_rather_than_ML) cout <<" LL= "<<LL<<" R= "<<avecorr;
	cout <<endl<<"  Model  fraction  mirror-fraction "<<endl;
	FOR_ALL_MODELS() {cout <<"  "<<ItoA(refno+1,5)<<" "<<FtoA(alpha_k[refno],10,7)<<" "<<FtoA(mirror_fraction[refno],10,7)<<endl;}
      }
    }
    
    if (converged) {
      if (verb>0) cerr <<" Optimization converged!"<<endl;
      break;
    }

  } // end loop iterations

  // Write out converged structures
  write_output_files(-1,SFa,DFf,sumw_allrefs,LL,avecorr,conv);

  if (write_docfile) {
    // Write out docfile with optimal transformation & references
    fn_img=fn_root+".doc";
    DFo.write(fn_img);
  }
  if (write_selfiles) {
    // Also write out selfiles of all experimental images, classified according to optimal reference image
    FOR_ALL_MODELS() {
      DFo.go_beginning();
      SFo.clear();
      for (int n=0; n<DFo.dataLineNo(); n++ ) {
	DFo.next();
	fn_img=((DFo.get_current_line()).get_text()).erase(0,3);
	DFo.adjust_to_data_line();
	if ((refno+1)==(int)DFo(5)) SFo.insert(fn_img,SelLine::ACTIVE);
      }
      fn_tmp=fn_root+"_ref";
      fn_tmp.compose(fn_tmp,refno+1,"sel");
      SFo.write(fn_tmp);
    }
  }

}

