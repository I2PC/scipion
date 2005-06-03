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
  } else {
    for (int i=1; i<argc; i++) {
      cline=cline+(string)argv[i]+" ";
    }
  }
    
  // Read command line
  if (check_param(argc,argv,"-show_all_options")) { usage(); extended_usage();}
  n_ref=AtoI(get_param(argc,argv,"-nref","0"));
  fn_ref=get_param(argc,argv,"-ref","");
  SF.read(get_param(argc,argv,"-i"));
  SF.ImgSize(dim,dim);
  fn_root=get_param(argc,argv,"-o","MLalign2D");
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
  fix_fractions=check_param(argc,argv,"-fix_fractions");
  fix_sigma_offset=check_param(argc,argv,"-fix_sigma_offset");
  fix_sigma_noise=check_param(argc,argv,"-fix_sigma_noise");
  verb=AtoI(get_param(argc,argv,"-verb","1"));
  fast_mode=check_param(argc,argv,"-fast");
  C_fast=AtoF(get_param(argc,argv,"-C","1e-12"));
  LSQ_rather_than_ML=check_param(argc,argv,"-LSQ");
  max_shift=AtoF(get_param(argc,argv,"-max_shift","5"));
  istart=AtoI(get_param(argc,argv,"-istart","1"));
  // Hidden arguments
  Paccept_fast=AtoF(get_param(argc,argv,"-Paccept","0.0"));

}

// Show ====================================================================
void Prog_MLalign2D_prm::show() {

  if (verb>0) {
    // To screen
    cerr << "--> Maximum-likelihood multi-reference refinement "<<endl; 
    cerr << "  Input images            : "<< SF.name()<<" ("<<SF.ImgNo()<<")"<<endl;
    cerr << "  Reference images        : "<< fn_ref<<" ("<<SFr.ImgNo()<<")"<<endl;
    cerr << "  Output rootname         : "<< fn_root<<endl;
    cerr << "  Number of iterations    : "<< Niter<<endl;
    cerr << "  Stopping criterium      : "<< eps<<endl;
    cerr << "  initial sigma noise     : "<< sigma_noise<<endl;
    cerr << "  initial sigma offset    : "<< sigma_offset<<endl;
    cerr << "  Psi sampling interval   : "<< psi_step<<endl;
    if (fn_frac!="") {
      cerr << "  -> Read initial model fractions from "<< fn_frac<<endl;
    }
    if (do_mirror) {
      cerr << "  -> Check mirror image of each reference as well."<<endl;
    }
    if (fast_mode) {
      cerr << "  -> Use fast, reduced search-space approach with C = "<<C_fast<<endl;
    }
    if (LSQ_rather_than_ML) {
      cerr << "  -> Use least-squares instead of maximum likelihood target."<<endl;
      if (max_shift>0.) 
	cerr << "  -> Use maximum shift criterion in translation search: max_shift = "<<max_shift<<endl;
    }
    if (write_docfile) {
      cerr << "  -> Write out docfile with most likely angles & translations. "<<endl;
    }
    if (write_selfiles) {
      cerr << "  -> Write out selfiles with most likely reference assignments. "<<endl;
    }

    // Hidden stuff
    if (!write_intermediate) {
      cerr << "  -> Do not write out images after each iteration."<<endl;
    }
    if (fix_fractions && !LSQ_rather_than_ML) {
      cerr << "  -> Do not update estimates of model fractions."<<endl;
    }
    if (fix_sigma_offset && !LSQ_rather_than_ML) {
      cerr << "  -> Do not update sigma-estimate of origin offsets."<<endl;
    }
    if (fix_sigma_noise && !LSQ_rather_than_ML) {
      cerr << "  -> Do not update sigma-estimate of noise."<<endl;
    }
    cerr << " -----------------------------------------------------------------"<<endl;
  }

} 

// Usage ===================================================================
void Prog_MLalign2D_prm::usage(bool ML3D) {
  if (!ML3D) {
    cerr << "Usage:  MLalign2D [options] "<<endl;
    cerr << "   -i <selfile>                : Selfile with input images \n";
    cerr      << "   -ref <selfile/image>        : Selfile with initial reference images/single reference image \n";
    cerr      << "      OR -nref <int>               OR number of bias-free references to generate automatically\n";
    cerr      << " [ -o <rootname=\"MLalign2D\"> ] : Output rootname \n";
  }
  cerr << " [ -noise <float=1> ]          : Expected standard deviation for pixel noise \n";
  cerr << " [ -offset <float=3> ]         : Expected standard deviation for origin offset [pix]\n";
  if (!ML3D) cerr      << " [ -mirror ]                   : Also check mirror image of each reference \n";
  cerr << " [ -output_docfile ]           : Write out docfile with most likely angles & translations \n";
  cerr << " [ -output_selfiles ]          : Write out selfiles with most likely reference assignments \n";
  cerr << " [ -fast ]                     : Use pre-centered images to pre-calculate significant orientations\n";
  if (!ML3D) cerr      << " [ -show_all_options ]         : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_MLalign2D_prm::extended_usage(bool ML3D) {
  cerr << "Additional options: "<<endl;
  cerr << " [ -eps <float=5e-5> ]         : Stopping criterium \n";
  cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
  cerr << " [ -psi_step <float=5> ]       : In-plane rotation sampling interval [deg]\n";
  cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
  cerr << " [ -C <double=1e-12> ]         : Significance criterion for fast approach \n";
  if (!ML3D) cerr << " [ -restart <logfile> ]        : restart a run with all parameters as in the logfile \n";
  if (!ML3D) cerr << " [ -istart <int> ]             : number of initial iteration \n";
  cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
  cerr << " [ -fix_sigma_offset]          : Do not re-estimate the standard deviation in the origin offsets \n";
  cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
  cerr << " [ -LSQ ]                      : Use maximum cross-correlation instead of maximum likelihood target \n";
  cerr << " [ -max_shift <float=5>]       : For LSQ only: maximum allowed shift [pix] \n";
  cerr << endl;
  exit(1);
}

// Read all input selfiles in memory 
void Prog_MLalign2D_prm::produce_Side_info() _THROW {

  FileName fn_img;
  ImageXmipp img;
  headerXmipp head;
  matrix1D<double> offsets(2);
  matrix2D<double> A(3,3);
  vector<matrix1D<double> > dum;
  int c,refno=0;
  float xx,yy;

  // Set nr_psi
  nr_psi=ROUND(90./psi_step);
  psi_step=90./nr_psi;
 
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

  // Read in all reference images in memory
  n_ref=0;
  SFr.go_beginning();
  while ((!SFr.eof())) {
    img.read(SFr.NextImg(),FALSE,FALSE,TRUE,FALSE);
    img().set_Xmipp_origin();
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

  if (LSQ_rather_than_ML) {
    if (max_shift<0) max_shift=dim/2;
    fix_sigma_noise=true;
    fix_sigma_offset=true;
    fix_fractions=true;
  }

  // For fast_mode: fill imgs_offsets vectors
  SF.go_beginning();
  c=0;
  offsets(0)=-999;
  offsets(1)=-999;
  while ((!SF.eof())) {
    SF.NextImg();
    imgs_offsets.push_back(dum);
    if (fast_mode) {
      for (int refno=0; refno<n_ref; refno++) {
	imgs_offsets[c].push_back(offsets);
	imgs_offsets[c].push_back(offsets);
      }
    }
    c++;
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

// Calculate probability density function of all in-plane transformations phi
void Prog_MLalign2D_prm::calculate_pdf_phi() _THROW {

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

}

// Rotate reference for all models and rotations and fill Fref vectors =============
void Prog_MLalign2D_prm::rotate_reference(vector<ImageXmipp> &Iref,bool &also_real_space, vector <vector< matrix2D<double> > > &Mref,
			vector <vector< matrix2D<complex<double> > > > &Fref) _THROW {

  double AA,stdAA,sumAA=0.,psi,dum,avg,mean_ref,stddev_ref,dummy;
  matrix2D<double> Maux;
  matrix2D<complex<double> > Faux;
  vector<matrix2D<complex <double> > > dumF;
  vector<matrix2D<double> > dumM;
  matrix2D<int> mask, omask;

  Maux.init_zeros(dim,dim);
  Maux.set_Xmipp_origin();
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  omask.resize(dim,dim);
  omask.set_Xmipp_origin();
  BinaryCircularMask(mask,dim/2,INNER_MASK);
  BinaryCircularMask(omask,dim/2,OUTSIDE_MASK);

  Fref.clear();
  Mref.clear();
  A2.clear();

  FOR_ALL_MODELS() {
    Mref.push_back(dumM);
    Fref.push_back(dumF);
    compute_stats_within_binary_mask(omask,Iref[refno](),dum,dum,avg,dum);
    FOR_ALL_ROTATIONS() {
      // Add arbitrary number (SMALLANGLE) to avoid 0-degree rotation (lacking interpolation effects)
      psi=(double)(ipsi*90./nr_psi)+SMALLANGLE; 
      Maux=Iref[refno]().rotate(psi,DONT_WRAP);
      apply_binary_mask(mask,Maux,Maux,avg);
      // Normalize the magnitude of the rotated references to 1st rot of that ref
      // This is necessary because interpolation due to rotation can lead to lower overall Fref 
      // This would result in lower probabilities for those rotations
      AA=Maux.sum2();      
      if (ipsi==0) {
	stdAA=AA;
	sumAA+=AA;
	A2.push_back(AA);
      }
      // Subtract mean_ref from image prior to FFT for maxCC
      if (LSQ_rather_than_ML) {
	Maux.compute_stats(mean_ref,stddev_ref,dummy,dummy);
	Maux-=mean_ref;
	if (ipsi==0) A2[refno]=stddev_ref;
      }
      if (also_real_space) {
	Mref[refno].push_back(Maux);
	if (AA>0) Mref[refno][ipsi]*=sqrt(stdAA/AA);
      }
      FourierTransformHalf(Maux,Faux);
      Faux*=dim*dim;
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux) {
	dMij(Faux,i,j)=conj(dMij(Faux,i,j));
      }
      Fref[refno].push_back(Faux);
      if (AA>0) Fref[refno][ipsi]*=sqrt(stdAA/AA);

    }
  }
  A2mean=sumAA/n_ref;

  // Check that the differences in powers are not to large 
  // for the numerical precision of the exp function...
  FOR_ALL_MODELS() {
    if (!LSQ_rather_than_ML && 
	ABS((A2mean-A2[refno])/(2*sigma_noise*sigma_noise))>1000.) {
      cerr <<"mean power= "<<A2mean<<" power reference "<<refno<<"= "<<A2[refno]<<endl;
      REPORT_ERROR(1,"Very large differences in powers of references give numerical problems, are the images normalized?");
    }
  }

}

// Collect all rotations and sum to update Iref() for all models ==========
void Prog_MLalign2D_prm::reverse_rotate_reference(vector <vector< matrix2D<complex<double> > > > &Fref, 
						  vector <vector< matrix2D<double> > > &Mref, bool &real_space, 
						  vector<matrix2D<double > > &Mnew ) _THROW {

  double psi,dum,avg,ang;
  int s;
  matrix2D<double> Maux,Maux2;
  matrix2D<int> mask, omask;
  Maux.resize(dim,dim);
  Maux2.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Maux2.set_Xmipp_origin();

  Mnew.clear();
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  omask.resize(dim,dim);
  omask.set_Xmipp_origin();
  BinaryCircularMask(mask,dim/2,INNER_MASK);
  BinaryCircularMask(omask,dim/2,OUTSIDE_MASK);

  FOR_ALL_MODELS() {
    Maux.init_zeros();
    Mnew.push_back(Maux);
    FOR_ALL_ROTATIONS() {
      // Add arbitrary number to avoid 0-degree rotation without interpolation effects
      psi=(double)(ipsi*90./nr_psi)+SMALLANGLE;
      if (real_space) {
	Maux=Mref[refno][ipsi];
      } else {
	InverseFourierTransformHalf(Fref[refno][ipsi],Maux,dim);
	Maux/=dim*dim;
	CenterFFT(Maux,true);
      }
      compute_stats_within_binary_mask(omask,Maux,dum,dum,avg,dum);
      Maux2=Maux.rotate(-psi,DONT_WRAP);
      apply_binary_mask(mask,Maux2,Maux2,avg);
      Mnew[refno]+=Maux2;
    }  
  }

}


// Pre-selection of significant refno and ipsi, based on current optimal translation ====================================
void Prog_MLalign2D_prm::preselect_significant_model_phi(matrix2D<double> &Mimg, vector<matrix1D<double> > &offsets, 
							 vector <vector< matrix2D<double > > > &Mref, 
							 matrix2D<int> &Msignificant) _THROW {

  matrix2D<double> Maux,Maux2,Mdsig(n_ref,nr_psi*nr_flip);
  double ropt,sigma_noise2,aux,fracpdf;
  double CC,maxCC=-99.e99,maxweight=-99.e99;
  int irot,irefmir;
  vector<double> maxw_ref(2*n_ref);
  matrix1D<double> trans(2);

  randomize_random_generator();

  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Maux2.resize(dim,dim);
  Maux2.set_Xmipp_origin();
  sigma_noise2=sigma_noise*sigma_noise;

  // Flip images and calculate correlations and maximum correlation
  FOR_ALL_MODELS() {
    FOR_ALL_FLIPS() {
      irefmir=FLOOR(iflip/4)*n_ref+refno;
      // Do not trust optimal offsets if they are larger than 3*sigma_offset: 
      ropt=sqrt(offsets[irefmir](0)*offsets[irefmir](0)+offsets[irefmir](1)*offsets[irefmir](1));
      if (ropt>3*sigma_offset)  {
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  dMij(Mdsig,refno,irot)=1.;
	}
      } else {
	Maux=Mimg.translate(offsets[irefmir],true);
	apply_geom(Maux2,F[iflip],Maux,IS_INV,WRAP);
	FOR_ALL_ROTATIONS() {
	  irot=iflip*nr_psi+ipsi;
	  CC=0.;
	  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Maux2) {
	    CC+=dMij(Maux2,i,j)*dMij(Mref[refno][ipsi],i,j);
	  }
	  dMij(Mdsig,refno,irot)=CC;
	  if (CC>maxCC) maxCC=CC;
	}
      }
    }
  }

  // Now that we have maxCC calculate the weighting matrices and maxweight
  FOR_ALL_MODELS() {
    FOR_ALL_FLIPS() {
      irefmir=FLOOR(iflip/4)*n_ref+refno;
      FOR_ALL_ROTATIONS() {
	irot=iflip*nr_psi+ipsi;
	if (iflip<4) fracpdf=alpha_k[refno]*(1.-mirror_fraction[refno]);
	else fracpdf=alpha_k[refno]*mirror_fraction[refno];
	aux=dMij(Mdsig,refno,irot)-maxCC;
	aux=exp(aux/sigma_noise2)*fracpdf;
	dMij(Mdsig,refno,irot)=aux;
	if (aux>maxw_ref[irefmir]) maxw_ref[irefmir]=aux;
      }
    }
  }

  // Now that we have maxweight calculate which weighting matrices are significant
  // Set Vweight to zero if insignificant
  FOR_ALL_MODELS() {
    FOR_ALL_FLIPS() {
      irefmir=FLOOR(iflip/4)*n_ref+refno;
      FOR_ALL_ROTATIONS() {
	irot=iflip*nr_psi+ipsi;
	if (dMij(Mdsig,refno,irot)>=C_fast*maxw_ref[irefmir]) dMij(Msignificant,refno,irot)=1;
	else {
	  // Allow a fraction Paccept_fast to have a complete search
	  if (rnd_unif(0,1)<Paccept_fast) dMij(Msignificant,refno,irot)=1;
	  else dMij(Msignificant,refno,irot)=0;
	}
      }
    }
  }

}


// Maximum Likelihood calculation for one image ============================================
// Integration over all translation, given  model and in-plane rotation
void Prog_MLalign2D_prm::ML_integrate_model_phi_trans(
     matrix2D<double> &Mimg, vector <vector< matrix2D<complex<double> > > > &Fref, 
     matrix2D<int> &Msignificant,
     vector <vector< matrix2D<complex<double> > > > &Fwsum_imgs, 
     double &wsum_sigma_noise, double &wsum_sigma_offset, 
     vector<double> &sumw, vector<double> &sumw_mirror, 
     double &LL, double &fracweight, int &opt_refno, double &opt_psi, 
     matrix1D<double> &opt_offsets, vector<matrix1D<double> > &opt_offsets_ref) _THROW {

  matrix2D<double> Maux,Mdzero;
  matrix2D<complex<double> > Fimg, Faux;
  vector<matrix2D<complex<double> > > Fimg_flip;
  vector<matrix2D<double> > dumM;
  vector<bool> dumb;
  vector <vector< matrix2D<double> > > Mweight;
  vector<double> refw(n_ref), refw_mirror(n_ref);
  double sigma_noise2,XiA,Xi2,aux,fracpdf,sigw,add,max,corrA2,maxc=-99.e99;
  double sum,wsum_corr=0., sum_refw=0., wsum_A2=0., maxweight=-99.e99;
  int irot,irefmir,sigdim,xmax,ymax;
  int ioptx=0,iopty=0,ioptpsi=0,ioptflip=0,imax=0;
  if (fast_mode) imax=n_ref*nr_flip/4;
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
  

  // Only translations smaller than 6 sigma_offset are considered!
  // This saves a lot of memory and CPU! (typically a factor 2, depending on sigma_offset vs. dim)
  sigdim=2*CEIL(sigma_offset*6);
  sigdim++; // (to get uneven number)
  sigdim=MIN(dim,sigdim);

  if (fast_mode) for (int i=0; i<imax; i++) maxw_ref[i]=-99.e99;
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Faux.resize(dim,dim);
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
	    MAT_ELEM(Mweight[refno][irot],i,j)=MAT_ELEM(Maux,i,j);
	  }
	  max=Mweight[refno][irot].compute_max();
	  if (max>maxc) maxc=max;
	} else {
	  Mweight[refno].push_back(Mdzero);
	}
      }
    }
  }

  // Now that we have maxc calculate the weighting matrices and maxweight
  FOR_ALL_MODELS() {
    refw[refno]=0.;
    refw_mirror[refno]=0.;
    // There is a numerical problem with the following line 
    // when the difference in power of the references is large, the exp becomes 0/inf.
    // For now, I dont know what to do... Put a warning in rotate_reference
    // If the images are properly normalized this does not seem to be so much of a problem?
    corrA2=exp((A2mean-A2[refno])/(2*sigma_noise2));
    FOR_ALL_ROTATIONS() {
      FOR_ALL_FLIPS() {
	irot=iflip*nr_psi+ipsi;
	irefmir=FLOOR(iflip/4)*n_ref+refno;
	if (dMij(Msignificant,refno,irot)) {
	  if (iflip<4) fracpdf=alpha_k[refno]*(1.-mirror_fraction[refno]);
	  else fracpdf=alpha_k[refno]*mirror_fraction[refno];
	  sum=0.;
	  FOR_ALL_ELEMENTS_IN_MATRIX2D(Mweight[refno][irot]) {
	    aux=MAT_ELEM(Mweight[refno][irot],i,j)-maxc;
	    aux=exp(aux/sigma_noise2)*fracpdf*MAT_ELEM(P_phi,i,j)*corrA2;
	    wsum_corr+=aux*MAT_ELEM(Mweight[refno][irot],i,j);
	    MAT_ELEM(Mweight[refno][irot],i,j)=aux;
	    sum+=aux;
	  }
	  wsum_A2+=sum*A2[refno];
	  if (iflip<4) refw[refno]+=sum;
	  else refw_mirror[refno]+=sum;
	  Mweight[refno][irot].max_index(ymax,xmax);
	  max=MAT_ELEM(Mweight[refno][irot],ymax,xmax);
	  if (max>maxweight) {
	    maxweight=max;
	    iopty=ymax;
	    ioptx=xmax;
	    ioptpsi=ipsi;
	    ioptflip=iflip;
	    opt_refno=refno;
	  }
	  if (fast_mode && max>maxw_ref[irefmir]) {
	    maxw_ref[irefmir]=max;
	    iopty_ref[irefmir]=ymax;
	    ioptx_ref[irefmir]=xmax;
	    ioptflip_ref[irefmir]=iflip;
	  } 
	}
      }
    }
    sum_refw+=refw[refno]+refw_mirror[refno];
  }

  // Normalize all weighted sums by sum_refw such that sum over all weights is one!
  // And accumulate the FT of the weighted, shifted images.
  wsum_sigma_noise+=(wsum_A2/sum_refw)+Xi2-(2*wsum_corr/sum_refw);
  FOR_ALL_MODELS() {
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

  // Calculate optimal transformation parameters
  if (fast_mode) {
    for (int i=0; i<imax;i++) {
      opt_offsets_ref[i](0)=-(double)ioptx_ref[i]*DIRECT_MAT_ELEM(F[ioptflip_ref[i]],0,0)-
	                     (double)iopty_ref[i]*DIRECT_MAT_ELEM(F[ioptflip_ref[i]],0,1);
      opt_offsets_ref[i](1)=-(double)ioptx_ref[i]*DIRECT_MAT_ELEM(F[ioptflip_ref[i]],1,0)-
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
  // 2nd term: for subtracting maxc
  // 3rd term: for only considering Xi*A instead of (A-Xi)^2
  // 4th term: for (sqrt(2pi)*sigma_noise)^-1 term in formula (12) Sigworth (1998)
  LL+= log(sum_refw) + maxc/sigma_noise2 - (A2mean+Xi2)/(2*sigma_noise2) 
       - dim*dim*log(2.50663*sigma_noise);

}


void Prog_MLalign2D_prm::LSQ_search_model_phi_trans(matrix2D<double> &Mimg, vector <vector< matrix2D<complex<double> > > > &Fref, 
						    double &max_shift, 
						    vector <vector< matrix2D<double> > > &Msum_imgs, 
						    vector<double> &sumw, vector<double> &sumw_mirror, 
						    double &maxCC, int &opt_refno, double &opt_psi, 
						    matrix1D<double> &opt_offsets) _THROW {

  matrix2D<double> Maux,Maux2;
  matrix2D<complex<double> > Fimg, Faux;
  double sigma_noise2,aux,avg,std,CC;
  int irot,sigdim,xmax,ymax;
  int ioptx=0,iopty=0,ioptpsi=0,ioptflip=0,imax=0;
  double stddev_img,mean_img,dummy;

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

  maxCC=-99.e99; 
  Maux.resize(dim,dim);
  Maux.set_Xmipp_origin();
  Maux2.resize(dim,dim);
  Maux2.set_Xmipp_origin();
  Faux.resize((int)(dim/2)+1,dim);
  sigma_noise2=sigma_noise*sigma_noise;
  matrix2D<int> shiftmask;

  if (max_shift>0.) {
    shiftmask.resize(dim,dim);
    shiftmask.set_Xmipp_origin();
    BinaryCircularMask(shiftmask,max_shift,INNER_MASK);
  }

  Mimg.compute_stats(mean_img,stddev_img,dummy,dummy);
  Maux2=Mimg;
  Maux2-=mean_img;

  // Flip images and calculate correlation matrices and maximum correlation
  FOR_ALL_FLIPS() {
    apply_geom(Maux,F[iflip],Maux2,IS_INV,WRAP);
    FourierTransformHalf(Maux,Fimg);
    Fimg*=dim*dim;
    FOR_ALL_MODELS() {
      FOR_ALL_ROTATIONS() {
	irot=iflip*nr_psi+ipsi;
	mul_elements(Fimg,Fref[refno][ipsi],Faux);
	InverseFourierTransformHalf(Faux,Maux,dim);
	Maux/=dim*dim;
	CenterFFT(Maux,true);
	if (max_shift>0.) apply_binary_mask(shiftmask,Maux,Maux,0.);
	Maux.max_index(ymax,xmax);
	CC=MAT_ELEM(Maux,ymax,xmax);
	CC/=A2[refno]*stddev_img; // A2[refno] now holds stddev_ref!
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
  maxCC/=dim*dim;

  // Calculate optimal transformation parameters
  opt_offsets(0)=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],0,0)-(double)iopty*DIRECT_MAT_ELEM(F[ioptflip],0,1);
  opt_offsets(1)=-(double)ioptx*DIRECT_MAT_ELEM(F[ioptflip],1,0)-(double)iopty*DIRECT_MAT_ELEM(F[ioptflip],1,1);
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
			  double &wsum_sigma_noise, double &wsum_sigma_offset, vector<double> &sumw, vector<double> &sumw_mirror) _THROW {


  ImageXmipp img;
  FileName fn_img;
  vector <vector< matrix2D<double> > > Mref,Msum_imgs;
  vector <vector< matrix2D<complex<double> > > > Fref,Fwsum_imgs;
  vector<matrix2D<complex <double> > > dum;
  vector<matrix2D<double> > dum2;
  matrix2D<complex <double> > Fdzero;
  matrix2D<double>  Mdzero;
  matrix2D<int> Msignificant(n_ref,nr_psi*nr_flip);
  matrix1D<double> dataline(8), opt_offsets(2);  
  double opt_psi,opt_flip,maxcorr;
  int c,nn,imgno,opt_refno;

  // Generate (FT of) each rotated version of all references 
  rotate_reference(Iref,fast_mode,Mref,Fref);

  // Initialize
  nn=SF.ImgNo();
  if (verb>0) init_progress_bar(nn);
  c=MAX(1,nn/60);
  Fwsum_imgs.clear();
  Msum_imgs.clear();
  sumw.clear();
  sumw_mirror.clear();
  Fdzero.resize((int)(dim/2)+1,dim);
  Mdzero.resize(dim,dim);
  Mdzero.set_Xmipp_origin();
  LL=0.;
  wsum_sigma_noise=0.;
  wsum_sigma_offset=0.;
  sumcorr=0.;
  FOR_ALL_MODELS() {
    sumw.push_back(0.);
    sumw_mirror.push_back(0.);
    if (LSQ_rather_than_ML) {
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
    img.read(fn_img,FALSE,FALSE,FALSE,FALSE);
    img().set_Xmipp_origin();

    // Perform the integration over all shifts for all significant models and phis
    if (LSQ_rather_than_ML) {

      LSQ_search_model_phi_trans(img(),Fref,max_shift,Msum_imgs,sumw,sumw_mirror,
      				 maxcorr,opt_refno,opt_psi,opt_offsets);

    } else {

      if (fast_mode) 
	// Pre-calculate which models and phis are significant at zero origin shifts
	preselect_significant_model_phi(img(),imgs_offsets[imgno],Mref,Msignificant);
      else  Msignificant.init_constant(1);

      ML_integrate_model_phi_trans(img(),Fref,Msignificant,
				   Fwsum_imgs,wsum_sigma_noise, wsum_sigma_offset,sumw,sumw_mirror, 
				   LL,maxcorr,opt_refno,opt_psi,opt_offsets,imgs_offsets[imgno]);
    }

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
      dataline(3)=opt_offsets(0);          // Xoff
      dataline(4)=opt_offsets(1);          // Yoff
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

  reverse_rotate_reference(Fwsum_imgs,Msum_imgs,LSQ_rather_than_ML,wsum_Mref);

}

// Update all model parameters
void Prog_MLalign2D_prm::update_parameters(vector<matrix2D<double> > &wsum_Mref,
					   double &wsum_sigma_noise, double &wsum_sigma_offset, 
					   vector<double> &sumw, vector<double> &sumw_mirror, 
					   double &sumcorr, double &sumw_allrefs) {

  // Pre-calculate sumw_allrefs
  sumw_allrefs=0.;
  FOR_ALL_MODELS() {
    sumw_allrefs+=sumw[refno];
  }

  FOR_ALL_MODELS() {
    if (sumw[refno]>0.) {
      Iref[refno]()=wsum_Mref[refno];
      Iref[refno]()/=sumw[refno];
      Iref[refno].weight()=sumw[refno];
      if (!fix_fractions) alpha_k[refno]=sumw[refno]/sumw_allrefs;
      if (!fix_fractions) mirror_fraction[refno]=sumw_mirror[refno]/sumw[refno];
    } else {
      Iref[refno].weight()=0.;
      Iref[refno]().init_zeros();
      alpha_k[refno]=0.;
      mirror_fraction[refno]=0.;
    }
  }
  if (!fix_sigma_offset) sigma_offset=sqrt(wsum_sigma_offset/(2*sumw_allrefs));
  if (!fix_sigma_noise)  sigma_noise=sqrt(wsum_sigma_noise/(sumw_allrefs*dim*dim));
  
  sumcorr/=sumw_allrefs;

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
    if (LSQ_rather_than_ML) cout <<"  iter "<<iter<<" <CC>= "+FtoA(sumcorr,10,5);
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

void Prog_MLalign2D_prm::write_output_files(const int iter, SelFile &SF, DocFile &DF, 
			 double &sumw_allrefs, double &LL, double &avecorr, vector<double> &conv ) _THROW {

  FileName fn_tmp;
  matrix1D<double> fracline(3);
  string comment;
 
  DF.clear();
  SF.clear();

  // Write out current reference images and fill sel & log-file
  FOR_ALL_MODELS() {
    fn_tmp=fn_root;
    if (iter>=0) {
      fn_tmp+="_it";
      fn_tmp.compose(fn_tmp,iter,"");
    }
    fn_tmp+="_ref";
    fn_tmp.compose(fn_tmp,refno+1,"");
    fn_tmp=fn_tmp+".xmp";
    Iref[refno].write(fn_tmp);
    SF.insert(fn_tmp,SelLine::ACTIVE);
    fracline(0)=alpha_k[refno];
    fracline(1)=mirror_fraction[refno];
    fracline(2)=1000*conv[refno]; // Output 1000x the change for precision
    DF.insert_comment(fn_tmp);
    DF.insert_data_line(fracline);
  }
  // Write out sel & log-file
  fn_tmp=fn_root;
  if (iter>=0) {
      fn_tmp+="_it";
      fn_tmp.compose(fn_tmp,iter,"sel");
  }
  SF.write(fn_tmp);

  DF.go_beginning();
  comment="MLalign2D-logfile: Number of images= "+FtoA(sumw_allrefs);
  if (LSQ_rather_than_ML) comment+=" <CC>= "+FtoA(avecorr,10,5);
  else { 
    comment+=" LL= "+FtoA(LL,10,5)+" <Pmax/sumP>= "+FtoA(avecorr,10,5);
    DF.insert_comment(comment);
    comment="-noise "+FtoA(sigma_noise,10,7)+" -offset "+FtoA(sigma_offset,10,7)+" -istart "+ItoA(iter+1);
  }
  DF.insert_comment(comment);
  DF.insert_comment(cline);
  DF.insert_comment("columns: model fraction (1); mirror fraction (2); 1000x signal change (3)");
  fn_tmp=fn_root;
  if (iter>=0) {
      fn_tmp+="_it";
      fn_tmp.compose(fn_tmp,iter,"log");
  } else fn_tmp+=".log";
  DF.write(fn_tmp);

}

