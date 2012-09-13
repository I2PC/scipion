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
#include "ml_align3d.h"

// Read arguments ==========================================================
void Prog_MLalign3D_prm::read(int argc, char **argv) {

  // Read command line
  if (checkParameter(argc,argv,"-more_options")) { usage(); extended_usage();}
  fn_ref=getParameter(argc,argv,"-ref");
  SF.read(getParameter(argc,argv,"-i"));
  nr_img=SF.ImgNo();
  fn_root=getParameter(argc,argv,"-o","MLalign3D");
  Niter=textToInteger(getParameter(argc,argv,"-iter","100"));
  fn_frac=getParameter(argc,argv,"-frac","");
  fix_fractions=checkParameter(argc,argv,"-fix_fractions");
  fix_sigma_noise=checkParameter(argc,argv,"-fix_sigma_noise");
  fix_sigma_offset=checkParameter(argc,argv,"-fix_sigma_offset");
  verb=textToInteger(getParameter(argc,argv,"-verb","1"));
  istart=textToInteger(getParameter(argc,argv,"-istart","1"));
  getThreeDoubleParams(argc,argv,"-rot",rot0,rotF,rot_step,0,0,1);
  getThreeDoubleParams(argc,argv,"-tilt",tilt0,tiltF,tilt_step,0,0,1);
  getThreeDoubleParams(argc,argv,"-psi",psi0,psiF,psi_step,0,0,1);
  sigma_offset=textToFloat(getParameter(argc,argv,"-offset","0"));
  max_shift=textToFloat(getParameter(argc,argv,"-max_shift","-1"));
  fn_doc=getParameter(argc,argv,"-doc","");
  fn_wlist=getParameter(argc,argv,"-wedge","");
  fn_sym=getParameter(argc,argv,"-sym","");
  ccf_mode=checkParameter(argc,argv,"-CCF");
  sigma_noise2=textToFloat(getParameter(argc,argv,"-noise","1"));
  sigma_noise2*=sigma_noise2;
  fn_solv=getParameter(argc,argv,"-solvent","");
  fn_solv2=getParameter(argc,argv,"-solvent2","");
  theta=textToFloat(getParameter(argc,argv,"-theta","0"));
  theta0=textToFloat(getParameter(argc,argv,"-theta0","0"));
  theta_step=textToFloat(getParameter(argc,argv,"-theta_step","0"));
  fn_misalign=getParameter(argc,argv,"-misalign","");
  // Hidden, for CCF-mode only...
  do_subtract_current=checkParameter(argc,argv,"-subtract_current");
  fn_mask=getParameter(argc,argv,"-mask","");

  if (ccf_mode) fix_sigma_offset=true;
  if (max_shift<0) max_shift=(double)(dim/2);

}

// Show ====================================================================
void Prog_MLalign3D_prm::show() {

  if (verb>0) {
    // To screen
    std::cerr << "--> Maximum-likelihood multi-reference refinement "<<std::endl; 
    std::cerr << "  Tomographic volumes     : "<< SF.name()<<" ("<<SF.ImgNo()<<")"<<std::endl;
    std::cerr << "  Reference volumes       : "<< fn_ref<<" ("<<SFr.ImgNo()<<")"<<std::endl;
    std::cerr << "  Output rootname         : "<< fn_root<<std::endl;
    std::cerr << "  Number of iterations    : "<< Niter<<std::endl;
    std::cerr << "  Rot search              : "<<rot0<<" "<<rotF<<" "<< rot_step<<std::endl;
    std::cerr << "  Tilt search             : "<<tilt0<<" "<<tiltF<<" "<< tilt_step<<std::endl;
    std::cerr << "  Psi search              : "<<psi0<<" "<<psiF<<" "<< psi_step<<std::endl;
    if (ccf_mode)
      std::cerr << "  -> Use maximum constrained CCF "<<std::endl;
    else {
      std::cerr << "  -> Use maximum-likelihood with a real-space error model "<<std::endl;
      std::cerr << "  Sigma noise             : "<<sqrt(sigma_noise2)<<std::endl;
    }
    if (fn_sym!="")
      std::cerr << "  Symmetry file           : "<< fn_sym<<std::endl;
    if (fn_solv!="")
      std::cerr << "  Solvent mask            : "<< fn_solv<<std::endl;
    if (fn_solv2!="")
      std::cerr << "  Second solvent mask     : "<< fn_solv2<<std::endl;
    if (fn_frac!="") {
      std::cerr << "  -> Read initial model fractions from "<< fn_frac<<std::endl;
    }
    if (fix_fractions) {
      std::cerr << "  -> Do not update estimates of model fractions."<<std::endl;
    }
    if (fix_sigma_noise) {
      std::cerr << "  -> Do not update sigma-estimate of noise."<<std::endl;
    }
    if (fix_sigma_offset) {
      std::cerr << "  -> Do not update sigma-estimate of origin offsets."<<std::endl;
    }
    std::cerr << " -----------------------------------------------------------------"<<std::endl;
  }

} 

// Usage ===================================================================
void Prog_MLalign3D_prm::usage() {
    std::cerr << "Usage:  MLalign3D [options] "<<std::endl;
    std::cerr << "   -i <selfile>                : Selfile with input images \n";
    std::cerr      << "   -ref <selfile/image>        : Selfile with initial reference images/single reference image \n";
    std::cerr      << " [ -o <rootname=\"MLalign3D\"> ] : Output rootname \n";
    std::cerr << " [ -offset <float=0>]          : Search radius for limited translations [pix] \n";
    std::cerr << " [ -max_shift <float=dim/2>]   : Maximum allowed shift [pix] \n";
    std::cerr << " [ -CCF]                       : Use maxCC target (default: max. likelihood) \n";
    std::cerr << "   Angular search parameters:        \n";                          
    std::cerr << " [ -rot  <rot0=0>  <rotF=0>  <step_rot=1> \n";
    std::cerr << " [ -tilt <tilt0=0> <tiltF=0> <step_tilt=1>\n";
    std::cerr << " [ -psi  <psi0=0>  <psiF=0>  <step_psi=1> \n";
    std::cerr << " [ -sym <symfile> ]            : Enforce symmetry \n";
    std::cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_MLalign3D_prm::extended_usage() {
  std::cerr << "Additional options: "<<std::endl;
  std::cerr << " [ -iter <int=100> ]           : Maximum number of iterations to perform \n";
  std::cerr << " [ -frac <docfile=\"\"> ]        : Docfile with expected model fractions (default: even distr.)\n";
  std::cerr << " [ -fix_sigma_noise]           : Do not re-estimate the standard deviation in the pixel noise \n";
  std::cerr << " [ -fix_sigma_offset]          : Do not re-estimate the standard deviation in the offsets \n";
  std::cerr << " [ -fix_fractions]             : Do not re-estimate the model fractions \n";
  std::cerr << std::endl;
  exit(1);
}

// Read all input selfiles in memory 
void Prog_MLalign3D_prm::produce_Side_info() {

  FileName fn_vol, fn_tmp;
  VolumeXmipp vol;
  headerXmipp head;
  Matrix1D<double> offsets(2);
  Matrix3D<double> Maux;
  int xdim,ydim,zdim,c,iaux,ifound,refno=0;
  float xx,yy;
  double dum,avg;

  nr_trans=0;
  nr_tilt=0;
  nr_rot_tilt=0;
  for (double tilt=tilt0; tilt<=tiltF; tilt+=tilt_step) {
    nr_tilt++;
    if (tilt==0) nr_rot_tilt+=1;
    else nr_rot_tilt+=CEIL(360.*sin(DEG2RAD(tilt))/rot_step);
  }
  nr_psi=0;
  for (double psi=psi0; psi<=psiF; psi+=psi_step) nr_psi++;

  // Get image size
  SF.go_beginning();
  vol.read(SF.NextImg());
  if (XSIZE(vol())!=YSIZE(vol())) REPORT_ERROR(1,"ERROR% unequal dimensions: Only cubic volumes are allowed!");
  if (XSIZE(vol())!=ZSIZE(vol())) REPORT_ERROR(1,"ERROR% unequal dimensions: Only cubic volumes are allowed!");
  dim=XSIZE(vol());
  bigdim=CEIL(sqrt(3.)*dim);
  dim3=(double)dim*dim*dim;

  // Make relevant masks 
  mask.resize(dim,dim,dim);
  mask.setXmippOrigin();
  outside_mask.resize(dim,dim,dim);
  outside_mask.setXmippOrigin();
  BinarySphericalMask(mask,dim/2,INNER_MASK);
  BinarySphericalMask(outside_mask,dim/2,OUTSIDE_MASK);

  // Read image- and reference- selfiles
  if (Is_VolumeXmipp(fn_ref)) {
    SFr.reserve(1);
    SFr.insert(fn_ref);
  } else {
    SFr.read(fn_ref);
  }

  // Read in all reference images in memory
  Maux.initZeros(dim,dim,dim);
  Maux.setXmippOrigin();
  nr_ref=0;
  SFr.go_beginning();
  while ((!SFr.eof())) {
    vol.read(SFr.NextImg());
    vol().setXmippOrigin();
    computeStats_within_binary_mask(outside_mask,vol(),dum,dum,avg,dum);
    apply_binary_mask(mask,vol(),vol(),avg);
    Maux+=vol();
    Iref.push_back(vol());
    // Default start is all equal model fractions
    alpha_k.push_back((double)1/SFr.ImgNo());
    nr_ref++;
  }

  // limited implementation of subtract_current
  if (do_subtract_current && nr_ref>1) 
    REPORT_ERROR(1,"ERROR% subtract_current for now only implemented with 1 reference!");
  if (do_subtract_current && !ccf_mode) 
    REPORT_ERROR(1,"ERROR% subtract_current for now only implemented for -CCF!");
  if (fn_mask!="") {
    if (!ccf_mode) REPORT_ERROR(1,"ERROR% mask option for now only implemented for -CCF!");
    else {
      corr_mask.read(fn_mask);
      corr_mask.setXmippOrigin();
    }
  }

  // Prepare for smoothing
  if (theta0>0) theta=theta0;
  if (theta>0) {
    double theta_corr=(double)nr_img/(double)(nr_ref*nr_ref);
    FOR_ALL_MODELS() {
      Iref[refno]+=(theta*theta_corr)*Maux;
      Iref[refno]/=(theta*theta_corr+1);
    }
  }

  // Read in symmetry information
  if (fn_sym!="") SL.read_sym_file(fn_sym);

  // read in model fractions if given on command line (else uniform distribution)
  if (fn_frac!="") {
    DocFile  DF;
    DocLine DL;
    DF.read(fn_frac);
    DF.go_first_data_line();
    double sumfrac=0.;
    for (refno=0; refno<nr_ref; refno++) {
      DL=DF.get_current_line();
      alpha_k[refno]=DL[0];
      sumfrac+=alpha_k[refno];
      DF.next_data_line();
    }
    if (ABS(sumfrac-1.)>1e-4) 
      std::cerr << " ->WARNING: Sum of all expected model fractions ("<<floatToString(sumfrac)<<") is not one!"<<std::endl;
    for (refno=0; refno<nr_ref; refno++) { alpha_k[refno]/=sumfrac; }
  }

  // Store angles for missing wedges
  nr_wedge=0;
  if (fn_wlist!="") {
    wedgelist ww;
    DocFile DF1;
    DF1.read(fn_wlist);
    DF1.go_beginning();
    while (!DF1.eof()) {
      ww.num=ROUND(DF1(0));
      ww.th0=(double)DF1(1);
      ww.thF=(double)DF1(2);
      wedges.push_back(ww);
      nr_wedge++;
      DF1.next();
    }
    DF1.clear();
  }

  // Store tomogram angles, offset vectors and missing wedge parameters
  if (fn_doc!="") {
    DocFile DF;
    DF.read(fn_doc);
    
    SF.go_beginning();
    while (!SF.eof()) {
      fn_vol=SF.NextImg();
      if (DF.search_comment(fn_vol)) {
	img_rot.push_back( DF(0));
	img_tilt.push_back(DF(1));
	img_psi.push_back( DF(2));
	img_xoff.push_back(DF(3));
	img_yoff.push_back(DF(4));
	img_zoff.push_back(DF(5));
	if (nr_wedge>0) {
	  img_wednr.push_back(DF(6));
	  iaux=ROUND(DF(6));
	  ifound=0;
	  for (int iw=0; iw<nr_wedge; iw++) {
	    if ( iaux==wedges[iw].num) {
	      img_th0.push_back(wedges[iw].th0);
	      img_thF.push_back(wedges[iw].thF);
	      ifound++;
	    }
	  }
	  if (ifound!=1) {
	    std::cerr <<ifound;
	    std::cerr << "ERROR% wedge "<<iaux<<" for tomogram "<<fn_vol<<" not found in wedge-list"<<std::endl;
	    exit(0);
	  }
	} else {
	  img_wednr.push_back(0.);
	  img_th0.push_back(0.);
	  img_thF.push_back(0.);
	}
      } else {
	std::cerr << "ERROR% "<<fn_vol<<" not found in document file"<<std::endl;
	exit(0);
      }
    }
  } else {
    SF.go_beginning();
    while (!SF.eof()) {
      SF.NextImg();
      img_rot.push_back( 0.);
      img_tilt.push_back(0.);
      img_psi.push_back( 0.);
      img_xoff.push_back(0.);
      img_yoff.push_back(0.);
      img_zoff.push_back(0.);
      img_wednr.push_back(0.);
      img_th0.push_back( 0.);
      img_thF.push_back( 0.);
    }
  }

}

// Calculate probability density function of the translations
void Prog_MLalign3D_prm::calculate_pdf_trans() {

  double r2,pdfpix,sum;
  pdf_trans.resize(dim,dim,dim);
  pdf_trans.setXmippOrigin();
  Mr2.resize(dim,dim,dim);
  Mr2.setXmippOrigin();

  FOR_ALL_ELEMENTS_IN_MATRIX3D(pdf_trans) {
    r2=(double)(j*j + i*i + k*k);
    if (sigma_offset>0.) {
      pdfpix=exp(-r2/(2*sigma_offset*sigma_offset));
      pdfpix/=(double)2*PI*sigma_offset*sigma_offset*nr_rot_tilt*nr_psi;
    } else {
      if (k==0 && j==0 && i==0) pdfpix=1.;
      else pdfpix=0.;
    }
    VOL_ELEM(pdf_trans,k,i,j)=pdfpix;
    VOL_ELEM(Mr2,k,i,j)=r2;
  }

}

void Prog_MLalign3D_prm::CCF_integrate(
	  Matrix3D<double> &Mimg, Matrix2D<double> &A_img, 
	  std::vector<Matrix3D<double> > &wsum_Mimgs,
	  std::vector<Matrix3D<double> > &wsum_Mwedge,
	  double &th0, double &thF, std::vector<double> &sumw, double &maxccf, 
	  int &opt_refno, double &opt_rot, double &opt_tilt, double &opt_psi,
	  double &opt_xoff, double &opt_yoff, double &opt_zoff) {


  Matrix3D<double> Maux, Mwedge, Mwedgebig, Msmall, Maux2;
  Matrix3D<std::complex<double> > Fimg, Faux, FFref;
  Matrix2D<double> A(4,4), A_rot(4,4), I(4,4);
  Matrix1D<double> opt_offsets(3);
  std::vector<double> refw(nr_ref);
  double aux,weight,diff,dum,avg;
  double sum_refw=0.;
  int ioptx,iopty,ioptz,xmax,ymax,zmax;
  double rot,tilt,psi,rot_sam,Xi2,A2,corrA2,max;

  maxccf=-999.;
  Maux.resize(dim,dim,dim);
  Maux2.resize(dim,dim,dim);
  Faux.resize(dim,dim,dim);
  Maux.setXmippOrigin();
  Maux2.setXmippOrigin();
  I.initIdentity();
  Mwedgebig.resize(bigdim,bigdim,bigdim);
  Mwedgebig.setXmippOrigin();
  if (nr_wedge>0) {
    Mwedge.resize(dim,dim,dim);
    Mwedge.setXmippOrigin();
  } else Mwedgebig.initConstant(1.);

  // Fourier Transform of the experimental image
  if (nr_wedge>0 || sigma_offset>0) {
    FourierTransform(Mimg,Fimg);
    Fimg*=dim3;
  }

  // Calculate Xi2 forstd::cing the supposed missing wedge
  if (nr_wedge>0) {
    BinaryWedgeMask(Mwedge,th0,thF,I);
    CenterFFT(Mwedge,false);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Fimg) {
      dVkij(Faux,k,i,j)=dVkij(Fimg,k,i,j)*(dVkij(Mwedge,k,i,j)/dim3);
    }
    InverseFourierTransform(Faux,Mimg);
    Mimg.setXmippOrigin();
  }
  Xi2=Mimg.sum2();

  if (do_subtract_current) {
    //Pre-orient current volume
    computeStats_within_binary_mask(outside_mask,Mimg,dum,dum,avg,dum);
    applyGeometryBSpline(Maux,A_img,Mimg,3,IS_NOT_INV,DONT_WRAP,avg);
    FOR_ALL_MODELS() {
      dum=1./nr_img;
      aux=nr_img/(nr_img-1.);
      Iold.push_back(Iref[refno]);
      Iref[refno]-=dum*Maux;
      Iref[refno]*=aux;
    }
  }

  // Calculate all CCFs
  FOR_ALL_TILT() {
    tilt=tilt0+(double)itilt*tilt_step; 
    // Make an even distribution of rot-tilt angles
    if (tilt==0) nr_rot=1;
    else nr_rot=CEIL(360.*sin(DEG2RAD(tilt))/rot_step);
    rot_sam=360./(double)nr_rot;
    FOR_ALL_ROT() {
      rot=rot0+(double)irot*rot_sam;
      FOR_ALL_PSI() {
	// local psi-search around -rot
	psi=psi0+(double)ipsi*psi_step; 
	A_rot=Euler_rotation3DMatrix(rot,tilt,psi);
	A=A_rot*A_img;
	FOR_ALL_MODELS() {
	  refw[refno]=0.;
	  computeStats_within_binary_mask(outside_mask,Iref[refno],dum,dum,avg,dum);
	  applyGeometryBSpline(Maux,A,Iref[refno],3,IS_INV,DONT_WRAP,avg);
	  // Only calculate cross-correlation in a limited mask region
	  if (fn_mask!="") {
	    applyGeometryBSpline(Maux2,A,corr_mask,3,IS_INV,DONT_WRAP,0.);
	    Maux*=Maux2;
	  }
	  if (nr_wedge>0 || sigma_offset>0) FourierTransform(Maux,FFref);
	  // compute A2
	  if (nr_wedge>0) {
	    A2=corrA2=0.;
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FFref) {
	      aux=dVkij(FFref,k,i,j).real()*dVkij(FFref,k,i,j).real()+
		dVkij(FFref,k,i,j).imag()*dVkij(FFref,k,i,j).imag();
	      corrA2+=aux;
	      if (dVkij(Mwedge,k,i,j)>0.) A2+=aux;
	    }
	    A2*=Maux.sum2()/(corrA2);
	  } else A2=Maux.sum2();
	  // compute CCFs
	  if (sigma_offset==0) {
	    // in direct space
	    if (nr_wedge>0) {
	      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FFref) {
		dVkij(FFref,k,i,j)*=dVkij(Mwedge,k,i,j);
	      }
	      InverseFourierTransform(FFref,Maux);
	      Maux.setXmippOrigin();
	    }
	    max=0.;
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Maux) {
	      max+=dVkij(Maux,k,i,j)*dVkij(Mimg,k,i,j);
	    }
	    max/=sqrt(A2)*sqrt(Xi2);
	    xmax=ymax=zmax=0;
	  } else {
	    // or in Fourier space
	    FFref*=dim3;
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FFref) {
	      dVkij(Faux,k,i,j)=dVkij(Fimg,k,i,j)*conj(dVkij(FFref,k,i,j));
	      if (nr_wedge>0) dVkij(Faux,k,i,j)*=dVkij(Mwedge,k,i,j);
	    }
	    InverseFourierTransform(Faux,Maux);
	    Maux/=dim3*sqrt(A2)*sqrt(Xi2);
	    Maux.setXmippOrigin();
	    CenterFFT(Maux,true);
	    Maux.maxIndex(zmax,ymax,xmax);
	    max=VOL_ELEM(Maux,zmax,ymax,xmax);
	  }
	  if (verb>0) std::cout <<rot<<" "<<tilt<<" "<<psi<<" "<<xmax<<" "<<ymax<<" "<<zmax<<" "<<max<<std::endl;
	  if (max>maxccf) {
	    maxccf=max;
	    opt_xoff=(double)-xmax;
	    opt_yoff=(double)-ymax;
	    opt_zoff=(double)-zmax;
	    opt_refno=refno;
	    opt_rot=rot;
	    opt_tilt=tilt;
	    opt_psi=psi;
	  }
	}
      }
    }
  }

  // If shift acceptable, add image and wedge to the sums
  double r2=opt_xoff*opt_xoff+opt_yoff*opt_yoff+opt_zoff*opt_zoff;
  // (for subtract_current: always accept shift!)
  if ( r2 <= max_shift*max_shift || do_subtract_current ) {
    opt_offsets(0)=opt_xoff;
    opt_offsets(1)=opt_yoff;
    opt_offsets(2)=opt_zoff;
    Mimg.selfTranslate(opt_offsets,WRAP);
    A_rot=Euler_rotation3DMatrix(opt_rot,opt_tilt,opt_psi);
    A=A_rot*A_img;
    computeStats_within_binary_mask(outside_mask,Mimg,dum,dum,avg,dum);
    applyGeometryBSpline(Maux,A,Mimg,3,IS_NOT_INV,DONT_WRAP,avg);
    Maux.setXmippOrigin();
    wsum_Mimgs[opt_refno]+=Maux;
    sumw[opt_refno]+=1.;
    if (nr_wedge>0) {
      BinaryWedgeMask(Mwedgebig,th0,thF,A);
      CenterFFT(Mwedgebig,false);
    }
    wsum_Mwedge[opt_refno]+=Mwedgebig;
    // optimal transformation
    A.resize(3,3);
    Euler_matrix2angles(A,opt_rot,opt_tilt,opt_psi);
  } else {
    opt_xoff=opt_yoff=opt_zoff=0.;
    opt_rot=opt_tilt=opt_psi=maxccf=0.;
    opt_refno=0;
  }

  if (do_subtract_current) {
    FOR_ALL_MODELS() {
      Iref[refno]=Iold[refno];
    }
    Iold.clear();
  }

}

// Maximum Likelihood calculation in real space for one image ======================
// Integration over all model,orientations and translations
void Prog_MLalign3D_prm::ML_integrate(
	  Matrix3D<double> &Mimg, Matrix2D<double> &A_img, 
	  std::vector<Matrix3D<double> > &wsum_Mimgs,
	  std::vector<Matrix3D<double> > &wsum_Mwedge,
	  double &wsum_sigma_noise2,  double &wsum_sigma_offset,  
	  double &th0, double &thF, 
	  std::vector<double> &sumw, double &LL, double &fracweight, 
	  int &opt_refno, double &opt_rot, double &opt_tilt, double &opt_psi,
	  double &opt_xoff, double &opt_yoff, double &opt_zoff) {

  Matrix3D<double> Maux, Mwedge, Msmall, Mwedgebig, Mimgori;
  std::vector<Matrix3D<double> > Mweight;
  Matrix3D<std::complex<double> > Fimg, Faux, FFref;
  Matrix2D<double> A(4,4), A_rot(4,4), I(4,4);
  Matrix1D<double> offsets(3);
  std::vector<double> refw(nr_ref);
  double sum_refw=0.,maxweight=-99.e99,mindiff2=99.e99,wsum_corr=0.;
  int iweight,ioptx,iopty,ioptz,xmax,ymax,zmax,smalldim;
  double rot,tilt,psi,rot_sam,Xi2,A2,corrA2,mind,sumweight,maxw,sum,aux,weight,dum,avg;

  Maux.resize(dim,dim,dim);
  Faux.resize(dim,dim,dim);
  Maux.setXmippOrigin();
  I.initIdentity();
  fracweight=0.;
  Mwedgebig.resize(bigdim,bigdim,bigdim);
  if (nr_wedge>0) {
    Mwedge.resize(dim,dim,dim);
    Mwedge.setXmippOrigin();
  } else Mwedgebig.initConstant(1.);
  Mwedgebig.setXmippOrigin();

  // Smaller matrices to store weights to save memory
  Mweight.resize(nr_ref*nr_rot_tilt*nr_psi);
  smalldim=6*CEIL(sigma_offset)+1;
  smalldim=MIN(dim,smalldim);
  smalldim=MAX(1,smalldim);
  if (!(smalldim%2)) smalldim--; //keep impair number!
  Msmall.resize(smalldim,smalldim,smalldim);
  Msmall.setXmippOrigin();

  // Fourier Transform of the experimental image
  if (nr_wedge>0 || sigma_offset>0) {
    FourierTransform(Mimg,Fimg);
    Fimg*=dim3;
  }

  // Calculate Xi2 forstd::cing the supposed missing wedge
  if (nr_wedge>0) {
    BinaryWedgeMask(Mwedge,th0,thF,I);
    CenterFFT(Mwedge,false);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Fimg) {
      dVkij(Faux,k,i,j)=dVkij(Fimg,k,i,j)*(dVkij(Mwedge,k,i,j)/dim3);
    }
    InverseFourierTransform(Faux,Mimg);
    Mimg.setXmippOrigin();
  }
  Xi2=Mimg.sum2();

  // Calculate all squared differences
  iweight=-1;
  FOR_ALL_TILT() {
    tilt=tilt0+(double)itilt*tilt_step; 
    // Make an even distribution of rot-tilt angles
    if (tilt==0) nr_rot=1;
    else nr_rot=CEIL(360.*sin(DEG2RAD(tilt))/rot_step);
    rot_sam=360./(double)nr_rot;
    FOR_ALL_ROT() {
      rot=rot0+(double)irot*rot_sam;
      FOR_ALL_PSI() {
	psi=psi0+(double)ipsi*psi_step; 
	A_rot=Euler_rotation3DMatrix(rot,tilt,psi);
	A=A_rot*A_img;
	FOR_ALL_MODELS() {
	  iweight++;
	  refw[refno]=0.;
	  computeStats_within_binary_mask(outside_mask,Iref[refno],dum,dum,avg,dum);
	  applyGeometryBSpline(Maux,A,Iref[refno],3,IS_INV,DONT_WRAP,avg);
	  if (nr_wedge>0 || sigma_offset>0) FourierTransform(Maux,FFref);
	  if (sigma_offset==0) {
	    if (nr_wedge>0) {
	      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FFref) {
		dVkij(FFref,k,i,j)*=dVkij(Mwedge,k,i,j);
	      }
	      InverseFourierTransform(FFref,Maux);
	      Maux.setXmippOrigin();
	    }
	    Maux-=Mimg;
	    mind=0.5*Maux.sum2();
	    dVkij(Msmall,0,0,0)=mind;
	  } else {
	    // Calculate A2
	    if (nr_wedge>0) {
	      A2=corrA2=0.;
	      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FFref) {
		aux=dVkij(FFref,k,i,j).real()*dVkij(FFref,k,i,j).real()+
		  dVkij(FFref,k,i,j).imag()*dVkij(FFref,k,i,j).imag();
		corrA2+=aux;
		A2+=aux*dVkij(Mwedge,k,i,j);
	      }
	      A2*=Maux.sum2()/(corrA2);
	    } else A2=Maux.sum2();
	    // Calculate correlation Matrix
	    FFref*=dim3;
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FFref) {
	      dVkij(Faux,k,i,j)=dVkij(Fimg,k,i,j)*conj(dVkij(FFref,k,i,j));
	      if (nr_wedge>0) dVkij(Faux,k,i,j)*=dVkij(Mwedge,k,i,j);
	    }
	    InverseFourierTransform(Faux,Maux);
	    Maux/=dim3;
	    Maux.setXmippOrigin();
	    CenterFFT(Maux,true);
	    FOR_ALL_ELEMENTS_IN_MATRIX3D(Msmall) {
	      VOL_ELEM(Msmall,k,i,j)=0.5*(A2+Xi2)-VOL_ELEM(Maux,k,i,j);
	    }
	    mind=Msmall.computeMin();
	  }
	  if (mind<mindiff2) mindiff2=mind;
	  Mweight[iweight]=Msmall;
	}
      }
    }
  }

  // Now that we have mindiff2 calculate all weights and maxweight
  iweight=-1;
  FOR_ALL_TILT() {
    tilt=tilt0+(double)itilt*tilt_step; 
    if (tilt==0) nr_rot=1;
    else nr_rot=CEIL(360.*sin(DEG2RAD(tilt))/rot_step);
    rot_sam=360./(double)nr_rot;
    FOR_ALL_ROT() {
      rot=rot0+(double)irot*rot_sam;
      FOR_ALL_PSI() {
	psi=psi0+(double)ipsi*psi_step; 
	FOR_ALL_MODELS() {
	  iweight++;
	  if (alpha_k[refno]>0.) {
	    sum=0.;
	    FOR_ALL_ELEMENTS_IN_MATRIX3D(Mweight[iweight]) {
	      aux=(VOL_ELEM(Mweight[iweight],k,i,j)-mindiff2)/(sigma_noise2);
	      // next line because of numerical precision of exp-function
	      if (aux>1000.) aux=0.;
	      else aux=exp(-aux)*alpha_k[refno]*VOL_ELEM(pdf_trans,k,i,j);
	      wsum_corr+=aux*VOL_ELEM(Mweight[iweight],k,i,j);
	      VOL_ELEM(Mweight[iweight],k,i,j)=aux;
	      sum+=aux;
	    }
	    refw[refno]+=sum;
	    sum_refw+=sum;
	    maxw=Mweight[iweight].computeMax();
	    if (maxw>maxweight) {
	      maxweight=maxw;
	      Mweight[iweight].maxIndex(ioptz,iopty,ioptx);
	      opt_refno=refno;
	      opt_rot=rot;
	      opt_tilt=tilt;
	      opt_psi=psi;
	    }
	  } else Mweight[iweight].initZeros();
	}
      }
    }
  }

  // Accumulate weighted sums of images, sigma2 and fraction-parameters
  // and normalize them by sum_refw, such that sum over all weights is one!
  FOR_ALL_MODELS() { sumw[refno]+=refw[refno]/sum_refw; }
  fracweight=maxweight/sum_refw; 
  wsum_sigma_noise2+=(2*wsum_corr/sum_refw);
  iweight=-1;
  FOR_ALL_TILT() {
    tilt=tilt0+(double)itilt*tilt_step; 
    if (tilt==0) nr_rot=1;
    else nr_rot=CEIL(360.*sin(DEG2RAD(tilt))/rot_step);
    rot_sam=360./(double)nr_rot;
    FOR_ALL_ROT() {
      rot=rot0+(double)irot*rot_sam;
      FOR_ALL_PSI() {
	// local psi-search around -rot
	psi=psi0+(double)ipsi*psi_step; 
	A_rot=Euler_rotation3DMatrix(rot,tilt,psi);
	A=A_rot*A_img;
	FOR_ALL_MODELS() {
	  iweight++;
	  sumweight=Mweight[iweight].sum();
	  if (sumweight>SIGNIFICANT_WEIGHT_LOW*maxweight) {
	    Mweight[iweight]/=sum_refw;
	    sumweight/=sum_refw;
	    if (sigma_offset==0) {
	      //computeStats_within_binary_mask(outside_mask,Mimg,dum,dum,avg,dum);
	      Maux.resize(dim,dim,dim); 
	      Maux.setXmippOrigin();
	      applyGeometryBSpline(Maux,A,Mimg,3,IS_NOT_INV,DONT_WRAP);
	      Maux*=sumweight;
	    } else {
	      // Use Maux, because Mweight is smaller than dim x dim x dim!
	      Maux.initZeros(dim,dim,dim);
	      Maux.setXmippOrigin();
	      FOR_ALL_ELEMENTS_IN_MATRIX3D(Mweight[iweight]) {
		weight=VOL_ELEM(Mweight[iweight],k,i,j);
		VOL_ELEM(Maux,k,i,j)=weight;
		wsum_sigma_offset+=weight*VOL_ELEM(Mr2,k,i,j);
	      }
	      FourierTransform(Maux,Faux);
	      Faux*=dim3;
	      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux) {
		dVkij(FFref,k,i,j)=conj(dVkij(Faux,k,i,j))*dVkij(Fimg,k,i,j);
	      }
	      InverseFourierTransform(FFref,Maux);
	      Maux/=dim3;
	      CenterFFT(Maux,true);
	      //computeStats_within_binary_mask(outside_mask,Maux,dum,dum,avg,dum);
	      // somehow providing avg on the next line goes MUCH worse...
	      applyGeometryBSpline(Maux,A,Maux,3,IS_NOT_INV,DONT_WRAP);
	      Maux.setXmippOrigin();
	    }
	    wsum_Mimgs[refno]+=Maux;
	    if (nr_wedge>0) {
	      BinaryWedgeMask(Mwedgebig,th0,thF,A);
	      CenterFFT(Mwedgebig,false);
	    }
	    Maux=sumweight*Mwedgebig;
	    wsum_Mwedge[refno]+=Maux;
	  }
	}
      }
    }
  }

  // Save the optimal (combined) angles and origin offsets
  A_rot=Euler_rotation3DMatrix(opt_rot,opt_tilt,opt_psi);
  A=A_rot*A_img;
  A.resize(3,3);
  Euler_matrix2angles(A,opt_rot,opt_tilt,opt_psi);
  opt_xoff=-(double)ioptx;
  opt_yoff=-(double)iopty;
  opt_zoff=-(double)ioptz;

  // Output the (log-likelihood) target function value
  LL+= log(sum_refw) - mindiff2/sigma_noise2 - dim3*log(2.50663*sqrt(sigma_noise2));
  Mweight.clear();

}

void Prog_MLalign3D_prm::ML_sum_over_all_images(SelFile &SF, std::vector<Matrix3D<double> > &Iref, 
			  double &LL, double &sumcorr, DocFile &DFo, 
			  std::vector<Matrix3D<double> > &wsum_Mref,
			  std::vector<Matrix3D<double> > &wsum_Mwedge, 
			  double &wsum_sigma_noise2, 
			  double &wsum_sigma_offset, std::vector<double> &sumw) {

  std::vector <Matrix3D<std::complex<double> > > Fref,wsum_Fimgs;
  Matrix3D<std::complex<double> > Fdzero;
  Matrix3D<double> Mdzero, Mdzerobig;
  Matrix1D<double> dataline(10), opt_offsets(3), mis_offsets(3);  
  Matrix2D<double> A_img(4,4);
  VolumeXmipp      img;
  FileName         fn_img;
  DocFile          DFmis;
  double           th0,thF,opt_rot,opt_tilt,opt_psi,maxcorr;
  double           opt_xoff,opt_yoff,opt_zoff;
  double           mis_rot,mis_tilt,mis_psi;
  int              nn,imgno,opt_refno;
  SelLine          line;

  // Initialize
  Mdzero.resize(dim,dim,dim);
  Mdzero.setXmippOrigin();
  Mdzerobig.resize(bigdim,bigdim,bigdim);
  Mdzerobig.setXmippOrigin();
  Fdzero.resize(dim,dim,dim);
  wsum_Fimgs.clear();
  sumw.clear();
  LL=0.;
  wsum_sigma_noise2=0.;
  wsum_sigma_offset=0.;
  wsum_Mwedge.clear();
  wsum_Mref.clear();
  sumcorr=0.;
  FOR_ALL_MODELS() {
    wsum_Mwedge.push_back(Mdzerobig);
    wsum_Mref.push_back(Mdzero);
    sumw.push_back(0.);
  }

  // Pre-calculate the prior of the origin offsets
  if (!ccf_mode) calculate_pdf_trans();

  // Read in document file with misalignment parameters
  if (fn_misalign!="") DFmis.read(fn_misalign);

  // Loop over all images
  nn=SF.ImgNo();
  if (verb>0) init_progress_bar(nn);
  imgno=0;
  SF.go_beginning();
  while ((!SF.eof())) {

    // get tomogram and geometric information
    fn_img=SF.NextImg();
    img.read(fn_img);
    img().setXmippOrigin();
    A_img=Euler_rotation3DMatrix(img_rot[imgno],img_tilt[imgno],img_psi[imgno]);
    opt_offsets(0)=ROUND(img_xoff[imgno]);
    opt_offsets(1)=ROUND(img_yoff[imgno]);
    opt_offsets(2)=ROUND(img_zoff[imgno]);
    if (nr_wedge>0) {
      th0=img_th0[imgno];
      thF=img_thF[imgno];
    } else th0=thF=0.;

    // planned misalignment 
    if (fn_misalign!="") misalign(A_img,DFmis,fn_img);
    // apply (wrapped around!) integer translation of pre-orientation
    img().selfTranslate(opt_offsets,WRAP);

    // Perform integration over all references, rotations and translations
    if (ccf_mode) {
      CCF_integrate(img(),A_img,wsum_Mref,wsum_Mwedge,
		    th0,thF,sumw,maxcorr,opt_refno,opt_rot,opt_tilt,opt_psi,
		    opt_xoff,opt_yoff,opt_zoff);

    } else {
      ML_integrate(img(),A_img,wsum_Mref,wsum_Mwedge,
		   wsum_sigma_noise2,wsum_sigma_offset,th0,thF,
		   sumw,LL,maxcorr,opt_refno,opt_rot,opt_tilt,opt_psi,
		   opt_xoff,opt_yoff,opt_zoff);

    }

    // Update translations in the docfile
    opt_xoff+=img_xoff[imgno];
    opt_yoff+=img_yoff[imgno];
    opt_zoff+=img_zoff[imgno];

    // Output to docfile
    sumcorr+=maxcorr;
    dataline(0)=opt_rot;                 // rot
    dataline(1)=opt_tilt;                // tilt
    dataline(2)=opt_psi;                 // psi
    dataline(3)=opt_xoff;                // Xoff
    dataline(4)=opt_yoff;                // Yoff
    dataline(5)=opt_zoff;                // Zoff
    if (nr_wedge>0) dataline(6)=img_wednr[imgno];        // missing wedge number
    else dataline(6)=0.;
    dataline(7)=(double)(opt_refno+1);   // Ref
    dataline(8)=maxcorr;                 // P_max/P_tot or Corr
    DFo.append_comment(img.name());
    DFo.append_data_line(dataline);

    if (do_subtract_current && ccf_mode) {
      img_xoff[imgno]=opt_xoff;
      img_yoff[imgno]=opt_yoff;
      img_zoff[imgno]=opt_zoff;
      img_rot[imgno]=opt_rot;
      img_tilt[imgno]=opt_tilt;
      img_psi[imgno]=opt_psi;
    }

    imgno++;
    if (verb>0) progress_bar(imgno);

  }
  if (verb>0) progress_bar(nn);


}

// Update all model parameters
void Prog_MLalign3D_prm::update_parameters(std::vector<Matrix3D<double> > &wsum_Mref,
                                           std::vector<Matrix3D<double> > &wsum_Mwedge, 
					   double &wsum_sigma_noise2, double &wsum_sigma_offset,
					   std::vector<double> &sumw, 
					   double &sumcorr, double &sumw_allrefs, int iter) {

  VolumeXmipp Vaux,Vaux2;
  Matrix1D<double> rmean_sigma2,rmean_wedge;
  Matrix1D<int> center(3),radial_count;
  Matrix3D<std::complex<double> > Faux,Fsum;
  Matrix3D<double> Maux,Msum,allsum_Mwedge;
  double rr,dum,avg,theta_corr,sum_ref=0.;

  // Pre-calculate sumw_allrefs
  allsum_Mwedge.initZeros(dim,dim,dim);
  allsum_Mwedge.setXmippOrigin();
  Msum.initZeros(dim,dim,dim);
  Msum.setXmippOrigin();
  sumw_allrefs=0.;
  FOR_ALL_MODELS() {
    sumw_allrefs+=sumw[refno];
  }
  sumcorr/=sumw_allrefs;

  // Symmetrize weighted sums
  if (fn_sym!="") {
    FOR_ALL_MODELS() {
      Vaux()=wsum_Mwedge[refno];
      CenterFFT(Vaux(),true);
      symmetrize_Bspline(SL,Vaux,Vaux2,3,false,false);
      CenterFFT(Vaux2(),false);
      wsum_Mwedge[refno]=Vaux2();
      Vaux()=wsum_Mref[refno];
      symmetrize_Bspline(SL,Vaux,Vaux2,3,false,true);
      wsum_Mref[refno]=Vaux2();
    }
  }

  // Now resize the big wsum_Mwedges
  FOR_ALL_MODELS() {
    CenterFFT(wsum_Mwedge[refno],true);
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Msum) {
      VOL_ELEM(Msum,k,i,j)=VOL_ELEM(wsum_Mwedge[refno],k,i,j);
    }
    wsum_Mwedge[refno]=Msum;
    CenterFFT(wsum_Mwedge[refno],false);
  }

  // Calculate reference images from weighted_sums and sum_weights
  Msum.initZeros();
  FOR_ALL_MODELS() {
    if (sumw[refno]>0.) {
      Iref[refno]=wsum_Mref[refno];
      FourierTransform(Iref[refno],Faux);
      allsum_Mwedge+=wsum_Mwedge[refno];
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux) {
	// sometimes, dividing by a number smaller than 1e-3 gives
	// strong artefacts...
	if (dVkij(wsum_Mwedge[refno],k,i,j)>1e-3) 
	  dVkij(Faux,k,i,j)/=dVkij(wsum_Mwedge[refno],k,i,j);
	else dVkij(Faux,k,i,j)=0.;
      }
      InverseFourierTransform(Faux,Iref[refno]);
      Msum+=Iref[refno];
      sum_ref+=1.;
    } else Iref[refno].initZeros();
  }

  // Smooth references...
  if (theta>0) {
    theta_corr=(double)nr_img/(double)(nr_ref*nr_ref);
    FourierTransform(Msum,Fsum);
    Fsum*=theta_corr*theta;
    FOR_ALL_MODELS() {
      FourierTransform(wsum_Mref[refno],Faux);
      Faux+=Fsum;
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux) {
	dVkij(Faux,k,i,j)/=dVkij(wsum_Mwedge[refno],k,i,j)+(double)sum_ref*theta*theta_corr;
      }
      InverseFourierTransform(Faux,Iref[refno]);
    }
    if (theta0>0) 
      theta=theta0*exp(-theta_step*iter);
    else {
      theta-=theta_step;
    }
    theta=MAX(0,theta);
  }

  // Post-processing the references: masking and solvent flattening
  FOR_ALL_MODELS() {
    if (fn_solv!="") solvent_flattening(fn_solv);
    if (fn_solv2!="") solvent_flattening(fn_solv2);
  }

  if (!fix_fractions) 
    FOR_ALL_MODELS() {alpha_k[refno]=sumw[refno]/sumw_allrefs;}

  if (!fix_sigma_offset) 
    sigma_offset=sqrt(wsum_sigma_offset/(3*sumw_allrefs));

  if (!fix_sigma_noise)  {
    if (theta==0) sigma_noise2=wsum_sigma_noise2/(sumw_allrefs*dim3);
    else {
      double diff2_ref=0.;
      FOR_ALL_MODELS() {
	for (int refno2=0; refno2<nr_ref; refno2++) {
	  if (refno!=refno2) {
	    Maux=Iref[refno]-Iref[refno2];
	    //apply_binary_mask(mask,Maux,Maux);
	    diff2_ref+=Maux.sum2();
	  }
	}
      }
      diff2_ref*=theta*theta_corr;
      sigma_noise2=(wsum_sigma_noise2+diff2_ref)/(sumw_allrefs*dim3);
    }
  }

}

// Modify reference volumes ======================================================
void Prog_MLalign3D_prm::solvent_flattening(FileName &fn_solvent) {

  VolumeXmipp solv;
  double solvavg,sumsolv;

  solv.read(fn_solvent);
  solv().setXmippOrigin();

  solvavg=0.,sumsolv=0.;
  FOR_ALL_MODELS() {
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(solv()) {
      solvavg+=dVkij(mask,k,i,j)*dVkij(solv(),k,i,j)*dVkij(Iref[refno],k,i,j);
      sumsolv+=dVkij(mask,k,i,j)*dVkij(solv(),k,i,j);
    }
  }
  solvavg/=sumsolv;
  FOR_ALL_MODELS() {
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(solv()) {
      dVkij(Iref[refno],k,i,j)-=dVkij(solv(),k,i,j)*(dVkij(Iref[refno],k,i,j)-solvavg);
    }
  }

}

void Prog_MLalign3D_prm::misalign(Matrix2D<double> &A_img, 
				  DocFile DFmis, FileName fn_img) {

  Matrix1D<double> mis_offsets(3);
  double           mis_rot,mis_tilt,mis_psi;
  Matrix2D<double> A(4,4),B(4,4);

  if (DFmis.search_comment(fn_img)) {
    mis_rot=DFmis(0);
    mis_tilt=DFmis(1);
    mis_psi=DFmis(2);
    mis_offsets(0)=DFmis(3);
    mis_offsets(1)=DFmis(4);
    mis_offsets(2)=DFmis(5);
  } else {
    std::cerr << "ERROR% "<<fn_img<<" not found in misalignment document file"<<std::endl;
    exit(0);
  }
  A=Euler_rotation3DMatrix(mis_rot,mis_tilt,mis_psi);
  A_img=A*A_img;
  A=translation3DMatrix(mis_offsets);
  A_img=A*A_img;

}

// Caculate symmetrized tomogram or missing wedge
void Prog_MLalign3D_prm::symmetrize_tomogram(Matrix3D<double> &Min, Matrix3D<double> &Mout, 
					     Matrix3D<double> &Mwedge, 
					     SymList &SL, Matrix2D<double> A, 
					     double th0, double thF,
					     bool do_inverse) {

  double dum,avg;
  Matrix2D<double> L(4,4), R(4,4); 
  Matrix3D<double> Maux,Mwed;
  Matrix3D<std::complex<double> > Faux;
  Maux.initZeros(dim,dim,dim);
  Maux.setXmippOrigin();
  Mwed.initZeros(dim,dim,dim);
  Mwed.setXmippOrigin();
  Mout.initZeros(dim,dim,dim);
  Mout.setXmippOrigin();
  Mwedge.initZeros(dim,dim,dim);
  Mwedge.setXmippOrigin();
  
  computeStats_within_binary_mask(outside_mask,Min,dum,dum,avg,dum);
  if (do_inverse) {
    R.initIdentity();
    BinaryWedgeMask(Mwedge,th0,thF,R);
    CenterFFT(Mwedge,false);
    Mout=Min;
  } else {
    BinaryWedgeMask(Mwedge,th0,thF,A);
    CenterFFT(Mwedge,false);
    applyGeometryBSpline(Mout,A,Min,3,IS_NOT_INV,DONT_WRAP,avg);
  }
  
  FourierTransform(Mout,Faux);
  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux) {
    dVkij(Faux,k,i,j)*=dVkij(Mwedge,k,i,j);
  }
  InverseFourierTransform(Faux,Mout);
  
  if (SL.SymsNo()==0) return;
  
  // Symmetrize
  for (int symno=0; symno<SL.SymsNo(); symno++) {
    SL.get_matrices(symno,L,R);
    if (do_inverse) R=A.inv()*R*A;
    else R=R*A;
    BinaryWedgeMask(Mwed,th0,thF,R);
    CenterFFT(Mwed,false);
    Mwedge+=Mwed;
    applyGeometryBSpline(Maux,R,Min,3,IS_NOT_INV,DONT_WRAP,avg);
    FourierTransform(Maux,Faux);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux) {
      dVkij(Faux,k,i,j)*=dVkij(Mwed,k,i,j);
    }
    InverseFourierTransform(Faux,Maux);
    Mout+=Maux;
  }
  //   CenterFFT(Mwedge,false);
  Mwedge/=(double)(SL.SymsNo()+1);
  Mout/=(double)(SL.SymsNo()+1);

}

void Prog_MLalign3D_prm::write_output_files(const int iter, SelFile &SF, DocFile &DF, 
					    double &sumw_allrefs, std::vector<double> &sumw, 
					    double &LL, double &avecorr) {

  FileName fn_tmp,fn_base;
  Matrix1D<double> fracline(1);
  string comment;
  VolumeXmipp tmpvol;
 
  DF.clear();

  fn_base=fn_root;
  if (iter>=0) {
    fn_base+="_it";
    fn_base.compose(fn_base,iter,"");
  }

  // Write out current reference images and fill log-file
  FOR_ALL_MODELS() {
    fn_tmp=fn_base+"_ref";
    fn_tmp.compose(fn_tmp,refno+1,"");
    fn_tmp=fn_tmp+".vol";
    tmpvol()=Iref[refno];
    tmpvol.write(fn_tmp);
    fracline(0)=sumw[refno]/sumw_allrefs;
    DF.insert_comment(fn_tmp);
    DF.insert_data_line(fracline);
  }

  // Write out log-file
  DF.go_beginning();
  comment="MLalign3D-logfile: Number of images= "+floatToString(sumw_allrefs);
  comment+=" LL= "+floatToString(LL,10,5)+" <Pmax/sumP>= "+floatToString(avecorr,10,5);
  DF.insert_comment(comment);
  comment="-noise "+floatToString(sqrt(sigma_noise2),10,7)+" -offset "+floatToString(sigma_offset,10,7)+" -istart "+integerToString(iter+1);
  if (theta>0) comment+=" -theta "+floatToString(theta,6,3)+" -theta_step "+floatToString(theta_step,6,3);
  DF.insert_comment(comment);
  DF.insert_comment("columns: model fraction (1); ");
  fn_tmp=fn_base+".log";
  DF.write(fn_tmp);

}

