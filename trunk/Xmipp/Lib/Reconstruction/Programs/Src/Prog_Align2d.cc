/***************************************************************************
 *
 * Authors:    Sjors Scheres                 (scheres@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
#include "../Prog_Align2d.hh"
#include <Reconstruction/Programs/Prog_FourierFilter.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippFilters.hh>
#include <XmippData/xmippMasks.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFuncs.hh>

// Read arguments ==========================================================
void Prog_align2d_prm::read(int argc, char **argv) _THROW  {

  fn_sel=get_param(argc,argv,"-i");
  SF.read(fn_sel);
  // Filename for output reference
  fn_ave=fn_sel.without_extension()+".xmp";
  fn_ave=fn_ave.insert_before_extension(".ref");
  // Reference image
  fn_ref = get_param(argc,argv,"-ref","");
  // Extension for output images (overwrite input if none)
  oext=get_param(argc,argv,"-oext","");
  // Write out document file?
  fn_doc=get_param(argc,argv,"-doc","");
  // Maximum shift (discard images that shift more in last iteration)
  max_shift=AtoF(get_param(argc,argv,"-max_shift","0"));
  // Maximum rotational change (discard images that rotate more in last iteration)
  max_rot=AtoF(get_param(argc,argv,"-max_rot","0"));
  // Inner and outer radii for rotational correlation
  Ri=AtoI(get_param(argc,argv,"-Ri","0"));
  Ro=AtoI(get_param(argc,argv,"-Ro","0"));
  int dim;
  SF.ImgSize(dim,dim);
  if (Ro==0) Ro=(int)dim/2;
  if (Ro<=Ri) REPORT_ERROR(1,"Align2D: Rout should be larger than Rin");
  // Expected resolution and sampling rate (for filtering)
  if (do_filter=check_param(argc,argv,"-filter")) {
    resol=AtoF(get_param(argc,argv,"-filter"));
    sam=AtoF(get_param(argc,argv,"-sampling"));
  }
  // Number of iterations
  Niter=AtoI(get_param(argc,argv,"-iter","4"));
  // Only translational/rotational
  do_rot=!check_param(argc,argv,"-only_trans");
  do_trans=!check_param(argc,argv,"-only_rot");
  do_acf=check_param(argc,argv,"-ACF");
  do_complete=check_param(argc,argv,"-complete");
  psi_interval=AtoF(get_param(argc,argv,"-psi_step","10"));
  shift_range=AtoI(get_param(argc,argv,"-shift_range","5"));
}

// Show ====================================================================
void Prog_align2d_prm::show() {
    cerr << " Input selfile         : "<<  SF.name()<<endl;
  if (oext!="") 
    cerr << " Output extension      : "<< oext<<endl;
    cerr << " Number of iterations  : "<< Niter<<endl;
  if (fn_ref!="") 
    cerr << " Alignment Reference   : "<<  fn_ref<<endl;
  else 
    cerr << " Alignment Reference   : piramidal combination of images"<< endl;

  if (do_filter) {
    cerr << " Low pass filter [Ang] : "<<  resol<<endl;
    cerr << " Sampling rate   [Ang] : "<<  sam<<endl;
  }
  if (Ri!=0 || Ro!= 0)
    cerr << " Inner radius          : "<< Ri  <<endl;
    cerr << " Outer radius          : "<< Ro  <<endl;
  if (max_shift!=0) 
    cerr << " Max. shift last iter. : "<<max_shift<<endl;
  if (max_rot!=0) 
    cerr << " Max. rotat. last iter.: "<<max_rot<<endl;
  if (fn_doc!="") 
    cerr << " Output document file  : "<<fn_doc<<endl;
  if (!do_rot) 
    cerr << "Skip rotational alignment "<<endl;
  if (!do_trans) 
    cerr << "Skip translational alignment "<<endl;
  if (do_acf) 
    cerr << "Use Auto-Correlation Function alignment" <<endl;
  if (do_complete) { 
    cerr << "Use complete-search alignment with:" <<endl;
    cerr << " Psi interval          : "<< psi_interval  <<endl;
    cerr << " Shift range           : "<< shift_range  <<endl;

  }
}

// usage ===================================================================
void Prog_align2d_prm::usage() {
  cerr << "Usage:  "<<endl;
  cerr << "  align2d [options]"<<endl;
  cerr << "   -i <selfile>             : Selfile containing images to be aligned \n"
       << " [ -ref <image> ]           : reference image; if none: piramidal combination of subset of images \n"
       << " [ -oext <extension> ]      : For output images & selfile; if none: input will be overwritten \n"
       << " [ -iter <int=4> ]          : Number of iterations to perform \n"
       << " [ -Ri <inner radius=0> ]   : Region between radii Ri and Ro will be considered \n"
       << " [ -Ro <outer radius=dim/2>]     for rotational correlation\n"
       << " [ -filter <resol.> ]       : Fourier-filter images to expected resolution [Ang] \n"
       << " [ -sampling <pix. size> ]  : Sampling rate, i.e. pixel size [Ang]; required for filtering. \n"
       << " [ -max_shift <float> ]     : Discard images that shift more in the last iteration [pix]\n"
       << " [ -max_rot <float> ]       : Discard images that rotate more in the last iteration [deg]\n"
       << " [ -doc <docfile> ]         : write output document file with rotations & translations \n"
       << " [ -only_trans ]            : Skip rotational alignment \n"
       << " [ -only_rot ]              : Skip translational alignment \n"
       << " [ -ACF ]                   : Use Auto-Correlation Function alignment \n"
       << " [ -complete ]              : Use complete-search alignment \n"
       << " [ -psi_step <float=10>]    : Sampling interval to search rotation [deg] \n"
       << " [ -shift_range <int=5> ]   : Shift range to be searched (from -range to +range) [pix] \n"

       << endl;
}

// Rotational alignment ========================================================
bool Prog_align2d_prm::align_rot(ImageXmipp &img, const matrix2D<double> &Mref, 
		 const float &max_rot, const float &Rin, const float &Rout) _THROW {

  matrix2D<double> Mimg,Maux,A;
  matrix1D<double> corr;
  matrix2D<int>    mask;
  int nstep;
  int i,i_maxcorr,avewidth;
  double psi_actual,psi_max_coarse,psi_max_fine,sumcorr,psi_coarse_step=15.;
  float psi;

  mask.resize(img().RowNo(),img().ColNo());
  mask.set_Xmipp_origin();
  if (Rout<=Rin)  REPORT_ERROR(1,"Align2d_rot: Rout <= Rin");
  BinaryCrownMask(mask,Rin,Rout,INNER_MASK);

  Mimg.resize(img());
  Mimg.set_Xmipp_origin();
  A=img.get_transformation_matrix();
  apply_geom(Mimg,A,img(),IS_INV,DONT_WRAP,0.);
  Maux.resize(Mimg);
  Maux.set_Xmipp_origin();

  // Optimize correlation in coarse steps
  nstep=(int)(360/psi_coarse_step);
  corr.resize(nstep);
  for (i = 0; i < nstep; i++) {
    psi_actual=(double)i*psi_coarse_step;
    Maux=Mimg.rotate(psi_actual,DONT_WRAP);
    corr(i)=correlation_index(Mref,Maux,&mask);
  }
  corr.max_index(i_maxcorr);
  psi_max_coarse=(double)i_maxcorr*psi_coarse_step;

  // Optimize correlation in fine steps
  nstep=(int)(2*psi_coarse_step-1);
  corr.resize(nstep);
  for (i = 0; i < nstep; i++) {
    psi_actual=psi_max_coarse-psi_coarse_step+1+(double)i;
    Maux=Mimg.rotate(psi_actual,DONT_WRAP);
    corr(i)=correlation_index(Mref,Maux,&mask);
  }
  corr.max_index(i_maxcorr);
  psi_max_fine=psi_max_coarse-psi_coarse_step+(double)i_maxcorr;

  psi=0.;
  sumcorr=0.;
  // Weighted average over neighbours
  avewidth=MIN(i_maxcorr,nstep-i_maxcorr);
  if (avewidth>0) {
    for (i = i_maxcorr-avewidth+1; i < i_maxcorr+avewidth; i++) {
      psi_actual=psi_max_coarse-psi_coarse_step+1+(double)i;
      psi+=corr(i)*psi_actual;
      sumcorr+=corr(i);
    }
    psi/=sumcorr;
  } else psi=psi_max_fine;

  Mimg.core_deallocate();
  Maux.core_deallocate();
  mask.core_deallocate();

  psi=realWRAP(psi,-180.,180.);
  if ((max_rot<XMIPP_EQUAL_ACCURACY) || (ABS(psi) < max_rot)) {
    // Store new rotation in the header of the image 
    // beware: for untilted images (phi+psi) is rotated, and phi can be non-zero!
    // Add new rotation to psi only!
    psi+=img.Psi();
    img.Psi()=realWRAP(psi,0.,360.);   
    return TRUE;
  } else return FALSE;
}

// translational alignment =====================================================
bool Prog_align2d_prm::align_trans(ImageXmipp &img, const matrix2D<double> &Mref,const float &max_shift) _THROW {

  matrix2D<double> Maux,Mcorr,A;
  int              dim,imax,jmax,i_actual,j_actual,dim2;
  double           max,xmax,ymax,sumcorr;
  float            xshift,yshift,shift;
 
  xshift=0.;
  yshift=0.;
  dim=img().RowNo();
  Maux.resize(img());
  Maux.set_Xmipp_origin();
  Mcorr.resize(img());
  Mcorr.set_Xmipp_origin();

  // Apply transformation already present in its header
  A=img.get_transformation_matrix();
  apply_geom(Maux,A,img(),IS_INV,DONT_WRAP,0.);

  // Calculate cross-correlation
  correlation_matrix(Maux,Mref,Mcorr);
  Mcorr.statistics_adjust(0.,1.);
  Mcorr.max_index(imax,jmax);
  max=MAT_ELEM(Mcorr,imax,jmax);

  int              n_max=-1;
  bool             neighbourhood=TRUE;
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
  xmax /= sumcorr; ymax /= sumcorr;
  xshift=(float)-xmax; yshift=(float)-ymax;

  Maux.core_deallocate();
  Mcorr.core_deallocate();

  shift=sqrt(xshift*xshift+yshift*yshift);
  if ((max_shift<XMIPP_EQUAL_ACCURACY) || (shift < max_shift)) {
    // Store shift in the header of the image
    img.Xoff()+=xshift*DIRECT_MAT_ELEM(A,0,0)+yshift*DIRECT_MAT_ELEM(A,0,1);
    img.Yoff()+=xshift*DIRECT_MAT_ELEM(A,1,0)+yshift*DIRECT_MAT_ELEM(A,1,1);
    img.Xoff()=realWRAP(img.Xoff(),(float)-dim/2.,(float)dim/2.);
    img.Yoff()=realWRAP(img.Yoff(),(float)-dim/2.,(float)dim/2.);
    return TRUE;
  } else return FALSE;

}

// ACF-alignment ===============================================================
bool Prog_align2d_prm::align_acf(ImageXmipp &img, const matrix2D<double> &Mref,const float &max_shift, 
                                 const float &max_rot, const float &Rin, const float &Rout) _THROW {

  matrix2D<double> Mimg,Mcorr,Mcorr2,Mimg2,Mref2,Maux,A;
  matrix1D<double> corr;
  matrix2D<int>    mask;
  int              dim,nstep,i,i_maxcorr,avewidth,xdim,ydim,imax,jmax,imax2,jmax2,i_actual,j_actual;
  double           max,max2,xmax,ymax,psi_actual,psi_max_coarse,psi_max_fine,sumcorr,psi_coarse_step=15.;
  float            xshift,yshift,psi,shift;
  bool             OK;

  // Apply transformation already present in header of the image
  Mimg=img();
  Mimg.set_Xmipp_origin();
  A=img.get_transformation_matrix();
  apply_geom(Mimg,A,img(),IS_INV,DONT_WRAP,0.);

  dim=img().RowNo();
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  if (Rout<=Rin)  REPORT_ERROR(1,"Align2d_acf: Rout <= Rin");
  BinaryCrownMask(mask,Rin,Rout,INNER_MASK);
  Maux.resize(img()); 
  Maux.set_Xmipp_origin();
  Mimg2.resize(img()); 
  Mimg2.set_Xmipp_origin();
  Mref2.resize(img()); 
  Mref2.set_Xmipp_origin();
  Mcorr.resize(img());
  Mcorr.set_Xmipp_origin();
  Mcorr2.resize(img());
  Mcorr2.set_Xmipp_origin();

  apply_binary_mask(mask,Mimg,Maux);
  auto_correlation_matrix(Maux,Mimg2);
  apply_binary_mask(mask,Mref,Maux);
  auto_correlation_matrix(Maux,Mref2);
  BinaryCircularMask(mask,Mref2.ColNo(),INNER_MASK);

  // Optimize rotation in coarse steps
  nstep=(int)(360/psi_coarse_step);
  corr.resize(nstep);
  for (i = 0; i < nstep; i++) {
    psi_actual=(double)i*psi_coarse_step;
    Maux=Mimg2.rotate(psi_actual,DONT_WRAP);
    corr(i)=correlation_index(Mref2,Maux,&mask);
  }
  corr.max_index(i_maxcorr);
  psi_max_coarse=(double)i_maxcorr*psi_coarse_step;

  // Optimize rotation in fine steps
  nstep=(int)(2*psi_coarse_step-1);
  corr.resize(nstep);
  for (i = 0; i < nstep; i++) {
    psi_actual=psi_max_coarse-psi_coarse_step+1+(double)i;
    Maux=Mimg2.rotate(psi_actual,DONT_WRAP);
    corr(i)=correlation_index(Mref2,Maux,&mask);
  }
  corr.max_index(i_maxcorr);
  psi_max_fine=psi_max_coarse-psi_coarse_step+(double)i_maxcorr;

  psi=0.;
  sumcorr=0.;
  // Weighted average over neighbours
  avewidth=MIN(i_maxcorr,nstep-i_maxcorr);
  if (avewidth>0) {
    for (i = i_maxcorr-avewidth+1; i < i_maxcorr+avewidth; i++) {
      psi_actual=psi_max_coarse-psi_coarse_step+1+(double)i;
      psi+=corr(i)*psi_actual;
      sumcorr+=corr(i);
    }
    psi/=sumcorr;
  } else psi=psi_max_fine;

  // Calculate cross-correlation
  // Try both psi and psi+180. for CCF

  Maux=Mref.rotate(-psi,DONT_WRAP);
  correlation_matrix(Mimg,Maux,Mcorr);
  Mcorr.max_index(imax,jmax);
  max=MAT_ELEM(Mcorr,imax,jmax);

  Maux=Mref.rotate(-psi+180.,DONT_WRAP);
  correlation_matrix(Mimg,Maux,Mcorr2);
  Mcorr2.max_index(imax2,jmax2);
  max2=MAT_ELEM(Mcorr2,imax2,jmax2);

  if (max2>max) {
    max=max2;
    imax=imax2;
    jmax=jmax2;
    Mcorr=Mcorr2;
    psi=realWRAP(psi+180.,0.,360.);
  }

  // Refine best option
  int              n_max=-1;
  bool             neighbourhood=TRUE;
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
  xmax /= sumcorr; ymax /= sumcorr;
  xshift=(float)-xmax; yshift=(float)-ymax;

  Mimg.core_deallocate();
  Mimg2.core_deallocate();
  Maux.core_deallocate();
  Mref2.core_deallocate();
  Mcorr.core_deallocate();
  Mcorr2.core_deallocate();

  OK=true;
  shift=sqrt(xshift*xshift+yshift*yshift);
  psi=realWRAP(psi,-180.,180.);
  if ((max_shift>XMIPP_EQUAL_ACCURACY) && (shift>max_shift)) OK=false;
  if ((max_rot>XMIPP_EQUAL_ACCURACY) && (ABS(psi)>max_rot)) OK=false;

  if (OK) {
    // Store rotation & translation in the header of the image
    psi+=img.Psi();
    img.Psi()=realWRAP(psi,0.,360.);
    A=img.get_transformation_matrix();
    img.Xoff()+=xshift*DIRECT_MAT_ELEM(A,0,0)+yshift*DIRECT_MAT_ELEM(A,0,1);
    img.Yoff()+=xshift*DIRECT_MAT_ELEM(A,1,0)+yshift*DIRECT_MAT_ELEM(A,1,1);
    img.Xoff()=realWRAP(img.Xoff(),(float)-dim/2.,(float)dim/2.);
    img.Yoff()=realWRAP(img.Yoff(),(float)-dim/2.,(float)dim/2.);
    return TRUE;
  } else return FALSE;

}

// Complete search alignment ========================================================
bool Prog_align2d_prm::align_complete_search(ImageXmipp &img, const matrix2D<double> &Mref, 
                       const float &max_shift, const float &max_rot, const float &psi_interval, 
                       const int &shift_range, const float &Rin, const float &Rout) _THROW{

  matrix2D<double> Mimg,Maux,Mref2,A;
  matrix2D<int>    mask;
  int dim,nstep;
  double psi_actual,corr,maxcorr,xshift,yshift,shift,psi_max;
  bool OK;

  Mref2=Mref;
  Mref2.set_Xmipp_origin();
  Maux.resize(Mimg);
  Maux.set_Xmipp_origin();
  Mimg=img();
  Mimg.set_Xmipp_origin();
  A=img.get_transformation_matrix();
  apply_geom(Mimg,A,img(),IS_INV,DONT_WRAP,0.);

  dim=img().RowNo();
  mask.resize(dim,dim);
  mask.set_Xmipp_origin();
  if (Rout<=Rin)  REPORT_ERROR(1,"Align2d: Rout <= Rin");
  BinaryCrownMask(mask,Rin,Rout,INNER_MASK);

  // Optimize correlation in coarse steps
  nstep=(int)(360/psi_interval);
  maxcorr=0.;
  xshift=0.;
  yshift=0.;
  psi_max=0.;
  for (int i = 0; i < nstep; i++) {
    psi_actual=(double)i*psi_interval;
    Maux=Mimg.rotate(psi_actual,DONT_WRAP);
    for (int x=-shift_range; x <= shift_range; x++) {
      for (int y=-shift_range; y <= shift_range; y++) {   
	corr=correlation(Mref2,Maux,&mask,y,x);
	if (corr>maxcorr) {
	  maxcorr=corr;
	  psi_max=psi_actual;
	  xshift=x;
	  yshift=y;
	}
      }
    }
  }

  OK=true;
  psi_max=realWRAP(psi_max,-180.,180.);
  shift=sqrt(xshift*xshift+yshift*yshift);
  if ((max_shift>XMIPP_EQUAL_ACCURACY) && (shift>max_shift)) OK=false;
  if ((max_rot>XMIPP_EQUAL_ACCURACY) && (ABS(psi_max)>max_rot)) OK=false;

  if (OK) {
    // Store rotation & translation in the header of the image
    psi_max+=img.Psi();
    img.Psi()=realWRAP(psi_max,0.,360.);
    A=img.get_transformation_matrix();
    img.Xoff()+=xshift*DIRECT_MAT_ELEM(A,0,0)+yshift*DIRECT_MAT_ELEM(A,0,1);
    img.Yoff()+=xshift*DIRECT_MAT_ELEM(A,1,0)+yshift*DIRECT_MAT_ELEM(A,1,1);
    img.Xoff()=realWRAP(img.Xoff(),(float)-dim/2.,(float)dim/2.);
    img.Yoff()=realWRAP(img.Yoff(),(float)-dim/2.,(float)dim/2.);
    return TRUE;
  } else return FALSE;

}

// PsPc piramidal combination of images ========================================
void Prog_align2d_prm::do_pspc() _THROW {

  int               imgno,nlev,n_piram,nlevimgs;
  float             xshift,yshift,psi,zero=0.;
  matrix2D<double>  Mref,Maux;

  // Set-up matrices, etc.
  Mref.resize(images[0]());
  Maux.resize(Mref);
  Mref.set_Xmipp_origin();
  Maux.set_Xmipp_origin();
 
  // Use piramidal combination of images to construct an initial reference
  nlev=SF.ImgNo();
  int i=0;
  while (nlev>1) {
    nlev=nlev/2;
    i++;
  }
  nlev=i;
  n_piram=(int)pow(2.,(double)nlev);
  
  // Copy n_piram to a new temporary array of ImageXmipp
  vector<ImageXmipp>  imgpspc;
  for (imgno=0;imgno<n_piram;imgno++) imgpspc.push_back(images[imgno]);

  cerr << "  Piramidal combination of "<<n_piram<<" images" <<endl;
  init_progress_bar(n_piram);
       
  imgno=0;
  for (int lev=nlev; lev>0; lev--) {
	 
    nlevimgs=(int)pow(2.,(double)lev);
    for (int j=0; j<nlevimgs/2; j++) {
	   
      apply_geom(Mref,imgpspc[2*j].get_transformation_matrix(),imgpspc[2*j](),IS_INV,DONT_WRAP);
      Mref.set_Xmipp_origin();
	   
      
      if (do_acf) {
	align_acf(imgpspc[2*j+1],Mref,zero,zero,Ri,Ro);
      } else if (do_complete) {
	success[imgno]=align_complete_search(images[imgno],Mref,zero,
					     zero,psi_interval,shift_range,Ri,Ro);
      } else {
	/** FIRST **/
	if (do_trans) align_trans(imgpspc[2*j+1],Mref,zero);
	if (do_rot)   align_rot(imgpspc[2*j+1],Mref,zero,Ri,Ro);
	/** SECOND **/
	if (do_trans) align_trans(imgpspc[2*j+1],Mref,zero);
	if (do_rot)   align_rot(imgpspc[2*j+1],Mref,zero,Ri,Ro);
      }

      // Re-calculate average image
      apply_geom(Maux,imgpspc[2*j+1].get_transformation_matrix(),imgpspc[2*j+1](),IS_INV,DONT_WRAP);
      Maux.set_Xmipp_origin();
      Mref=(Mref+Maux)/2;
      imgpspc[j]()=Mref;
      imgpspc[j].Psi()=0.;
      imgpspc[j].Xoff()=0.;
      imgpspc[j].Yoff()=0.;
	   
      imgno++;
      if (imgno%100==0) progress_bar(imgno);
      
    }//loop over images
	 
  }//loop over levels
  progress_bar(n_piram);
  imgpspc.clear();

  // Write out inter-mediate reference
  Iref()=Mref;
  Iref.write(fn_ave);

  Maux.core_deallocate();
  Mref.core_deallocate();
}

// Alignment of all images by iterative refinement  ========================================
void Prog_align2d_prm::refinement() _THROW {

  int               n_refined;
  float             curr_max_shift,curr_max_rot,xshift,yshift,psi,zero=0.;
  matrix2D<double>  Mref,Maux,Msum;

  // Set-up matrices, etc.
  Mref=Iref();
  Mref.set_Xmipp_origin();
  Maux.resize(Mref);
  Msum.resize(Mref);
  Maux.set_Xmipp_origin();
  Msum.set_Xmipp_origin();
  Msum.init_zeros();

  
  n_refined=0;
  cerr << "  Alignment:  iteration " << 1 <<" of "<< Niter << " (with "<<n_images<<" images)" <<endl;
  for (int iter=0; iter < Niter; iter++) {

    if (iter>0) cerr<<"  Refinement: iteration "<<iter+1<<" of "<<Niter<<endl;
    if (iter==(Niter-1)) {
      curr_max_rot=max_rot; 
      curr_max_shift=max_shift;
    } else {
      curr_max_rot=zero; 
      curr_max_shift=zero;
    }

    init_progress_bar(n_images);
    for (int imgno=0; imgno < n_images; imgno++) {

      if (iter!=0) { 
	// Subtract current image from the reference
	apply_geom(Maux,images[imgno].get_transformation_matrix(),images[imgno](),IS_INV,DONT_WRAP);
	Maux.set_Xmipp_origin();
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mref) {
	  dMij(Mref,i,j)*=n_refined;
	  dMij(Mref,i,j)-=dMij(Maux,i,j);
	  dMij(Mref,i,j)/=(n_refined-1);
	}
      }
      
      // Align translationally and rotationally
      if (do_acf) {
	success[imgno]=align_acf(images[imgno],Mref,curr_max_shift,curr_max_rot,Ri,Ro);
      } else if (do_complete) {
	success[imgno]=align_complete_search(images[imgno],Mref,curr_max_shift,
					     curr_max_rot,psi_interval,shift_range,Ri,Ro);
      } else {
	if (do_trans) success[imgno]=align_trans(images[imgno],Mref,curr_max_shift);
	if (do_rot && success[imgno]) success[imgno]=align_rot(images[imgno],Mref,curr_max_rot,Ri,Ro);
      }
      
      if (!success[imgno]) {
	n_refined--;
      } else {
	if (iter==0) {
	  // Add refined images to form a new reference
	  n_refined++;
	  apply_geom(Maux,images[imgno].get_transformation_matrix(),images[imgno](),IS_INV,DONT_WRAP);
	  Maux.set_Xmipp_origin();
	  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Msum) {
	    dMij(Msum,i,j)+=dMij(Maux,i,j);
	  }
	} else {
	  // Add refined image to reference again
	  apply_geom(Maux,images[imgno].get_transformation_matrix(),images[imgno](),IS_INV,DONT_WRAP);
	  Maux.set_Xmipp_origin();
	  FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mref) {
	    dMij(Mref,i,j)*=(n_refined-1);
	    dMij(Mref,i,j)+=dMij(Maux,i,j);
	    dMij(Mref,i,j)/=n_refined;
	  }
	} 
      }
      if (imgno%100==0) progress_bar(imgno);

    } // loop over all images
    progress_bar(n_images);

    if (n_images>n_refined) cerr << "  Discarded "<<n_images-n_refined<<" images." <<endl;
    if (iter==0) FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Msum) {
      dMij(Mref,i,j)=dMij(Msum,i,j)/(n_refined);}

    // Write out inter-mediate reference
    Mref.set_Xmipp_origin();
    Iref()=Mref;
    Iref.write(fn_ave);

  } // loop over iterations
  
  Maux.core_deallocate();
  Mref.core_deallocate();
  Msum.core_deallocate();
}

// Write out results  ========================================================
void Prog_align2d_prm::align2d() _THROW {

  // Read in all images
  FileName fn_img;
  ImageXmipp Itmp;
  n_images=0;
  SF.go_beginning();
  while ((!SF.eof())) {
    fn_img=SF.NextImg();
    if (fn_img == "") continue;  
    Itmp.read(fn_img);
    Itmp().set_Xmipp_origin();
    images.push_back(Itmp);
    success.push_back(true);
    n_images++;
  }

  // Filter if necessary

  if (do_filter) {
    FourierMask fmask;
    fmask.w1=sam/resol; 
    fmask.raised_w=0.1;
    fmask.FilterShape=RAISED_COSINE;
    fmask.FilterBand=LOWPASS;
    matrix2D< complex<double> > fft;
    FourierTransform(images[0](), fft);
    fmask.generate_mask(fft);
    for (int imgno=0;imgno<n_images;imgno++) {
      FourierTransform(images[imgno](), fft);
      fmask.apply_mask_Fourier(fft);
      InverseFourierTransform(fft,images[imgno]());
    }
  }

  // Get Reference (either from file or from piramidal combination of images)
  if (fn_ref!="") {
    Iref.read(fn_ref);
    if (Iref.Psi()!=0. || Iref.Xoff()!=0. || Iref.Yoff()!=0.) 
      cerr << "WARNING!! Transformation in header of reference will NOT be applied"<<endl;
  } else do_pspc();

  // Circular mask around reference image
  matrix2D<int> mask;
  mask.resize(Iref().RowNo(),Iref().ColNo());
  mask.set_Xmipp_origin();
  BinaryCrownMask(mask,0,(int)(Iref().ColNo()/2),INNER_MASK);
  Iref().set_Xmipp_origin();
  apply_binary_mask(mask,Iref(),Iref());

  // Do actual alignment & iteratively refine the reference 
  refinement();
  
  // Write out docfile
  if (fn_doc!="") {
    DocFile           DFo;
    DFo.reserve(n_images);
    DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5)");
    matrix1D<double>  dataline(5);
    for (int imgno=0; imgno < n_images; imgno++) {
      dataline(0)=images[imgno].Phi();
      dataline(1)=images[imgno].Theta();
      dataline(2)=images[imgno].Psi();
      dataline(3)=images[imgno].Xoff();
      dataline(4)=images[imgno].Yoff();
      DFo.append_comment(images[imgno].name());
      DFo.append_data_line(dataline);
    }
    DFo.write(fn_doc);
  }
     
  // Write out images & selfile
  FileName          fn_out;
  SelFile           SFo;
  SFo.reserve(n_images);
  for (int imgno=0; imgno < n_images; imgno++) {
    fn_img=images[imgno].name();
    if (oext!="") fn_out=fn_img.without_extension()+"."+oext;
    else fn_out=fn_img;
    if (do_filter) {
      Itmp.read(fn_img);
      Itmp.Xoff()=images[imgno].Xoff();
      Itmp.Yoff()=images[imgno].Yoff();
      Itmp.Psi()=images[imgno].Psi();
      Itmp.write(fn_out);
    } else images[imgno].write(fn_out);
    if (success[imgno]) SFo.insert(fn_out,SelLine::ACTIVE);
    else SFo.insert(fn_out,SelLine::DISCARDED);
  }
  fn_out=fn_sel;
  if (oext!="") fn_out=fn_out.insert_before_extension("_"+oext);
  SFo.write(fn_out);
  
  images.clear();
}
