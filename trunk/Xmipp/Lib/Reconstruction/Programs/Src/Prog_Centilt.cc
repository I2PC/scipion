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
#include "../Prog_Centilt.hh"

// Read arguments ==========================================================
void Prog_centilt_prm::read(int argc, char **argv) _THROW  {

  FileName fn_sel;

  // Selfile with untilted images
  fn_sel=get_param(argc,argv,"-u");
  SFu.read(fn_sel);
  // Selfile with tilted images
  fn_sel=get_param(argc,argv,"-t");
  SFt.read(fn_sel);
  if (SFu.ImgNo()!=SFt.ImgNo()) 
    REPORT_ERROR(1,"Unequal number of active images in untilted and tilted selfiles");
  // Extension if not to overwrite input images
  oext=get_param(argc,argv,"-oext","");
  // Write out document file?
  fn_doc=get_param(argc,argv,"-doc","");
  // Maximum shift (discard images that shift more in last iteration)
  max_shift=AtoF(get_param(argc,argv,"-max_shift","0"));
  // Force x-shift to be zero?
  force_x_zero=check_param(argc,argv,"-force_x_zero");
  // Perform centering?
  do_center=!check_param(argc,argv,"-skip_centering");
  // Perform cosine stretching?
  do_stretch=!check_param(argc,argv,"-skip_stretching");
}

// Show ====================================================================
void Prog_centilt_prm::show() {
    cerr << " Selfile untilted images              : "<<  SFu.name()<<endl;
    cerr << " Selfile tilted images                : "<<  SFt.name()<<endl;
  if (oext!="") 
    cerr << " Output extension for tilted images   : "<< oext<<endl;
  if (max_shift!=0) 
    cerr << " Discard images that shift more than  : "<<max_shift<<endl;
  if (fn_doc!="") 
    cerr << " Output document file (tilted images) : "<<fn_doc<<endl;
  if (force_x_zero)
    cerr << " Force x-shift to be zero "<<endl;
  if (!do_stretch) 
    cerr << " Skip cosine stretching "<<endl;
  if (!do_center) 
    cerr << " Skip centering "<<endl;

}

// usage ===================================================================
void Prog_centilt_prm::usage() {
  cerr << "Usage:  "<<endl;
  cerr << "  centilt [options]"<<endl;
  cerr << "   -u <selfile>             : Selfile containing untilted images \n"
       << "   -t <selfile>             : Selfile containing tilted images \n"
       << " [ -oext <extension> ]      : For output tilted images; if not to overwrite input\n"
       << " [ -max_shift <float> ]     : Discard images that shift more [pix]\n"
       << " [ -doc <docfile> ]         : write output document file with rotations & translations \n"
       << " [ -force_x_zero ]          : Force x-shift to be zero \n"
       << " [ -skip_stretching ]       : do not perform cosine stretching \n"
       << " [ -skip_centering ]        : do not perform centering, i.e. only modify angles \n"
       << endl;
}

// Center one tilted image  =====================================================
bool Prog_centilt_prm::center_tilted_image(const ImageXmipp &Iu, ImageXmipp &It, double &ccf) _THROW {

  matrix2D<double> A(3,3),Maux(It()),Mcorr(It());
  int              n_max=-1;
  bool             neighbourhood=TRUE;
  int              imax,jmax,i_actual,j_actual,x_zero=1;
  double           maxcorr,xmax,ymax,sumcorr;
  float            xshift,yshift,shift;

  if (do_stretch) {
    // Cosine stretching, store stretched image in Maux
    A.init_identity();
    A(0,0)=COSD(It.Theta());
    Maux.init_zeros();
    apply_geom(Maux,A,It(),IS_INV,DONT_WRAP);
  } else Maux=It();

  // Calculate cross-correlation
  correlation_matrix(Maux,Iu(),Mcorr);
  Mcorr.statistics_adjust(0.,1.);
  if (force_x_zero) {
    x_zero=0;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr) {
      if (j!=0) MAT_ELEM(Mcorr,i,j)=0.;
    }
  }
  Mcorr.max_index(imax,jmax);
  maxcorr=MAT_ELEM(Mcorr,imax,jmax);

  while (neighbourhood) {  
    n_max ++;
    for (int i=-n_max; i <= n_max; i++)
      for (int j=-n_max*x_zero; j <= n_max*x_zero; j++) {   
	i_actual = i+imax;
	j_actual = j+jmax;
	if (i_actual < Mcorr.startingY()  || j_actual < Mcorr.startingX() && 
	    i_actual > Mcorr.finishingY() || j_actual > Mcorr.finishingX() ) 
	  neighbourhood=FALSE;
	else if (maxcorr/1.414 > MAT_ELEM(Mcorr,i_actual,j_actual))
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

  // Calculate correlation coefficient 
  A.init_identity();
  A(0,2)=-xshift;
  A(1,2)=-yshift;
  apply_geom(Maux,A,Maux,IS_INV,DONT_WRAP);
  Maux.set_Xmipp_origin();
  ccf=correlation_index(Iu(),Maux);

  if (do_stretch) xshift*=COSD(It.Theta());
  shift=sqrt(xshift*xshift+yshift*yshift);
  if ((max_shift<XMIPP_EQUAL_ACCURACY) || (shift < max_shift)) {
    // Store shift in the header of the image
    It.Xoff()+=xshift;
    It.Yoff()+=yshift;
    It.Phi()=Iu.Psi();
    return TRUE;
  } else return FALSE;

  Maux.core_deallocate();
  Mcorr.core_deallocate();

}

// Main program  ===============================================================
void Prog_centilt_prm::centilt() _THROW {

  DocFile           DFo;
  FileName          fn_img;
  ImageXmipp        Iu, It;
  matrix2D<double>  Maux, A(3,3);
  matrix1D<double>  dataline(6);
  double            ccf,outside;
  bool              OK;
  int               imgno,barf,n_images,n_discarded;

  if (fn_doc!="") {
    DFo.reserve(SFt.ImgNo()*2+1);
    DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Corr (6)");
  }

  n_images=SFt.ImgNo();
  n_discarded=0;
  cerr << "  Centering of "<<n_images<<" tilted images"<<endl;
  init_progress_bar(n_images);
  barf=MAX(1,(int)(1+(n_images/60)));
  imgno=0;
  while (imgno<n_images) {

    SFu.go_beginning();
    SFt.go_beginning();
    SFu.jump(imgno);
    SFt.jump(imgno-n_discarded);
    // Read in untilted image and apply shifts (center) and Phi (align tilt-axis with y-axis) 
    Iu.read(SFu.get_current_file());
    Iu().set_Xmipp_origin();
    Euler_angles2matrix (Iu.Phi(),0.,0.,A);
    A(0,2)=-Iu.Xoff();
    A(1,2)=-Iu.Yoff();
    outside=dMij(Iu(),0,0);
    Iu().self_apply_geom(A,IS_INV,DONT_WRAP,outside);
    // Read in tilted image and apply Psi (align tilt-axis with y-axis) and shifts if present
    It.read(SFt.get_current_file());
    // Store original matrix for later output
    Maux.resize(It());
    Maux=It();
    Euler_angles2matrix (It.Psi(),0.,0.,A);
    A(0,2)=-It.Xoff();
    A(1,2)=-It.Yoff();
    outside=dMij(It(),0,0);
    It().self_apply_geom(A,IS_INV,DONT_WRAP,outside);
    It().set_Xmipp_origin();

    if (do_center) OK=center_tilted_image(Iu,It,ccf);
    else {
      OK=true;
      ccf=1.;
      It.Phi()=Iu.Psi();
    }
    if (OK) {
      fn_img=It.name();
      if (oext!="") {
	fn_img=fn_img.without_extension()+"."+oext;
      }
      // Add line to document file
      if (fn_doc!="") {
	dataline(0)=It.Phi();
	dataline(1)=It.Theta();
	dataline(2)=It.Psi();
	dataline(3)=It.Xoff();
	dataline(4)=It.Yoff();
	dataline(5)=ccf;
	DFo.append_comment(fn_img);
	DFo.append_data_line(dataline);
      }
      SFt.set_current_filename(fn_img);
      // Re-store original matrix & write out tilted image
      It()=Maux;
      It.write(fn_img);
    } else {
      SFt.remove(It.name());
      n_discarded++;
    }

    imgno++;
    if (imgno%barf==0) progress_bar(imgno);
  }
  progress_bar(n_images);
  if (max_shift>0) cerr << "  Discarded "<<n_discarded<<" tilted images that shifted too much"<<endl;

  // Write out document file
  if (fn_doc!="") DFo.write(fn_doc);

  // Write out selfile
  fn_img=SFt.name();
  if (oext!="") fn_img=fn_img.insert_before_extension("_"+oext);
  SFt.write(fn_img);


}

