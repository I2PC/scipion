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
#include "../Prog_WBP.hh"

// Read arguments ==========================================================
void Prog_WBP_prm::read(int argc, char **argv) _THROW  {

  fn_sel=get_param(argc,argv,"-i");
  SF.read(fn_sel);
  SF.ImgSize(dim,dim);
  apply_shifts=!check_param(argc,argv,"-dont_apply_shifts");
  fn_out =  get_param(argc,argv,"-o","wbp.vol"); 
  fn_sym =  get_param(argc,argv,"-sym",""); 
  threshold=AtoF(get_param(argc,argv,"-threshold","1.667"));
  diameter=2*AtoI(get_param(argc,argv,"-radius","0"));
  if (diameter==0) diameter=dim;
  sampling=AtoF(get_param(argc,argv,"-filsam","4"));
  do_all_matrices=check_param(argc,argv,"-use_each_image");
  // Hidden
  verb=AtoI(get_param(argc,argv,"-verb","1"));

}

// Show ====================================================================
void Prog_WBP_prm::show() {

  if (verb>0) {
    // To screen
    cerr << " ================================================================="<<endl;
    cerr <<" Weighted-back projection (arbitrary geometry) "<<endl;
    cerr << " ================================================================="<<endl;
    cerr << " Input selfile            : "<< fn_sel<<endl;
    cerr << " Output volume            : "<< fn_out<<endl;
    cerr << " Reconstruction radius    : "<< diameter/2<<endl;
    cerr << " Filter threshold         : "<< threshold<<endl;
    if (fn_sym!="")
    cerr << " Symmetry file:           : "<< fn_sym<<endl;
    if (!apply_shifts)
    cerr << " --> Do not apply shifts upon reading the images"<<endl;
    if (do_all_matrices)
    cerr << " --> Use all projection directions in arbitrary geometry filter"<<endl;
    else
    cerr << " --> Use sampled directions for filter, sampling = " <<sampling<<endl;
    cerr << " -----------------------------------------------------------------"<<endl;
  }

}

// Usage ====================================================================
void Prog_WBP_prm::usage() {

  // To screen
  cerr << "  Usage:\n";
  cerr << "  WBP <options>\n";
  cerr << "   -i <input selfile>          : selection file with input images \n";
  cerr << " [ -o <name=\"wbp.vol\">         : filename for output volume \n";
  cerr << " [ -radius <int=dim/2> ]       : Reconstruction radius \n";
  cerr << " [ -sym <symfile> ]            : Enforce symmetry \n";
  cerr << " [ -threshold <float=1.66> ]   : Lower threshold for filter values \n";
  cerr << " [ -filsam <float=4> ]         : Angular sampling rate for geometry filter \n";
  cerr << " [ -use_each_image]            : Use each image instead of sampled representatives for filter \n";
  cerr << " [ -dont_apply_shifts ]        : dont apply origin offsets as stored in the image headers\n";
  cerr << " -----------------------------------------------------------------"<<endl;

}

void Prog_WBP_prm::produce_Side_info() {

  // Fill arrays of transformation matrices
  if (do_all_matrices) get_all_matrices(SF);
  else get_sampled_matrices(SF);

}

void Prog_WBP_prm::get_sampled_matrices(SelFile &SF) {
  
  DocFile           DFlib;
  SymList           SL;
  headerXmipp       head;
  matrix2D<double>  A(3,3);
  matrix2D<double>  L(4,4), R(4,4);
  double            newrot,newtilt,newpsi;
  int               NN,dir,optdir;
  vector<int>       count_imgs;

  if (verb>0) cerr <<"--> Sampling the filter ..."<<endl;

  // Create an even projection direction distribution
  make_even_distribution(DFlib,sampling,SL,false);

  // Only now read symmetry file
  NN=DFlib.LineNo();
  if (fn_sym!="") {
    SL.read_sym_file(fn_sym);
    NN*=(SL.SymsNo()+1);
  } 
  count_imgs.resize(NN+1);

  // Each experimental image contributes to the nearest of these directions
  SF.go_beginning();
  while (!SF.eof()) {
    head.read(SF.NextImg());
    count_imgs[find_nearest_direction(head.Phi(),head.Theta(),DFlib,0,1)]+=1;
    if (fn_sym!="") 
      // Also add symmetry-related projection directions
      for (int i=0; i<SL.SymsNo(); i++) {
	SL.get_matrices(i,L,R);
	L.resize(3,3);R.resize(3,3);
	Euler_apply_transf(L,R,head.Phi(),head.Theta(),head.Psi(),newrot,newtilt,newpsi);
	count_imgs[find_nearest_direction(newrot,newtilt,DFlib,0,1)]+=1;
      }
  }

  no_mats=0;
  for (int i=0; i<NN; i++) if (count_imgs[i]>0.) no_mats++;
  mat_g=(column*)malloc(no_mats*sizeof(column));
  no_mats=0;
  for (int i=1; i<NN+1; i++) {
    if (count_imgs[i]>0.) {
      DFlib.get_angles(i,newrot,newtilt,newpsi,"rot","tilt","psi");
      Euler_angles2matrix (newrot,-newtilt,newpsi,A);
      mat_g[no_mats].zero=A(2,0);
      mat_g[no_mats].one=A(2,1);
      mat_g[no_mats].two=A(2,2);
      mat_g[no_mats].count=count_imgs[i];
      no_mats++;
    }
  }

}

// Fill array with transformation matrices needed for arbitrary geometry filter
void Prog_WBP_prm::get_all_matrices(SelFile &SF) {

  headerXmipp      head;
  matrix2D<double> A(3,3);
  SymList          SL;
  matrix2D<double> L(4,4), R(4,4);
  double           newrot,newtilt,newpsi;
  int              NN;

  SF.go_beginning();
  no_mats=0;

  NN=SF.ImgNo();
  if (fn_sym!="") {
    SL.read_sym_file(fn_sym);
    NN*=(SL.SymsNo()+1);
  } 
  mat_g=(column*)malloc(NN*sizeof(column));

  while (!SF.eof()) {
    head.read(SF.NextImg());
    Euler_angles2matrix (head.Phi(),-head.Theta(),head.Psi(),A);
    mat_g[no_mats].zero=A(2,0);
    mat_g[no_mats].one=A(2,1);
    mat_g[no_mats].two=A(2,2);
    mat_g[no_mats].count=1.;
    no_mats++;
    if (fn_sym!="") {
      // Also add symmetry-related projection directions
      for (int i=0; i<SL.SymsNo(); i++) {
	SL.get_matrices(i,L,R);
	L.resize(3,3);R.resize(3,3);
	Euler_apply_transf(L,R,head.Phi(),-head.Theta(),head.Psi(),newrot,newtilt,newpsi);
	Euler_angles2matrix (newrot,newtilt,newpsi,A);
	mat_g[no_mats].zero=A(2,0);
	mat_g[no_mats].one=A(2,1);
	mat_g[no_mats].two=A(2,2);
	mat_g[no_mats].count=1.;
	no_mats++;
      }
    }
  }

}

// Simple backprojection of a single image
void Prog_WBP_prm::simple_backprojection(Projection &img, VolumeXmipp &vol, 
					 int diameter) {
  int i, j, k, l, m;
  matrix2D<double> A(3,3);
  float dim2, x, y, z, xp, yp;
  float value1, value2, scalex, scaley, scale1, value;
  float radius2, x2, y2, z2, z2_plus_y2;

  // Use minus-tilt, because code copied from OldXmipp
  Euler_angles2matrix(img.rot(),-img.tilt(),img.psi(),A);
  A=A.inv();

  radius2 = diameter/2.;
  radius2 = radius2*radius2;
  dim2 = dim/2;

  for (i = 0; i < dim; i++) {
    z = -i + dim2;   /*** Z points upwards ***/
    z2 = z*z;
    for (j = 0; j < dim; j++) {
      y = j - dim2;
      y2 = y*y;
      z2_plus_y2 = z2 + y2;
      x = 0 - dim2;   /***** X for k == 0 *****/
      xp = x*A(0,0) + y*A(1,0) + z*A(2,0) + dim2;
      yp = x*A(0,1) + y*A(1,1) + z*A(2,1) + dim2;
      for (k = 0; k < dim; k++, xp += A(0,0), yp += A(0,1), x++) {
	x2 = x*x;
	if (x2 + z2_plus_y2 > radius2)
	  continue;
	if ((xp >= (dim-1) || xp < 0) || (yp >= (dim-1) || yp < 0))
	  continue;
	
	/**** interpolation ****/
	l = (int)yp;
	m = (int)xp;
	scalex = xp - m;
	scaley = yp - l;
	scale1 = 1. - scalex;
	value1 = scalex*dMij(img(),l,m+1) + scale1*dMij(img(),l,m);
	value2 = scalex*dMij(img(),l+1,m+1) + scale1*dMij(img(),l+1,m);
	value  = scaley*value2 + (1.-scaley)*value1;
	dVkij(vol(),i,j,k) += value;

      }
    }
  }


}

// Calculate the filter in 2D and apply ======================================
void Prog_WBP_prm::apply_2Dfilter_arbitrary_geometry(SelFile &SF, VolumeXmipp &vol)  {

  int               c,nn,imgno;
  float             factor, argum, weight, x, y;
  Projection        proj;
  FourierImageXmipp IMG;
  matrix2D<double>  A(3,3);

  factor=(float)diameter;
  vol().resize(dim,dim,dim);
  vol().set_Xmipp_origin();
  vol().init_zeros();

  // Tabulated sinc
  tabsinc TSINC(0.001,dim);

   // Initialize time bar
  if (verb>0) cerr <<"--> Back-projecting ..."<<endl;
  nn=SF.ImgNo();
  if (verb>0) init_progress_bar(nn);
  c=MAX(1,nn/60);

  mat_f=(column*)malloc(no_mats*sizeof(column));

  SF.go_beginning();
  imgno=0;
  while (!SF.eof()) {
    proj.read(SF.NextImg(),apply_shifts);
    proj().set_Xmipp_origin();
    Euler_angles2matrix (proj.rot(),-proj.tilt(),proj.psi(),A);
    A=A.inv();
    FourierTransform( proj(), IMG());
    CenterFFT(IMG(),true);

    // loop over all transformation matrices
    for (int k=0; k<no_mats; k++) {
      mat_f[k].zero=A(0,0)*mat_g[k].zero+
	            A(1,0)*mat_g[k].one+
	            A(2,0)*mat_g[k].two; 
      mat_f[k].one=A(0,1)*mat_g[k].zero+
                   A(1,1)*mat_g[k].one+
	           A(2,1)*mat_g[k].two;
    }

    FOR_ALL_ELEMENTS_IN_MATRIX2D(IMG()) {
      y=(float)i; 
      x=(float)j;
      weight = 0.;
      for (int k=0; k<no_mats; k++) {
	argum = diameter/(float)dim*
	  (x*mat_f[k].zero + y*mat_f[k].one);
	// The following line is the most expensive of all...
	weight += mat_g[k].count*TSINC(argum);
      }
      if (fabs(weight) < threshold) {
	MAT_ELEM(IMG(),i,j) /= (threshold*factor);
      } else {
	MAT_ELEM(IMG(),i,j) /= (weight*factor);
      }

    }

    // Calculate back-projection with the filtered projection
    CenterFFT(IMG(),false);
    InverseFourierTransform(IMG(), proj());
    simple_backprojection(proj, vol, diameter);

    if (verb>0) if (imgno%c==0) progress_bar(imgno);
    imgno++;

  }
  if (verb>0) progress_bar(nn);

  // Symmetrize if necessary
  if (fn_sym!="") {
    SymList SL;
    VolumeXmipp aux;
    aux().resize(vol());
    SL.read_sym_file(fn_sym);
    symmetrize(SL, vol, aux);
    vol()=aux();
    aux.clear();
  }

  // free memory
  free(mat_g);
  free(mat_f);

}

