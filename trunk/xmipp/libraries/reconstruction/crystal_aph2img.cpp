/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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

#include <data/args.h>
#include <data/geometry.h>
#include <interface/aph.h>

#include "crystal_aph2img.h"
#include "art_crystal.h"

#define GCC_VERSION (__GNUC__ * 10000 \
   + __GNUC_MINOR__ * 100 \
   + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
   #include <sstream>
#else
   #include <strstream.h>
#endif

#define twoPI 2*M_PI

void Spot2RealSpace2D_Parameters::read_from_file(const FileName &fnprm)
{
   FILE *fh_param;
   if ((fh_param = fopen(fnprm.c_str(), "r")) == NULL)
      REPORT_ERROR(3005,
         (string)"Spot2RealSpace2D_Parameters::read: There is a problem "
         "opening the file "+fnprm);

   try {
#ifdef REMOVE
      fnaph_ref=get_param(fh_param,"APH Ref file",0,"",
         3007,"Spot2RealSpace2D_Parameters::: APH Ref File not found");
#endif
      fnaph_in=get_param(fh_param,"APH file",0,NULL,
         3007,"Spot2RealSpace2D_Parameters::read: APH File not found");
      fn_out=get_param(fh_param,"Output file",0,NULL,
         3007,"Spot2RealSpace2D_Parameters::read: Output File not found");
      maxIQ=AtoI(get_param(fh_param,"maxIQ",0,"6"));
      string aux_str=get_param(fh_param,"Output cells no",0,"1 1");
      #if GCC_VERSION < 30300
	 istrstream *is=NULL;
         is = new istrstream(aux_str.c_str());
      #else
	 istringstream *is=NULL;
         is = new istringstream(aux_str.c_str());
      #endif
      try {*is >> NoCells_X >> NoCells_Y;}
      catch (...) {REPORT_ERROR(3007,
              "Spot2RealSpace2D_Parameters::read: Cannot read Output Cell no");}
      a_mag=AtoF(get_param(fh_param,"a_mag",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: a_mag not found"));
      b_mag=AtoF(get_param(fh_param,"b_mag",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: b_mag not found"));
      a_b_ang=DEG2RAD(AtoF(get_param(fh_param,"a_b_ang",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: a_b_ang not found")));
      taxa=DEG2RAD(AtoF(get_param(fh_param,"taxa",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: Taxa not found")));
      mrc_tilt=DEG2RAD(AtoF(get_param(fh_param,"tilt",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: tilt not found")));
      SamplingRate=AtoF(get_param(fh_param,"sampling rate",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: sampling rate not found"));
      Phase_Shift.resize(2);
      aux_str=get_param(fh_param,"Phase Shift",0,"180.0 180.0");
      delete is;
      #if GCC_VERSION < 30300
         is = new istrstream(aux_str.c_str());
      #else
         is = new istringstream(aux_str.c_str());
      #endif
      try {*is >> XX(Phase_Shift) >> YY(Phase_Shift);}
      catch (...) {REPORT_ERROR(3007,
              "Spot2RealSpace2D_Parameters::read: Cannot read Phase Shift");}
      mrc_label=AtoI(get_param(fh_param,"mrc_label",0,"-1"));
      KeepContrast=AtoI(get_param(fh_param,"Keep Contrast",0,"1"));
      aux_str=get_param(fh_param,"Cell Size");
      delete is;
      #if GCC_VERSION < 30300
         is = new istrstream(aux_str.c_str());
      #else
         is = new istringstream(aux_str.c_str());
      #endif
      try {*is >> CellXdim >> CellYdim;}
      catch (...) {REPORT_ERROR(3007,
               "Spot2RealSpace2D_Parameters::read: Cannot read Cell Size");}
      delete is;
      Scale_Factor=AtoF(get_param(fh_param,"Scale Factor",0,"-1"));
      generate_symmetrical_reflections=
         check_param(fh_param,"Generate symmetrical reflections");
      str_symmetry_group=
         get_param(fh_param,"Symmetry group",0,"P1");
//      align_a_axis_with_x=check_param(fh_param,"Align A axis");
   } catch (Xmipp_error XE) {
      cout << XE << endl;
      REPORT_ERROR(3007,(string)"There is an error reading "+fnprm);
   }
   fclose(fh_param);
}

/* Produce Side info ------------------------------------------------------- */
void Spot2RealSpace2D_Parameters::produce_SideInfo() {

    // Translate symmetry group
    if      (str_symmetry_group=="P1")      {symmetry_group = sym_P1;
                                             a_b_ang = 0;}//not sure about this
    else if (str_symmetry_group=="P2")      symmetry_group = sym_P2;
    else if (str_symmetry_group=="P21")     symmetry_group = sym_P2_1;
    else if (str_symmetry_group=="C2")      symmetry_group = sym_C2;
    else if (str_symmetry_group=="P222")    symmetry_group = sym_P222;
    else if (str_symmetry_group=="P2_122")  symmetry_group = sym_P2_122;
    else if (str_symmetry_group=="P22_12")  symmetry_group = sym_P22_12;
    else if (str_symmetry_group=="P22121")  symmetry_group = sym_P22_12_1;
    else if (str_symmetry_group=="P4")      symmetry_group = sym_P4;
    else if (str_symmetry_group=="P422")    symmetry_group = sym_P422;
    else if (str_symmetry_group=="P4212")   symmetry_group = sym_P42_12;
    else if (str_symmetry_group=="P3")      symmetry_group = sym_P3;
    else if (str_symmetry_group=="P312")    symmetry_group = sym_P312;
    else if (str_symmetry_group=="P6")      symmetry_group = sym_P6;
    else if (str_symmetry_group=="P622")    symmetry_group = sym_P622;
    else
       REPORT_ERROR(1,(string)"Spot2RealSpace2D_Parameters::produce_SideInfo:"
          " Unknown symmetry group:"+str_symmetry_group);

    // convert taxa, tilt to Euler angles
    // psi will be set to 0
    MRC_to_Euler(taxa, mrc_tilt, &rot, &tilt, &psi);
#ifdef REMOVE
    // Undo moving spots to asymmetryc unit
    aph_file.unasymmetrization( a_mag,  b_mag, taxa,  mrc_tilt,
			       symmetry_group);
#endif
#ifdef DEBUG
    ImageXmipp save;
    save()=aph_file.spots_abs; save.write("Spots_ABS0.xmp");
    save()=aph_file.spots_arg; save.write("Spots_ARG0.xmp");
#endif

}

/* Show parameters --------------------------------------------------------- */
ostream& operator << (ostream &o, const Spot2RealSpace2D_Parameters &prm) {
   o << "APH Input          : " << prm.fnaph_in << endl
#ifdef REMOVE
     << "APH Reference      : " << prm.fnaph_ref << endl
#endif
     << "Output Spider File : " << prm.fn_out << endl
     << "Maximum IQ         : " << prm.maxIQ << endl
     << "a_mag(A)           : " << prm.a_mag << endl
     << "b_mag(A)           : " << prm.b_mag << endl
     << "a->b angle         : " << RAD2DEG(prm.a_b_ang) << endl
     << "Rotational angle   : " << RAD2DEG(prm.rot) << endl
     << "Tilt angle         : " << RAD2DEG(prm.tilt) << endl
     << "psi angle          : " << RAD2DEG(prm.psi) << endl
     << "taxa angle         : " << RAD2DEG(prm.taxa) << endl
     << "mrc_tilt angle     : " << RAD2DEG(prm.mrc_tilt) << endl
     << "mrc_label          : " << prm.mrc_label << endl
     << "sampling rate      : " << prm.SamplingRate << endl
#ifdef REMOVE
     << "Lattice a          : " << prm.proj_latt_vec.Col(0).transpose() << endl
     << "Lattice b          : " << prm.proj_latt_vec.Col(1).transpose() << endl
#endif
     << "Ouput Layout (x,y) : " << prm.NoCells_X
                       << " x " << prm.NoCells_Y << endl
     << "Phase_Shift (x,y)  : " << XX(prm.Phase_Shift) << ","
                                << YY(prm.Phase_Shift)  << endl
     << "Keep Contrast      : " << prm.KeepContrast << endl
     << "Cell Size (x,y)    : " << prm.CellXdim << "," << prm.CellYdim << endl
#ifdef REMOVE
     << "Crystaledge        : " << prm.aph_file.Xdim << endl
     << "Crystaledge(ref)   : " << prm.aph_ref.Xdim << endl
#endif
     << "Scale_Factor       : " << prm.Scale_Factor << endl
#ifdef REMOVE
     << "Sampling Scale*    : " << prm.SamplingScale << endl
#endif
     << "3D Symmetry group  : " << prm.symmetry_group << endl;
   return o;
}

/* DFT^-1 ------------------------------------------------------------------ */
void IDFT(const matrix2D< complex<double> > &FT, matrix2D<double> &I,
   int ydim, int xdim) {
   I.init_zeros(ydim,xdim);
   time_config();
   init_progress_bar(YSIZE(I));
   for (int k=0; k<YSIZE(I); k++) {
      double wk= -1.*twoPI*k/YSIZE(I);
      progress_bar(k);
      for (int l=0; l<XSIZE(I); l++) {
         double wl= -1.*twoPI*l/XSIZE(I);
         for (int K=STARTINGY(FT); K<=FINISHINGY(FT); K++) {
            double wK=wk*K;
            for (int L=STARTINGX(FT); L<=FINISHINGX(FT); L++) {
               if (MAT_ELEM(FT,K,L)==0.0) continue;
	       MAT_ELEM(I,k,l) +=
	          real(MAT_ELEM(FT,K,L))*cos(wK+wl*L)
	         +imag(MAT_ELEM(FT,K,L))*sin(wK+wl*L);
            }
         }
      }
   }
   progress_bar(YSIZE(I));
   I /= (double)(xdim*ydim);
}

/* Main routine for transforming ------------------------------------------- */
void ROUT_Spots2RealSpace(Spot2RealSpace2D_Parameters &prm,
   Projection &prj) {
   int h,k;
   //symmery should be reforced after APPLYYING SHIFT
   prm.produce_SideInfo();
   cout << prm;
   // Read APH files, vector size must be in Amstrongs
   (prm.aph_file).read(prm.fnaph_in,prm.mrc_label);
   // Now we should be able to check if the points does belong to the
   // plane (only for aph with z*) or not
   // if data is P1, all point should be OK if not several should
   // not belong to the plane.

   // Apply phase shift
   // Up to here angles are in degrees
   // When copied to FT convert degrees to radians
   double phase_shift=(prm.KeepContrast)? 0.:180.;
   if(phase_shift!=0 || (prm.Phase_Shift).module()!=0)
      for(int i=0;i< ((prm.aph_file).aph_data_vector).size();i++) {
	  h = (prm.aph_file.aph_data_vector[i]).h;
	  k = (prm.aph_file.aph_data_vector[i]).k;
	  ((prm.aph_file).aph_data_vector[i]).phase  =
              ((prm.aph_file).aph_data_vector[i]).phase + phase_shift +
              k*(YY(prm.Phase_Shift)) +
              h*(XX(prm.Phase_Shift));
      }
   // Move all the reflections to the same plane
   #define DEBUG
   #ifdef DEBUG
   (prm.aph_file).write("IN1");
   #endif
   #undef DEBUG



   if (prm.generate_symmetrical_reflections)
      ((prm.aph_file).unasymmetrization)( prm.a_mag,    prm.b_mag,
                                          prm.taxa,     prm.mrc_tilt,
			                  prm.a_b_ang,  prm.symmetry_group,
					  prm.Counter);
   #define DEBUG
   #ifdef DEBUG
   (prm.aph_file).write("OUT1");
   #endif
   #undef DEBUG
   //alloc a 2D matrix and transfer data
   // Resize and symmetrize Fourier Transform
   #define DEBUG
   #ifdef DEGUG
   cout << "prm.a_mag: "          << prm.a_mag             << endl;
   cout << "prm.b_mag: "          << prm.b_mag             << endl;
   cout << "prm.taxa: "           << RAD2DEG(prm.taxa)     << endl;
   cout << "prm.mrc_tilt: "       << RAD2DEG(prm.mrc_tilt) << endl;
   cout << "prm.symmetry_group: " << prm.symmetry_group    << endl;
   #endif
   #undef DEBUG
   int ksize=MAX(ABS((prm.aph_file).min_k),ABS((prm.aph_file).max_k));
   int hsize=MAX(ABS((prm.aph_file).min_h),ABS((prm.aph_file).max_h));

   //When data is copied to FT move it from degrees to radians.
   //It is done inside the polar function
   matrix2D< complex<double> > FT;
   FT.init_zeros(2*ksize+1,2*hsize+1);
   STARTINGY(FT)=-ksize;
   STARTINGX(FT)=-hsize;

   for(int i=0;i< prm.aph_file.aph_data_vector.size();i++) {
      if (ABS((prm.aph_file.aph_data_vector[i]).IQ)<=prm.maxIQ &&
              (prm.aph_file.aph_data_vector[i]).amp!=0) {
	 h = (prm.aph_file.aph_data_vector[i]).h;
	 k = (prm.aph_file.aph_data_vector[i]).k;
	
	 MAT_ELEM(FT, k, h) = polar((prm.aph_file.aph_data_vector[i]).amp,
                            DEG2RAD((prm.aph_file.aph_data_vector[i]).phase));
 	 MAT_ELEM(FT,-k,-h) = conj(MAT_ELEM(FT, k, h));
      }
   }

   // Create projection
   if (prm.CellXdim==-1 && prm.CellYdim==-1) {
      prm.CellYdim=YSIZE(FT);
      prm.CellXdim=XSIZE(FT);
   }
   IDFT(FT,prj(),prm.CellYdim,prm.CellXdim);
   prj() *= prm.Scale_Factor;
//ADD HERE ANY OTHER FACTOR, MAY BE AS OPTION

   // Any repetition?
   prj().resize(prm.CellYdim*prm.NoCells_Y, prm.CellXdim*prm.NoCells_X);
   for (int I=0; I<prm.NoCells_Y; I++)
      for (int J=0; J<prm.NoCells_X; J++) {
         if (I==0 && J==0) continue;
         for (int i=0; i<prm.CellYdim; i++)
	    for (int j=0; j<prm.CellXdim; j++)
	      IMGPIXEL(prj,i+I*prm.CellYdim,j+J*prm.CellXdim)=IMGPIXEL(prj,i,j);
	 }

   // Set Euler angles and Save Image
   prj.set_angles(RAD2DEG(prm.rot), RAD2DEG(prm.tilt), RAD2DEG(prm.psi));
   prj.write(prm.fn_out);
}

/* Main routine for transforming ------------------------------------------- */
void ROUT_RealSpace2Spots(RealSpace2Spots2D_Parameters &prm,
   Projection &prj) {
   //show parameters
   cout << prm;
   //loop through the sel file
   matrix2D< complex<double> > FT;
   // CALL DFT
   prj.read(prm.fn_in);
//   #define DEBUG
   #ifdef DEBUG
   cout << prj();
   #endif
   #undef DEBUG
     DFT(prj(),FT);
     // SAve aph format
     (prm.aph_file).clear();
     //prm.mrc_labelinitial value in parameter file
     (prm.aph_file).read_mrc_label=prm.mrc_label;
     spot tmp_spot;
     int my_counter_spot=1;
     int my_counter_film=1;
	
     matrix1D<double> normal_plane(3);
     matrix1D<double> aux_vector(3);
     int DoesIntersect;
     matrix1D<double> intersection_point;
     matrix2D<double>  A(3,3);
//     double mrc_tilt, mrc_taxa;
     matrix2D<double> E,E2D, vp(2,2), Vp(2,2), v(2,2);
     //compute a* and b* in projection coordinate system
     //we will call them ap* and bp*
     Euler_angles2matrix (prm.taxa, prm.mrc_tilt, 0., E);
     E2D=E;
     E2D.resize(2,2);
     v.setCol(0,prm.vector_a);
     v.setCol(1,prm.vector_b);
     vp=E2D*v;
     Vp=(vp.transpose()).inv();//ap*,bp* column vector

//     normal_plane=E2D.Col(2);

     FOR_ALL_ELEMENTS_IN_MATRIX2D(FT) {
//     if(j==0 && i==1) {cout << "EXTRA" << abs(MAT_ELEM(FT,i,j));}
//     else continue;
	 tmp_spot.amp	     = abs(MAT_ELEM(FT,i,j));
         if(tmp_spot.amp<0.0001) continue;
         tmp_spot.h	     = j; //x
         tmp_spot.k	     = i; //y
         if(tmp_spot.h==0 && tmp_spot.k< 0) continue;
	 tmp_spot.phase      = RAD2DEG(arg(MAT_ELEM(FT,i,j)));
	 tmp_spot.FILM       = (prm.aph_file).read_mrc_label;
	 tmp_spot.IQ	     = 1;
	 tmp_spot.FLMWGT     = 1;
	 tmp_spot.BACK	     = 1;
 	 tmp_spot.myCTF	     = 1;

	 // compute z*
	 // from euler angles get normal vector
	 // compute intersección central plane with drid line
//tenemos que calcular los vectores de red en espacio real y quitar  el 128	
//	 DoesIntersect=line_plane_intersection(normal_plane,
//                           vector_R3(0.,0.,1.),
//                           intersection_point,
//                           vector_R3((double)tmp_spot.h/128.,
//			   (double)tmp_spot.k/128.,0.),
//			   0.
//			   );
//	 if(DoesIntersect!=0)
//            {
//	     cerr << "ERROR: plane and line does not intersec" << endl;
//	     exit(1);
//	     }
//         tmp_spot.zstar      = ZZ(intersection_point);
         aux_vector = vector_R3((double)tmp_spot.h,
	                        (double)tmp_spot.k,0.);//h
	 Vp.resize(3,3);
	 tmp_spot.zstar=ZZ((E.inv())*Vp*aux_vector);			
	 prm.aph_file.aph_data_vector.push_back(tmp_spot);
         #undef DEBUG
//	 #define DEBUG
	 #ifdef DEBUG
         fprintf(stdout," %d  %d %f %f %f %d %d %f %f %f\n",
			       tmp_spot.h,	
			       tmp_spot.k,	
			       tmp_spot.zstar,
			       tmp_spot.amp,	
			       tmp_spot.phase,
			       tmp_spot.FILM,
			       tmp_spot.IQ,	
			       tmp_spot.FLMWGT,
			       tmp_spot.BACK,
			       tmp_spot.myCTF);
         fflush(stdout);
         #endif
     }//rellenar las estructura de aph
      // sort values
      sort(prm.aph_file.aph_data_vector.begin(),
           prm.aph_file.aph_data_vector.end());
     //salvar el fichero aph
     (prm.aph_file).write(prm.fnaph_out);

}

void RealSpace2Spots2D_Parameters::read_from_file(const FileName &fnprm)
{
   FILE *fh_param;
   if ((fh_param = fopen(fnprm.c_str(), "r")) == NULL)
      REPORT_ERROR(3005,
         (string)"RealSpace2Spots2D_Parameters::read: There is a problem "
         "opening the file "+fnprm);

   try {
      fn_in=get_param(fh_param,"INput file",0,NULL,
         3007,"RealSpace2Spots2D_Parameters::read: Input File not found");
      fnaph_out=get_param(fh_param,"APH file",0,NULL,
         3007,"RealSpace2Spots2D_Parameters::read: Output File not found");
      {
      string aux_str=get_param(fh_param,"lattice_a",0,NULL,
         3007,"RealSpace2Spots2D_Parameters::read: crystal vector a not found");
      #if GCC_VERSION < 30300
	 istrstream *is=NULL;
         is = new istrstream(aux_str.c_str());
      #else
	 istringstream *is=NULL;
         is = new istringstream(aux_str.c_str());
      #endif
      try {*is >> vector_a(0) >> vector_a(1);}
      catch (...) {REPORT_ERROR(3007,
              "RealSpace2Spots2D_Parameters::read: crystal vector a");}
      delete is;
      }
      {
      string aux_str=get_param(fh_param,"lattice_b",0,NULL,
         3007,"RealSpace2Spots2D_Parameters::read: crystal vector b not found");
      #if GCC_VERSION < 30300
	 istrstream *is=NULL;
         is = new istrstream(aux_str.c_str());
      #else
	 istringstream *is=NULL;
         is = new istringstream(aux_str.c_str());
      #endif
      try {*is >> vector_b(0) >> vector_b(1);}
      catch (...) {REPORT_ERROR(3007,
              "RealSpace2Spots2D_Parameters::read: crystal vector b Angles");}
      delete is;
      }
      mrc_label=AtoI(get_param(fh_param,"mrc_label",0,"-1"));
      taxa=(AtoF(get_param(fh_param,"taxa",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: Taxa not found")));
      mrc_tilt=(AtoF(get_param(fh_param,"tilt",0,NULL,
	 3007,"Spot2RealSpace2D_Parameters::read: tilt not found")));
   } catch (Xmipp_error XE) {
      cout << XE << endl;
      REPORT_ERROR(3007,(string)"There is an error reading "+fnprm);
   }
   fclose(fh_param);
}
     #define SHOW(str,v,i) x=v(0,i); y=v(1,i); cout << str << "("<< x \
              << "," << y << ") (" << sqrt(x*x+y*y) << ")\n";

/* Show parameters --------------------------------------------------------- */
ostream& operator << (ostream &o, const RealSpace2Spots2D_Parameters &prm) {
   o << "Input XMIPP File   : " << prm.fn_in << endl
     << "Output APH File    : " << prm.fnaph_out << endl
     << "taxa               : " << (prm.taxa) << endl
     << "Tilt angle         : " << (prm.mrc_tilt) << endl
     << "mrc_label          : " << prm.mrc_label << endl
     << "lattice_a          : " << prm.vector_a << endl
     << "lattice_b          : " << prm.vector_b << endl
     ;
   return o;
}
/* DFT ------------------------------------------------------------------ */
void DFT(const matrix2D<double> &I, matrix2D< complex<double> > &FT)
 {
//   matrix2D<double> I; I.resize(3,3);
//   I(0,0)=1;I(0,1)=2;I(0,2)=3;
//   I(1,0)=4;I(1,1)=5;I(1,2)=6;
//   I(2,0)=7;I(2,1)=8;I(2,2)=9;
   int xdim=I.ColNo();
   int ydim=I.RowNo();
   double myreal,myimag;
   FT.init_zeros(ydim,xdim);
   FT.set_Xmipp_origin();

   time_config();
   init_progress_bar(YSIZE(I));
   for (int k=0; k<YSIZE(I); k++) {
      double wk=twoPI*k/YSIZE(I);
      progress_bar(k);
      for (int l=0; l<XSIZE(I); l++) {
         double wl=twoPI*l/XSIZE(I);
         for (int K=STARTINGY(FT); K<=FINISHINGY(FT); K++) {
            double wK=wk*K;
            for (int L=0/*STARTINGX(FT)*/; L<=FINISHINGX(FT); L++) {
               if (MAT_ELEM(I,k,l)==0.0) continue;
	       myreal = MAT_ELEM(I,k,l)*cos(wK+wl*L);
	       myimag = MAT_ELEM(I,k,l)*sin(wK+wl*L);
  	       MAT_ELEM(FT,K,L) += complex<double>(myreal,0);
	       MAT_ELEM(FT,K,L) -= complex<double>(0,myimag);
//cout << k << " " << l << "  " << K << " " << L
//     << " " << MAT_ELEM(I,k,l) << " "
//     << myreal << " " << myimag << endl;
//cout << MAT_ELEM(FT,K,L) << endl;
             }
         }
      }
   }
//cout << I;
//cout << FT;
   progress_bar(YSIZE(I));

 }
