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

#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippGeometry.hh>
#include <XmippInterface/xmippAPH.hh>
#include "../Prog_Spots2RealSpace2D.hh"
#include "../Prog_art_crystal.hh"

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
      fnaph_ref=get_param(fh_param,"APH Ref file",0,"",
         3007,"Spot2RealSpace2D_Parameters::read: APH Ref File not found");
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
              "Spot2RealSpace2D_Parameters::read: Cannot read Phase Shift");}
//    NOT longer needed
//      tilt_sign=AtoI(get_param(fh_param,"tilt sign",0,NULL,
//         3007,"Spot2RealSpace2D_Parameters::read: Tilt Sign not found"));
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
void Spot2RealSpace2D_Parameters::produce_SideInfo() _THROW {
   // Read APH files
   aph_file.read(fnaph_in);
   if (fnaph_ref!="") aph_ref.read(fnaph_ref);
    ImageXmipp save;
    save()=aph_file.spots_abs; save.write("Spots_ABS0.xmp");
    save()=aph_file.spots_arg; save.write("Spots_ARG0.xmp");

   // Compute Reference lattice vectors in both spaces
   aph_ref=(fnaph_ref!="")?aph_ref:aph_file;   
   // Reference
   matrix2D<double> Vref(2,2), Vfile(2,2),Proj_Matrix(2,2);
   Vref.setCol(0,aph_ref.astar);
   Vref.setCol(1,aph_ref.bstar);
   vol_latt_vec=aph_ref.Xdim*(Vref.transpose()).inv();
   
   //FACTOR
   {
    double gammarefstar=atan2(YY(aph_ref.bstar),XX(aph_ref.bstar)) -
                        atan2(YY(aph_ref.astar),XX(aph_ref.astar));
    double gammastar   =atan2(YY(aph_file.bstar),XX(aph_file.bstar)) -
                        atan2(YY(aph_file.astar),XX(aph_file.astar));
    double A =aph_ref.astar.module();
    double B =aph_ref.bstar.module();
    double AT=aph_file.astar.module();
    double BT=aph_file.bstar.module();
    double C1 = (A*A)/(AT*AT);
    double C2 = (B*B)/(BT*BT);
    double C3 = (A*B)/(AT*BT*cos(gammastar));
    double C4 = (A*B*cos(gammarefstar))/(AT*BT*cos(gammastar));
    double AX = C2*(C1/C3)*(C1/C3) - C1;
    double BX = C2 + 2.0*C2*(C1/C3)*((C1-C4)/C3) - C1;
    double CX = C2*((C1-C4)/C3)*((C1-C4)/C3);
    double DISC = BX*BX - 4.0*AX*CX;
    double PSQ1 = (-BX + sqrt(DISC))/(2.0*AX);
    double PSQ2 = (-BX - sqrt(DISC))/(2.0*AX);
    double PP=0;
    if(PSQ1 < 0.0 && PSQ2 < 0.0) {
                    cout << "HORROR: DISCRIMINANT LESS THAN ZERO" <<endl;
                    exit(1);
                    }  
    if(PSQ1 > 0.0) PP = sqrt(PSQ1);
    if(PSQ2 > 0.0) PP = sqrt(PSQ2);
    double QQ = (C1/C3)*PP + ((C1-C4)/C3)*(1.0/PP);
    SamplingScale = sqrt(C1*(PP*PP + 1.0));
    }
    
   // actual projection
   Vfile.setCol(0,aph_file.astar);
   Vfile.setCol(1,aph_file.bstar);

   aph_file.astar *= SamplingScale;
   aph_file.bstar *= SamplingScale;
   Vfile.setCol(0,aph_file.astar);
   Vfile.setCol(1,aph_file.bstar);
   proj_latt_vec=aph_ref.Xdim*(Vfile.transpose()).inv();
   Proj_Matrix=proj_latt_vec * vol_latt_vec.inv();
   //-----------------------------------------------
   //Tilt
   //-----------------------------------------------
   char all_done='n';

   if(Proj_Matrix.det()> 1.0)
       {
       cerr << "\nHORROR: Matrix determinant greater than 1,\n"
               "It is equal to:" << Proj_Matrix.det() <<
               "\nIf the value is almost 1 and tilt should be\n"
	       " arround 0, then it may be OK (I will make tilt=0)\n"
	       "\nPRESS anykey AND enter to continue\n";
       cin >> all_done; 
       tilt = 0.;
       }
   else    
       tilt = /* tilt_sign * */ acos(Proj_Matrix.det());
   //sign is given by the user 

   //-----------------------------------------------
   // ROT  and PSI
   //-----------------------------------------------
   rot=psi=0.;
   double v1=Proj_Matrix(0,0); double v2=Proj_Matrix(0,1);
   double v3=Proj_Matrix(1,0); double v4=Proj_Matrix(1,1);
   all_done='n';
   if(v2< 0.01 && v3 <0.01)// special case
        {
         cerr << "\nInput Files: " <<fnaph_in << " " << fnaph_ref;
	 
	 if(v4<v1)//shortening along y direction
	     {
             cerr << "\nWARMING, posible rounding error. Since"
                 "\nit seems  to be a shortening along the y axis"
                 "\nshould I make rot=270 and psi=90?\n"
                 "\n y= set, n= continue (y/n)";
             do {
                 cin >> all_done;     
	         rot=PI*3./2.;psi=PI/2.;
                } while ( !((all_done !='y') || 
                            (all_done !='n' ))   );
	     
	     }
	 else if(v4< -0.99 && v1< -0.99)//shortening along y direction
	     {
             cerr << "\nWARMING, posible rounding error. "
                 "\nshould I make rot=180 and psi=0? (if tilt is small this"
                 "\nis not very important)?\n\n"
                 "\n y= set, n= continue (y/n)";
             do {
                 cin >> all_done;     
	         rot=PI;psi=0;
                } while ( !((all_done !='y') || 
                            (all_done !='n' ))   );
	     
	     }
	  else
	    {   
             cerr << "\nWARMING, posible rounding error. Should I set"
                     "\nrot  and psi to =0 (if tilt is small this"
                     "\nis not very important)?\n"
                     "\n y= set to zero, n= continue (y/n)";
        	do {
                    cin >> all_done;     
                    rot=psi=0.;
        	   } while ( !((all_done !='y') || 
                               (all_done !='n' ))   );
	      }	       
        }           
   if(fabs(tilt)<0.0873 && all_done == 'n') //tilt zero or small 
       { rot = atan2(-v3,v4); psi=0.; 
         cerr << "WARMING, rot calculated using tilt=0\n"; 
        }
   else if(all_done == 'n')
       { 

         rot = atan2(v2*v2+v4*v4-1.,v1*v2+v3*v4);    
         psi = atan2( -(v1*v3) -(v2*v4), v1*v1+v2*v2-1.);

       //Our matrix v1, v2, v3, v4 is only a guess of the REAL euler Matrix
       //if any of the values in the atan2 functions is close to zero 
       //the sign of the angle may be wrong. A way to overcome this
       //is to recalculate the crystal vectors from the projection vectors and
       //if there is a sign error to invert the signs.
       
         matrix2D<double> E2D;
         Euler_angles2matrix(RAD2DEG(rot),
                             RAD2DEG(tilt),
                             RAD2DEG(psi),E2D); E2D.resize(2,2);
         if(v1-E2D(0,0)>0.1 || v2-E2D(0,1)>0.1 || 
            v3-E2D(1,0)>0.1 || v4-E2D(1,1)>0.1) psi += M_PI;
         Euler_angles2matrix(RAD2DEG(rot),
                             RAD2DEG(tilt),
                             RAD2DEG(psi),E2D); E2D.resize(2,2);
         if(v1-E2D(0,0)>0.1 || v2-E2D(0,1)>0.1 || 
            v3-E2D(1,0)>0.1 || v4-E2D(1,1)>0.1) 
         cerr << "\nOh My! spots2realspace can't find the angles\n";    
                 
        }
   //Ask user what angles want to use.
   EULER_CLIPPING_RAD(rot,tilt,psi);
//   if(all_done == 'n')//ask always
   {
double mrc_taxa, mrc_tilt;        
       cout << "\n\nThey are two solutions to the system"
            << "\nSelect the solution you think is the correct"
            <<  "\n(it will be stored in the image header )"
            << "\n\t (1) Rotational tilt and psi angle: " << 
                   RAD2DEG(rot) << ", " << 
                   RAD2DEG(tilt)<< ", " << 
                   RAD2DEG(psi) 
            << "\n\t                                    " << 
                   RAD2DEG(rot) << ", " << 
                   RAD2DEG(tilt)+180.<< ", " << 
                 -(RAD2DEG(psi)+180.) 
            << "\n\t                                    " << 
                   RAD2DEG(rot)+180. << ", " << 
                  -RAD2DEG(tilt)<< ", " << 
                   RAD2DEG(psi)-180. 
            << "\n\t                                    " << 
                   RAD2DEG(rot)+180. << ", " << 
                  -RAD2DEG(tilt)+180.<< ", " << 
                  -(RAD2DEG(psi)-180+180.);
       Euler_to_MRC( rot,  tilt,  psi, &mrc_tilt, &mrc_taxa);
       cout << "\n\tMRC (taxa,tilt): " << RAD2DEG(mrc_taxa) << " " 
                                     << RAD2DEG(mrc_tilt) << endl;
       cout << "\n\t (2) Rotational tilt and psi angle: "<< 
                   RAD2DEG(rot) << ", " << 
                   RAD2DEG(-tilt)<< ", " << 
                   RAD2DEG(psi) 
            << "\n\t                                    " << 
                   RAD2DEG(rot) << ", " << 
                   RAD2DEG(-tilt)+180.<< ", " << 
                 -(RAD2DEG(psi)+180.) 
            << "\n\t                                    " << 
                   RAD2DEG(rot)+180. << ", " << 
                  -RAD2DEG(-tilt)<< ", " << 
                   RAD2DEG(psi)-180. 
            << "\n\t                                    " << 
                   RAD2DEG(rot)+180. << ", " << 
                  -RAD2DEG(-tilt)+180.<< ", " << 
                  -(RAD2DEG(psi)-180+180.);
       Euler_to_MRC( rot, -tilt,  psi, &mrc_tilt, &mrc_taxa);
       cout << "\n\tMRC (taxa,tilt): " << RAD2DEG(mrc_taxa) << " " 
                                       << RAD2DEG(mrc_tilt) << endl;
		   
       cout<<  "\n(type 1, 2 <enter>)";
       int answer;
       // only two cases are different but it make live easier 
       // for the user to have the four
       do {
          cin >> answer;     
	  
          if (answer ==2 ) {tilt = -tilt; break;}
          else if (answer == 1) break;
          else cout << "\nType 1, 2!";
          } while (TRUE );
	  
     }
          
    //I like tilt between -90 and 90
    if(tilt > M_PI) tilt -= 2*M_PI;   

    // Translate symmetry group
    if      (str_symmetry_group=="P1")       symmetry_group = sym_P1;
    else if (str_symmetry_group=="P2")       symmetry_group = sym_P2;
    else if (str_symmetry_group=="P21")      symmetry_group = sym_P2_1;
    else if (str_symmetry_group=="C2")       symmetry_group = sym_C2;
    else if (str_symmetry_group=="P222")     symmetry_group = sym_P222;
    else if (str_symmetry_group=="P2221")    symmetry_group = sym_P222_1;
    else if (str_symmetry_group=="P22121")   symmetry_group = sym_P22_12_1;
    else if (str_symmetry_group=="P4")       symmetry_group = sym_P4;
    else if (str_symmetry_group=="P422")     symmetry_group = sym_P422;
    else if (str_symmetry_group=="P4212")    symmetry_group = sym_P42_12;
    else if (str_symmetry_group=="P3")       symmetry_group = sym_P3;
    else if (str_symmetry_group=="P312")     symmetry_group = sym_P312;
    else if (str_symmetry_group=="P6")       symmetry_group = sym_P6;
    else if (str_symmetry_group=="P622")     symmetry_group = sym_P622;
    else
       REPORT_ERROR(1,(string)"Spot2RealSpace2D_Parameters::produce_SideInfo:"
          " Unknown symmetry group:"+str_symmetry_group);
    
    // Generate symmetrical reflections
    if (generate_symmetrical_reflections)
       aph_file.generate_symmetrical_reflections(symmetry_group);
   
    save()=aph_file.spots_abs; save.write("Spots_ABS1.xmp");
    save()=aph_file.spots_arg; save.write("Spots_ARG1.xmp");
}

/* Show parameters --------------------------------------------------------- */
ostream& operator << (ostream &o, const Spot2RealSpace2D_Parameters &prm) {
   o << "APH Input          : " << prm.fnaph_in << endl
     << "APH Reference      : " << prm.fnaph_ref << endl
     << "Output Spider File : " << prm.fn_out << endl
     << "Maximum IQ         : " << prm.maxIQ << endl
     << "Rotational angle   : " << RAD2DEG(prm.rot) << endl
     << "Tilt angle         : " << RAD2DEG(prm.tilt) << endl
     << "psi angle          : " << RAD2DEG(prm.psi) << endl
     << "Lattice a          : " << prm.proj_latt_vec.Col(0).transpose() << endl
     << "Lattice b          : " << prm.proj_latt_vec.Col(1).transpose() << endl
     << "Ouput Layout (x,y) : " << prm.NoCells_X
                       << " x " << prm.NoCells_Y << endl
     << "Phase_Shift (x,y)  : " << XX(prm.Phase_Shift) << ","
                                << YY(prm.Phase_Shift)  << endl
     << "Keep Contrast      : " << prm.KeepContrast << endl
     << "Cell Size (x,y)    : " << prm.CellXdim << "," << prm.CellYdim << endl
     << "Crystaledge        : " << prm.aph_file.Xdim << endl
     << "Crystaledge(ref)   : " << prm.aph_ref.Xdim << endl
     << "Scale_Factor       : " << prm.Scale_Factor << endl
     << "Sampling Scale*    : " << prm.SamplingScale << endl;
   return o;
}

/* DFT^-1 ------------------------------------------------------------------ */
void IDFT(const matrix2D< complex<double> > &FT, matrix2D<double> &I,
   int ydim, int xdim) {
   I.init_zeros(ydim,xdim);
   time_config();
   init_progress_bar(YSIZE(I));
   for (int k=0; k<YSIZE(I); k++) {
      double wk=twoPI*k/YSIZE(I);
      progress_bar(k);
      for (int l=0; l<XSIZE(I); l++) {
         double wl=twoPI*l/XSIZE(I);
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
}

/* Main routine for transforming ------------------------------------------- */
void ROUT_Spots2RealSpace(Spot2RealSpace2D_Parameters &prm, 
   Projection &prj) {
   prm.produce_SideInfo();
   cout << prm;
   
   // Apply phase shift
   double phase_shift=(prm.KeepContrast)? 0:180;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(prm.aph_file.spots_arg) {
       MAT_ELEM(prm.aph_file.spots_arg,i,j)  =  (double) DEG2RAD 
           (MAT_ELEM(prm.aph_file.spots_arg,i,j) + phase_shift + 
           i*(YY(prm.Phase_Shift)/* +prm.mirror_phase_Y*/) +
           j*(XX(prm.Phase_Shift)/* +prm.mirror_phase_X*/));
   }

   // Resize and symmetrize Fourier Transform
   int kmin=STARTINGY (prm.aph_file.spots_abs);
   int kmax=FINISHINGY(prm.aph_file.spots_abs);
   int hmin=STARTINGX (prm.aph_file.spots_abs);
   int hmax=FINISHINGX(prm.aph_file.spots_abs);

   int ksize=MAX(ABS(kmin),ABS(kmax));
   int hsize=MAX(ABS(hmin),ABS(hmax));
   matrix2D< complex<double> > FT;
   FT.init_zeros(2*ksize+1,2*hsize+1);
   STARTINGY(FT)=-ksize;
   STARTINGX(FT)=-hsize;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(prm.aph_file.spots_abs) {
      if (ABS(MAT_ELEM(prm.aph_file.IQ,i,j))<=prm.maxIQ &&
          MAT_ELEM(prm.aph_file.IQ,i,j)!=0) {
	 MAT_ELEM(FT, i, j) = polar(MAT_ELEM(prm.aph_file.spots_abs,i,j),
                                    MAT_ELEM(prm.aph_file.spots_arg,i,j) );
 	 MAT_ELEM(FT,-i,-j) = conj(MAT_ELEM(FT, i, j));
      }
   }

   // Create projection
   if (prm.CellXdim==-1 && prm.CellYdim==-1) {
      prm.CellYdim=YSIZE(FT);
      prm.CellXdim=XSIZE(FT);
   }
   IDFT(FT,prj(),prm.CellYdim,prm.CellXdim);
   prj() *= prm.Scale_Factor/sqrt((double)prm.aph_file.Xdim*prm.aph_file.Ydim);
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
   
#ifdef DEBUGROUT_Spots2RealSpace
    cout << "\no_taxa " << prm.taxa << "  o_tilt " << prm.tilt << endl;
    cout << "taxa "     << new_rot  << "  tilt "   << new_tilt << endl;
#endif
}
#undef DEBUGROUT_Spots2RealSpace

