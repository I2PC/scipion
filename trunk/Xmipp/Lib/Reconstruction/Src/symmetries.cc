/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include <stdio.h>
#include "../symmetries.hh"

// Read Symmetry file ======================================================
void SymList::read_sym_file(FileName fn_sym, double accuracy) {
   int i,j,k,l;
   FILE *fpoii;
   char line[80];
   char *auxstr;
   double ang_incr, rot_ang;
   int  fold;
   matrix2D<double> L(4,4), R(4,4);
   matrix1D<double> axis(3), shift(3);
   
   // Open file ---------------------------------------------------------
   if ((fpoii = fopen (fn_sym.c_str(), "r")) == NULL)
      REPORT_ERROR(3005,(string)"SymList::read_sym_file:Can't open file: "
         +fn_sym);
   //reset space_group
   space_group=0;
   // Count the number of symmetries ------------------------------------
   true_symNo=0;
   // count number of axis and mirror planes. It will help to identify
   // the crystallographic symmetry
   
   int no_axis, no_mirror_planes;
   no_axis= no_mirror_planes= 0;
   
   while (fgets (line, 79, fpoii) != NULL) {
      if (line[0]==';' || line[0]=='#' || line[0]=='\0') continue;
      auxstr = first_token(line);
      if (auxstr==NULL) {
         cout << line;
         cout << "Wrong line in symmetry file, the line is skipped\n";
         continue;
      }
      if (strcmp(auxstr,"rot_axis")==0) {
         auxstr=next_token(); fold=AtoI(auxstr); 
	 true_symNo += (fold-1); no_axis++;
      } else if (strcmp(auxstr,"mirror_plane")==0) {true_symNo++;
                                                    no_mirror_planes++;}
      else if (strcmp(auxstr,"P4212")==0) true_symNo+=7;
   }

   fseek (fpoii, 0L, SEEK_SET);

   // Ask for memory
   __L.resize(4*true_symNo,4);
   __R.resize(4*true_symNo,4);
   __shift.resize(true_symNo,3);
   __chain_length.resize(true_symNo);
   __chain_length.init_constant(1);

   // Read symmetry parameters
   i=0;
   shift.init_zeros();
   while (fgets (line, 79, fpoii) != NULL) {
      if (line[0]==';' || line[0]=='#' || line[0]=='\0') continue;
      auxstr = first_token(line);
      // Rotational axis ---------------------------------------------------
      if (strcmp(auxstr,"rot_axis")==0) {
         auxstr=next_token(); fold=AtoI(auxstr);
         auxstr=next_token(); axis.X()=AtoF(auxstr);
         auxstr=next_token(); axis.Y()=AtoF(auxstr);
         auxstr=next_token(); axis.Z()=AtoF(auxstr);
         ang_incr=360./fold;
         L.init_identity();
         for (j=1, rot_ang=ang_incr; j<fold; j++, rot_ang+=ang_incr) {
            R=rot3D_matrix(rot_ang,axis);
	    set_shift(i,shift);
            set_matrices(i++,L,R.transpose());
         }
         __sym_elements++;
      // Mirror plane ------------------------------------------------------
      } else if (strcmp(auxstr,"mirror_plane")==0) {
         auxstr=next_token(); axis.X()=AtoF(auxstr);
         auxstr=next_token(); axis.Y()=AtoF(auxstr);
         auxstr=next_token(); axis.Z()=AtoF(auxstr);
         L.init_identity(); L(2,2)=-1;
         matrix2D<double> A=align_with_Z(axis); A=A.transpose();
         R=A*L*A.inv();
	 set_shift(i,shift);
         set_matrices(i++,L,R);
         __sym_elements++;
      // P4212 -------------------------------------------------------------
      } else if (strcmp(auxstr,"P4212")==0) {
         space_group=sym_P42_12;
         accuracy=-1; // Do not compute subgroup
	 L.init_identity();

      	 // With 0 shift
	 R.init_zeros(); R(3,3)=1;
	 R(0,0)=R(1,1)=-1; R(2,2)=1; set_shift(i,shift); set_matrices(i++,L,R);
	 R.init_zeros(); R(3,3)=1;
	 R(2,2)=-1; R(0,1)=R(1,0)=1; set_shift(i,shift); set_matrices(i++,L,R);
	 R.init_zeros(); R(3,3)=1;
	 R(2,2)=R(0,1)=R(1,0)=-1;    set_shift(i,shift); set_matrices(i++,L,R);

      	 // With 1/2 shift
	 VECTOR_R3(shift,0.5,0.5,0);
	 R.init_zeros(); R(3,3)=1;
	 R(0,1)=-1; R(1,0)=R(2,2)=1; set_shift(i,shift); set_matrices(i++,L,R);
	 R.init_zeros(); R(3,3)=1;
	 R(1,0)=-1; R(0,1)=R(2,2)=1; set_shift(i,shift); set_matrices(i++,L,R);
	 R.init_zeros(); R(3,3)=1;
	 R(0,0)=R(2,2)=-1; R(1,1)=1; set_shift(i,shift); set_matrices(i++,L,R);
	 R.init_zeros(); R(3,3)=1;
	 R(1,1)=R(2,2)=-1; R(0,0)=1; set_shift(i,shift); set_matrices(i++,L,R);

         __sym_elements++;
      }
   }
   fclose(fpoii);
   if (accuracy>0) compute_subgroup(accuracy);

   //possible crystallographic symmetry
   if(no_axis==0 && no_mirror_planes== 0 && 
      true_symNo==7 && space_group==sym_P42_12)
         space_group=sym_P42_12;
   // P4 and P6	 
   else if (no_axis==1 && no_mirror_planes== 0 && 
            fabs(R(2,2)-1.) < XMIPP_EQUAL_ACCURACY &&
	    fabs(R(0,0) - R(1,1)) < XMIPP_EQUAL_ACCURACY &&
	    fabs(R(0,1) + R(1,0)) < XMIPP_EQUAL_ACCURACY )
	    {
	    switch(true_symNo){
	       case(5):
                  space_group=sym_P6;
		  break;
	       case(3):
                  space_group=sym_P4;
		  break;
              default: 
                 space_group=sym_undefined;
		 break;
	       }//switch end
	    }//end else if (no_axis==1 && no_mirror_planes== 0 
   else if (no_axis==0 && no_mirror_planes== 0) 
         space_group=sym_P1;
   else 
        space_group=sym_undefined;	 
   
}

// Get matrix ==============================================================
void SymList::get_matrices(int i, matrix2D<double> &L, matrix2D<double> &R)
   const {
   int k, l;
   L.init_zeros(4,4);
   R.init_zeros(4,4);
   for (k=4*i; k<4*i+4; k++)
      for (l=0; l<4; l++) {
         DIRECT_MAT_ELEM(L,k-4*i,l)=DIRECT_MAT_ELEM(__L,k,l);
         DIRECT_MAT_ELEM(R,k-4*i,l)=DIRECT_MAT_ELEM(__R,k,l);
      }
}

// Set matrix ==============================================================
void SymList::set_matrices(int i, const matrix2D<double> &L,
   const matrix2D<double> &R) {
   int k, l;
   for (k=4*i; k<4*i+4; k++)
      for (l=0; l<4; l++) {
         DIRECT_MAT_ELEM(__L,k,l)=DIRECT_MAT_ELEM(L,k-4*i,l);
         DIRECT_MAT_ELEM(__R,k,l)=DIRECT_MAT_ELEM(R,k-4*i,l);
      }
}

// Get/Set shift ===========================================================
void SymList::get_shift(int i, matrix1D<double> &shift) const {
   shift.resize(3);
   XX(shift)=DIRECT_MAT_ELEM(__shift,i,0);
   YY(shift)=DIRECT_MAT_ELEM(__shift,i,1);
   ZZ(shift)=DIRECT_MAT_ELEM(__shift,i,2);
}

void SymList::set_shift(int i, const matrix1D<double> &shift) _THROW {
   if (XSIZE(shift)!=3)
      REPORT_ERROR(1002,"SymList::add_shift: Shift vector is not 3x1");
   DIRECT_MAT_ELEM(__shift,i,0)=XX(shift);
   DIRECT_MAT_ELEM(__shift,i,1)=YY(shift);
   DIRECT_MAT_ELEM(__shift,i,2)=ZZ(shift);
}

void SymList::add_shift(const matrix1D<double> &shift) _THROW {
   if (XSIZE(shift)!=3)
      REPORT_ERROR(1002,"SymList::add_shift: Shift vector is not 3x1");
   int i=YSIZE(__shift);
   __shift.resize(i+1,3);
   set_shift(i,shift);
}

// Add matrix ==============================================================
void SymList::add_matrices(const matrix2D<double> &L, const matrix2D<double> &R,
   int chain_length) _THROW {
   if (XSIZE(L)!=4 || YSIZE(L)!=4 || XSIZE(R)!=4 || YSIZE(R)!=4 )
      REPORT_ERROR(1002,"SymList::add_matrix: Transformation matrix is not 4x4");
   if (TrueSymsNo()==SymsNo()) {
      __L.resize(__L.ydim+4,4);
      __R.resize(__R.ydim+4,4);
      __chain_length.resize(__chain_length.xdim+1);
   }

   set_matrices(true_symNo,L,R);
   __chain_length(XSIZE(__chain_length)-1)=chain_length;
   true_symNo++;
}

// Compute subgroup ========================================================
bool found_not_tried(const matrix2D<int> &tried, int &i, int &j,
   int true_symNo) {
   i=j=0;
   int n=0;
   while (n!=YSIZE(tried)) {
      if (MAT_ELEM(tried,i,j)==0 && !(i>=true_symNo && j>=true_symNo))
         return TRUE;
      if (i!=n) {
         // Move downwards
         i++;
      } else {
         // Move leftwards
         j--; if (j==-1) {n++; j=n; i=0;}
      }
   }
   return FALSE;
}

//#define DEBUG
void SymList::compute_subgroup(double accuracy) {
   matrix2D<double> I(4,4); I.init_identity();
   matrix2D<double> L1(4,4), R1(4,4), L2(4,4), R2(4,4), newL(4,4), newR(4,4);
   matrix2D<int>    tried(true_symNo, true_symNo);
   matrix1D<double> shift(3); shift.init_zeros();
   int i,j;
   int new_chain_length;
   while (found_not_tried(tried,i,j,true_symNo)) {
      tried(i,j)=1;

      // Form new symmetry matrices
      // if (__chain_length(i)+__chain_length(j)>__sym_elements+2) continue;
      
      get_matrices(i,L1,R1);
      get_matrices(j,L2,R2);
      newL=L1*L2;
      newR=R1*R2;
      new_chain_length=__chain_length(i)+__chain_length(j);

      if (newL.IsIdent() && newR.IsIdent()) continue;

      // Try to find it in current ones
      bool found;
      found=FALSE;
      for (int l=0; l<SymsNo(); l++) {
         get_matrices(l,L1,R1);
         if (newL.equal(L1,accuracy) && newR.equal(R1,accuracy))
            {found=TRUE; break;}
      }
      
      if (!found) {
         #ifdef DEBUG
            cout << "Matrix size " << XSIZE(tried) << " "
                 << "trying " << i << " " << j << " "
                 << "chain length=" << new_chain_length << endl;
            cout << "Result R\n" << newR;
         #endif
         add_matrices(newL,newR,new_chain_length);
	 add_shift(shift);
         tried.resize(YSIZE(tried)+1,XSIZE(tried)+1);
      }
   }
}
   /** Guess Crystallographic space group.  
       Return the  \URL[space group]{
       http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen} number. So
       far it has only been implemented for P1 (1), P4 (75), P4212 (90) and P6
       (168) */

int  SymList::crystallographic_space_group (double mag_a,double mag_b,
				     double ang_a2b_deg) const
{
    
    switch(space_group){
       case sym_undefined: 
       case sym_P1: return(space_group);
       case sym_P4:
	 if ( fabs((mag_a - mag_b)) >XMIPP_EQUAL_ACCURACY ||
	      fabs(ang_a2b_deg -90) >XMIPP_EQUAL_ACCURACY)
         cerr << "\nWARNING: P42 but mag_a != mag_b\n"
	      << " or ang_a2b !=90" << endl;
	 return(space_group);
	 break;
       case sym_P42_12:
	 if ( fabs((mag_a - mag_b)) >XMIPP_EQUAL_ACCURACY ||
	      fabs(ang_a2b_deg -90) >XMIPP_EQUAL_ACCURACY)
         cerr << "\nWARNING: P42_12 but mag_a != mag_b\n"
	      << " or ang_a2b !=90" << endl;
	 return(space_group);
	 break;
       case sym_P6:
	 if ( fabs((mag_a - mag_b)) >XMIPP_EQUAL_ACCURACY ||
	      fabs(ang_a2b_deg -120.) >XMIPP_EQUAL_ACCURACY)
           {
	   cerr << "\nWARNING: marked as P6 but mag_a != mag_b\n"
	        << "or ang_a2b !=120" << endl;
           cerr << "\nWARNING: P1 is assumed\n";
	   return(sym_P1);
	   }
	 else  return(space_group); 
	 break;
       default:
         cerr << "\n Congratulations: you have found a bug in the\n"
	      << "routine crystallographic_space_group or\n"
	      << "You have called to this rotuine BEFORE reading\n"
	      << "the symmetry info" << endl;
	 exit(0);
	 break;
       }//switch(space_group)  end 

}//crystallographic_space_group end

// Symmetrize_crystal_vectors==========================================
//A note: Please realize that we are not repeating code here.
//The class SymList deals with symmetries when expressed in 
//Cartesian space, that is the basis is orthonormal. Here
//we describe symmetries in the crystallographic way
//that is, the basis and the crystal vectors are the same.
//For same symmetries both representations are almost the same
//but in general they are rather different.

//IMPORTANT: matrix orden should match the one used in "read_sym_file"
//if not the wrong angles are assigned to the different matrices

void symmetrize_crystal_vectors(matrix1D<double> &aint, 
			      matrix1D<double> &bint,
			      matrix1D<double> &shift,
			      int space_group,
			      int sym_no,
			      const matrix1D<double> &eprm_aint,
			      const matrix1D<double> &eprm_bint){
//Notice that we should use R.inv and not R to relate eprm.aint and aint
   shift.init_zeros();//I think this init is OK even the vector dim=0
   switch(space_group)
     {
     case(sym_undefined):
     case(sym_P1):
               XX(aint)=   XX(eprm_aint);
               YY(aint)=                   YY(eprm_aint);
               XX(bint)=   XX(eprm_bint);
               YY(bint)=                   YY(eprm_bint);
        break;
     case(sym_P2):       cerr << "\n Group P2 not implemented\n"; exit(1);
        break;
     case(sym_P2_1):     cerr << "\n Group P2_1 not implemented\n"; exit(1);
        break;
     case(sym_C2):       cerr << "\n Group C2 not implemented\n"; exit(1);
        break;
     case(sym_P222):     cerr << "\n Group P222 not implemented\n"; exit(1);
        break;
     case(sym_P222_1):   cerr << "\n Group P222_1 not implemented\n"; exit(1);
        break;
     case(sym_P22_12_1): cerr << "\n Group P22_12_1 not implemented\n"; exit(1);
        break;
     case(sym_P4):       
           switch(sym_no)
	   {
	    case(-1): XX(aint)=   XX(eprm_aint);
                      YY(aint)=                   YY(eprm_aint);
                      XX(bint)=   XX(eprm_bint);
                      YY(bint)=                   YY(eprm_bint);
              break;
	    case(0):  XX(aint)=                 - YY(eprm_aint);  
                      YY(aint)=   XX(eprm_aint);
		      XX(bint)=                 - YY(eprm_bint);  
                      YY(bint)=   XX(eprm_bint);
              break;
	    case(1):  XX(aint)= - XX(eprm_aint);  
                      YY(aint)=                 - YY(eprm_aint);
		      XX(bint)= - XX(eprm_bint);   
                      YY(bint)=                 - YY(eprm_bint);
              break;
	    case(2):  XX(aint)=                   YY(eprm_aint);
                      YY(aint)= - XX(eprm_aint);
                      XX(bint)=                   YY(eprm_bint);
                      YY(bint)= - XX(eprm_bint);
              break;
	   }//switch P4 end   
        break;
     case(sym_P422):     cerr << "\n Group P422 not implemented\n"; exit(1);
        break;
     case(sym_P42_12):
           switch(sym_no)
	   {
            case(-1): XX(aint)=   XX(eprm_aint);
                      YY(aint)=                   YY(eprm_aint);
                      XX(bint)=   XX(eprm_bint);
                      YY(bint)=                   YY(eprm_bint);
              break;
	    case(0):  XX(aint)= - XX(eprm_aint);  
                      YY(aint)=                 - YY(eprm_aint);
		      XX(bint)= - XX(eprm_bint);   
                      YY(bint)=                 - YY(eprm_bint);
              break;
	    case(1):  XX(aint)=                 + YY(eprm_aint);
                      YY(aint)= + XX(eprm_aint);
		      XX(bint)=                 + YY(eprm_bint);  
                      YY(bint)= + XX(eprm_bint);
              break;
	    case(2):  XX(aint)=                 - YY(eprm_aint);  
                      YY(aint)= - XX(eprm_aint);
		      XX(bint)=                 - YY(eprm_bint);   
                      YY(bint)= - XX(eprm_bint);
              break;
	    case(3):  XX(aint)=                 + YY(eprm_aint);  
                      YY(aint)= - XX(eprm_aint);
		      XX(bint)=                 + YY(eprm_bint);  
                      YY(bint)= - XX(eprm_bint);
		    VECTOR_R3(shift, 0.5, 0.5,0);
            break;
	    case(4):  XX(aint)=		        - YY(eprm_aint);
		      YY(aint)= + XX(eprm_aint);
		      XX(bint)= 		- YY(eprm_bint);
		      YY(bint)= + XX(eprm_bint);
		    VECTOR_R3(shift,0.5,0.5,0);
	      break;
	    case(5):  XX(aint)= - XX(eprm_aint);  
		      YY(aint)=		        + YY(eprm_aint);
		      XX(bint)= - XX(eprm_bint);   
		      YY(bint)=		        + YY(eprm_bint);
		    VECTOR_R3(shift,0.5,0.5,0);
	    break;
	    case(6):  XX(aint)= + XX(eprm_aint);  
		      YY(aint)=		        - YY(eprm_aint);
		      XX(bint)= + XX(eprm_bint);   
		      YY(bint)=		        - YY(eprm_bint);
		  VECTOR_R3(shift,0.5,0.5,0);
	    break;
	    default:
                 cout<<"\n Wrong symmetry number "
		       "in symmetrize_crystal_vectors, bye"<<endl; exit(1);
      break;

    
	   }//switch P6 end   
        break;
     case(sym_P3):       cerr << "\n Group P3 not implemented\n"; exit(1);
        break;
     case(sym_P312):     cerr << "\n Group P312 not implemented\n"; exit(1);
        break;
     case(sym_P6): 
           switch(sym_no)
	   {
	    case(-1): XX(aint)=   XX(eprm_aint);
                      YY(aint)=                   YY(eprm_aint);
                      XX(bint)=   XX(eprm_bint);
                      YY(bint)=                   YY(eprm_bint);
              break;
	    case(0):  XX(aint)=   XX(eprm_aint) - YY(eprm_aint);  
                      YY(aint)=   XX(eprm_aint);
		      XX(bint)=   XX(eprm_bint) - YY(eprm_bint);  
                      YY(bint)=   XX(eprm_bint);
              break;
	    case(1):  XX(aint)=                 - YY(eprm_aint);  
                      YY(aint)=   XX(eprm_aint) - YY(eprm_aint);
		      XX(bint)=                 - YY(eprm_bint);   
                      YY(bint)=   XX(eprm_bint) - YY(eprm_bint);
              break;
	    case(2):  XX(aint)= - XX(eprm_aint);
                      YY(aint)=                 - YY(eprm_aint);
                      XX(bint)= - XX(eprm_bint);
                      YY(bint)=                 - YY(eprm_bint);
              break;
	    case(3):  XX(aint)= - XX(eprm_aint) + YY(eprm_aint);  
                      YY(aint)= - XX(eprm_aint);
		      XX(bint)= - XX(eprm_bint) + YY(eprm_bint);  
                      YY(bint)= - XX(eprm_bint);
              break;
	    case(4):  XX(aint)=                 + YY(eprm_aint);  
                      YY(aint)= - XX(eprm_aint) + YY(eprm_aint);
		      XX(bint)=                 + YY(eprm_bint);   
                      YY(bint)= - XX(eprm_bint) + YY(eprm_bint);
              break;
	   }//switch P6 end   
        break;

     case(sym_P622):     cerr << "\n Group P622 not implemented\n"; exit(1);
        break;
     }
			      
}//symmetrize_crystal_vectors end			         

#define Symmetrize_Vol(X) {\
             for (int i=0; i<vol_in.VolumesNo(); i++)\
                  X(vol_in(i),vol_in.grid(i),eprm_aint,eprm_bint,mask,i, \
		  grid_type);\
		  }

// Symmetrize_crystal_volume==========================================
//IMPORTANT: matrix orden should match the one used in "read_sym_file"
//if not the wrong angles are assigned to the different matrices
void symmetrize_crystal_volume(GridVolume &vol_in,
                              const matrix1D<double> &eprm_aint, 
			      const matrix1D<double> &eprm_bint,
			      int eprm_space_group,
			      const matrix2D<int> &mask, int grid_type){

   switch(eprm_space_group)
     {
     case(sym_undefined):
     case(sym_P1):
        break;
     case(sym_P2):       cerr << "\n Group P2 not implemented\n"; exit(1);
        break;
     case(sym_P2_1):     cerr << "\n Group P2_1 not implemented\n"; exit(1);
        break;
     case(sym_C2):       cerr << "\n Group C2 not implemented\n"; exit(1);
        break;
     case(sym_P222):     cerr << "\n Group P222 not implemented\n"; exit(1);
        break;
     case(sym_P222_1):   cerr << "\n Group P222_1 not implemented\n"; exit(1);
        break;
     case(sym_P22_12_1): cerr << "\n Group P22_12_1 not implemented\n"; exit(1);
        break;
     case(sym_P4):       
	Symmetrize_Vol(symmetry_P4)//already has ;
        break;
     case(sym_P422):     cerr << "\n Group P422 not implemented\n"; exit(1);
        break;
     case(sym_P42_12):
	Symmetrize_Vol(symmetry_P42_12)//already has ;
        break;
     case(sym_P3):       cerr << "\n Group P3 not implemented\n"; exit(1);
        break;
     case(sym_P312):     cerr << "\n Group P312 not implemented\n"; exit(1);
        break;
     case(sym_P6): 
	Symmetrize_Vol(symmetry_P6)//already has ;
        break;
     case(sym_P622):     cerr << "\n Group P622 not implemented\n"; exit(1);
        break;
     }
    		  

}//symmetrize_crystal_vectors end			         
#define put_inside(j,j_min,j_max,jint)  \
   if( (j) < (j_min) ) { (j) = (j) + (jint);}\
   else if( (j) > (j_max) ) { (j) = (j) - (jint);};

/* Symmetrizes a simple grid with P4  symmetry --------------------------*/
     void symmetry_P4(Volume &vol, const SimpleGrid &grid,
                      const matrix1D<double> &eprm_aint, 
		      const matrix1D<double> &eprm_bint,
		      const matrix2D<int> &mask, int volume_no, int grid_type)
{
   int ZZ_lowest =(int) ZZ(grid.lowest);
   int YY_lowest =STARTINGY(mask);
   int XX_lowest =STARTINGX(mask);
   int ZZ_highest=(int) ZZ(grid.highest);
   int YY_highest=FINISHINGY(mask);
   int XX_highest=FINISHINGX(mask);
   
   while(1)
     {
      if(mask(0,XX_lowest)==0) XX_lowest++;
      else
         break;
      if(XX_lowest==XX_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(0,XX_highest)==0) XX_highest--;
      else
         break;
      if(XX_lowest==XX_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(YY_lowest,0)==0) YY_lowest++;
      else
         break;
      if(YY_lowest==YY_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(YY_highest,0)==0) YY_highest--;
      else
         break;
      if(YY_lowest==YY_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
	 
     }

//   int ZZ_lowest =(int) ZZ(grid.lowest);
//   int YY_lowest =MAX((int) YY(grid.lowest),STARTINGY(mask));
//   int XX_lowest =MAX((int) XX(grid.lowest),STARTINGX(mask));
//   int ZZ_highest=(int) ZZ(grid.highest);
//   int YY_highest=MIN((int) YY(grid.highest),FINISHINGY(mask));
//   int XX_highest=MIN((int) XX(grid.highest),FINISHINGX(mask));

   int maxZ, maxY, maxX, minZ, minY, minX;
   int x,y,z;

   int x0,y0,z0;
   int x1,y1,z1;
   int x2,y2,z2;

   int XXaint, YYbint;
   XXaint = (int) XX(eprm_aint); YYbint = (int)YY(eprm_bint);

   int XXaint_2, YYbint_2;
   XXaint_2 = XXaint/2; YYbint_2 = YYbint/2;
   int xx,yy,zz;
   
   if( ABS(XX_lowest) > ABS(XX_highest)) 
       { minX=XX_lowest; maxX=0;}
   else 
       { minX=0; maxX=XX_highest;}
   if( ABS(YY_lowest) > ABS(YY_highest)) 
       { minY=YY_lowest; maxY=0;}
   else 
       { minY=0; maxY=YY_highest;}

   minZ=ZZ_lowest;
   maxZ=ZZ_highest;
   //FCC non supported yet
   if(volume_no==1 && grid_type==FCC)   
      {cerr<< "\nSimetries using FCC not implemented\n";exit(1);}
   for (z=minZ;z<=maxZ;z++) 
     for (y=minY;y<=maxY;y++) 
       for (x=minX;x<=maxX;x++) {
	   //sym=-1---------------------------------------------------------
           if (!MAT_ELEM(mask,y,x) || z < ZZ_lowest || z > ZZ_highest)
             continue;

	   //sym=0 ---------------------------------------------------------
	   xx =-y; yy= x; zz=z;
           //only the first simple grid center and the origen is the same 
	   //point. 
	   if(volume_no==1) {xx--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }	     
           x0=xx;y0=yy;z0=zz;
	   
	   //sym=1----------------------------------------------------------
	   xx = -x; yy= -y; zz= z;
	   if(volume_no==1) {xx--;yy--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=1 end	     
           x1=xx;y1=yy;z1=zz;

	   //sym=2----------------------------------------------------------
	   xx = y; yy= -x; zz= z;
	   if(volume_no==1) {yy--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=2 end	     
           x2=xx;y2=yy;z2=zz;

           VOLVOXEL(vol,z ,y ,x ) = VOLVOXEL(vol,z0,y0,x0) =  
           VOLVOXEL(vol,z1,y1,x1) = VOLVOXEL(vol,z2,y2,x2) = 
	  (VOLVOXEL(vol,z ,y ,x ) + VOLVOXEL(vol,z0,y0,x0) +
	   VOLVOXEL(vol,z1,y1,x1) + VOLVOXEL(vol,z2,y2,x2) )/4.0;
	   }//for end

}				
/* Symmetrizes a simple grid with P4212 symmetry--------------------------*/

     void symmetry_P42_12(Volume &vol, const SimpleGrid &grid,
				const matrix1D<double> &eprm_aint, 
			        const matrix1D<double> &eprm_bint,
				const matrix2D<int> &mask, int volume_no,
				int grid_type)
{

   int ZZ_lowest =(int) ZZ(grid.lowest);
   int YY_lowest =STARTINGY(mask);
   int XX_lowest =STARTINGX(mask);
   int ZZ_highest=(int) ZZ(grid.highest);
   int YY_highest=FINISHINGY(mask);
   int XX_highest=FINISHINGX(mask);

   //if there is an extra slice in the z direction there is no way
   //to calculate the -z slice
   if(ABS(ZZ_lowest) >ABS(ZZ_highest) ) ZZ_lowest= -(ZZ_highest);
   else ZZ_highest= ABS(ZZ_lowest);
   
   while(1)
     {
      if(mask(0,XX_lowest)==0) XX_lowest++;
      else
         break;
      if(XX_lowest==XX_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(0,XX_highest)==0) XX_highest--;
      else
         break;
      if(XX_lowest==XX_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(YY_lowest,0)==0) YY_lowest++;
      else
         break;
      if(YY_lowest==YY_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(YY_highest,0)==0) YY_highest--;
      else
         break;
      if(YY_lowest==YY_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
	 
     }

//   int ZZ_lowest =(int) ZZ(grid.lowest);
//   int YY_lowest =MAX((int) YY(grid.lowest),STARTINGY(mask));
//   int XX_lowest =MAX((int) XX(grid.lowest),STARTINGX(mask));
//   int ZZ_highest=(int) ZZ(grid.highest);
//   int YY_highest=MIN((int) YY(grid.highest),FINISHINGY(mask));
//   int XX_highest=MIN((int) XX(grid.highest),FINISHINGX(mask));

   int maxZ, maxY, maxX, minZ, minY, minX;
   int x,y,z;

   int x0,y0,z0;
   int x1,y1,z1;
   int x2,y2,z2;
   int x3,y3,z3;
   int x4,y4,z4;
   int x5,y5,z5;
   int x6,y6,z6;

   int XXaint, YYbint;
   XXaint = (int) XX(eprm_aint); YYbint = (int)YY(eprm_bint);

   int XXaint_2, YYbint_2;
   XXaint_2 = XXaint/2; YYbint_2 = YYbint/2;
   int xx,yy,zz;
   
   if( ABS(XX_lowest) > ABS(XX_highest)) 
       { minX=XX_lowest; maxX=0;}
   else 
       { minX=0; maxX=XX_highest;}
   if( ABS(YY_lowest) > ABS(YY_highest)) 
       { minY=YY_lowest; maxY=0;}
   else 
       { minY=0; maxY=YY_highest;}
   if( ABS(ZZ_lowest) > ABS(ZZ_highest)) 
       { minZ=ZZ_lowest; maxZ=0;}
   else 
       { minZ=0; maxZ=ZZ_highest;}

   //FCC non supported yet
   if(volume_no==1 && grid_type==FCC)   
      {cerr<< "\nSimetries using FCC not implemented\n";exit(1);}
       
   for (z=minZ;z<=maxZ;z++) 
     for (y=minY;y<=maxY;y++) 
       for (x=minX;x<=maxX;x++) {
	   //sym=-1---------------------------------------------------------
           if (!MAT_ELEM(mask,y,x) || z < ZZ_lowest || z > ZZ_highest)
             continue;

	   //sym=0 ---------------------------------------------------------
	   xx =-x; yy=-y; zz=z;
	   if(volume_no==1) { xx--;yy--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }	     
           x0=xx;y0=yy;z0=zz;
	   
	   //sym=1----------------------------------------------------------
	   xx =y; yy=x; zz= -z;//I think z-- is always inside the grid
	                       //we do not need to check
	   if(volume_no==1) { zz--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=1 end	     
           x1=xx;y1=yy;z1=zz;

	   //sym=2----------------------------------------------------------
	   xx = -y; yy= -x; zz= -z;
	   if(volume_no==1) { xx--;yy--;zz--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=2 end	     
           x2=xx;y2=yy;z2=zz;

	   //sym=3----------------------------------------------------------
	   xx = y+XXaint_2; yy= -x+YYbint_2; zz=z;
	   if(volume_no==1) { yy--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=3 end	     
           x3=xx;y3=yy;z3=zz;

	   //sym=4----------------------------------------------------------
	   xx = -y+XXaint_2; yy= +x+YYbint_2; zz=z;
	   if(volume_no==1) { xx--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=4 end	     
           x4=xx;y4=yy;z4=zz;
	   //sym=5----------------------------------------------------------
	   xx = -x+XXaint_2; yy= +y+YYbint_2; zz= -z;
	   if(volume_no==1) { xx--;zz--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=5 end	     
           x5=xx;y5=yy;z5=zz;

	   //sym=6----------------------------------------------------------
	   xx = +x+XXaint_2; yy= -y+YYbint_2; zz= -z;
	   if(volume_no==1) { yy--;zz--;}
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=6 end	      	      	       
           x6=xx;y6=yy;z6=zz;

           //only the first simple grid center and the origen is the same 
	   //point. 
//	   if(volume_no==1)
//	   {
// 	    switch (grid_type) {
//		case FCC: cerr<< "\nSimetries using FCC not implemented\n";break;
//		case BCC: x0--;y0--;               z1--; 
//		          x2--;y2--;z2--;     y3--;
//			  x4--;          x5--;     z5--;
//			       y6--;z6--;
//		          break;
//		case CC:  break;
//	     }
//	   }

           VOLVOXEL(vol,z ,y ,x ) = VOLVOXEL(vol,z0,y0,x0) =  
           VOLVOXEL(vol,z1,y1,x1) = VOLVOXEL(vol,z2,y2,x2) = 
           VOLVOXEL(vol,z3,y3,x3) = VOLVOXEL(vol,z4,y4,x4) = 
           VOLVOXEL(vol,z5,y5,x5) = VOLVOXEL(vol,z6,y6,x6) = 
	  (VOLVOXEL(vol,z ,y ,x ) + VOLVOXEL(vol,z0,y0,x0) +
	   VOLVOXEL(vol,z1,y1,x1) + VOLVOXEL(vol,z2,y2,x2) +
	   VOLVOXEL(vol,z3,y3,x3) + VOLVOXEL(vol,z4,y4,x4) +
	   VOLVOXEL(vol,z5,y5,x5) + VOLVOXEL(vol,z6,y6,x6))/8.0;
	   }//for end
}//symmetryP42_12 end				
/* Symmetrizes a simple grid with P6 symmetry-----------------------------*/
     void symmetry_P6(Volume &vol, const SimpleGrid &grid,
                                const matrix1D<double> &eprm_aint, 
			        const matrix1D<double> &eprm_bint,
				const matrix2D<int> &mask, int volume_no,
				int grid_type)
{

   int ZZ_lowest =(int) ZZ(grid.lowest);
   int YY_lowest =STARTINGY(mask);
   int XX_lowest =STARTINGX(mask);
   int ZZ_highest=(int) ZZ(grid.highest);
   int YY_highest=FINISHINGY(mask);
   int XX_highest=FINISHINGX(mask);

   
   while(1)
     {
      if(mask(0,XX_lowest)==0) XX_lowest++;
      else
         break;
      if(XX_lowest==XX_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(0,XX_highest)==0) XX_highest--;
      else
         break;
      if(XX_lowest==XX_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(YY_lowest,0)==0) YY_lowest++;
      else
         break;
      if(YY_lowest==YY_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
     }
   while(1)
     {
      if(mask(YY_highest,0)==0) YY_highest--;
      else
         break;
      if(YY_lowest==YY_highest) 
         {
	 cerr << "Error in symmetry_P42_12, while(1)" << endl;
	 exit(0);
	 }
	 
     }

//   int ZZ_lowest =(int) ZZ(grid.lowest);
//   int YY_lowest =MAX((int) YY(grid.lowest),STARTINGY(mask));
//   int XX_lowest =MAX((int) XX(grid.lowest),STARTINGX(mask));
//   int ZZ_highest=(int) ZZ(grid.highest);
//   int YY_highest=MIN((int) YY(grid.highest),FINISHINGY(mask));
//   int XX_highest=MIN((int) XX(grid.highest),FINISHINGX(mask));

   int maxZ, maxY, maxX, minZ, minY, minX;
   int x,y,z;

   int x0,y0,z0;
   int x1,y1,z1;
   int x2,y2,z2;
   int x3,y3,z3;
   int x4,y4,z4;

   int XXaint, YYbint;
   XXaint = (int) XX(eprm_aint); YYbint = (int)YY(eprm_bint);

   int XXaint_2, YYbint_2;
   XXaint_2 = XXaint/2; YYbint_2 = YYbint/2;
   int xx,yy,zz;
   
   if( ABS(XX_lowest) > ABS(XX_highest)) 
       { minX=XX_lowest; maxX=0;}
   else 
       { minX=0; maxX=XX_highest;}
   //P6 is tricky. I have decide to apply it to half the volume
   //instead of to 1 sizth. I think the amount of ifs that I save
   //are worth this larger loop
   
   minY=YY_lowest;
   maxY=YY_highest;

   minZ=ZZ_lowest; 
   maxZ=ZZ_highest;

   //FCC non supported yet
   if(volume_no==1 && grid_type==FCC)   
      {cerr<< "\nSimetries using FCC not implemented\n";exit(1);}
       
   for (z=minZ;z<=maxZ;z++) 
     for (y=minY;y<=maxY;y++) 
       for (x=minX;x<=maxX;x++) {
	   //sym=-1---------------------------------------------------------
           if (!MAT_ELEM(mask,y,x) || z < ZZ_lowest || z > ZZ_highest)
             continue;

	   //sym=0 ---------------------------------------------------------
	   xx = x - y; yy= x; zz=z;
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }	     
           x0=xx;y0=yy;z0=zz;
	   
	   //sym=1----------------------------------------------------------
	   xx = -y; yy= x - y; zz= z;
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=1 end	     
           x1=xx;y1=yy;z1=zz;

	   //sym=2----------------------------------------------------------
	   xx = -x; yy= -y; zz= z;
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=2 end	     
           x2=xx;y2=yy;z2=zz;

	   //sym=3----------------------------------------------------------
	   xx = -x + y; yy= -x; zz=z;
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=3 end	     
           x3=xx;y3=yy;z3=zz;

	   //sym=4----------------------------------------------------------
	   xx = +y; yy= -x+y; zz=z;
	   if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx) )
	      {
	      put_inside(xx,XX_lowest,XX_highest,XXaint)
	      put_inside(yy,YY_lowest,YY_highest,YYbint)
	      if (!MAT_ELEM(mask,yy,xx) || mask.outside(yy,xx))
	          cerr << "ERROR in symmetry_P function"
		       << "after correction spot is still"
		       << "outside mask\a"<<endl;
	      }//sym=4 end	     
           x4=xx;y4=yy;z4=zz;

	   if(volume_no==1)
	   {
 	    switch (grid_type) {
		case FCC: cerr<< "\nSimetries using FCC not implemented\n";break;
	  //there is no way to reinforce P6 in the second grid without 
	  //interpolation. This is the best we can do.
		case BCC: x1=x0=x;y1=y0=y;
		          x2=x3=x4= -x-1;y2=y3=y4= -y-1;
		          break;
		case CC:  break;
	     }
	   }


           VOLVOXEL(vol,z ,y ,x ) = VOLVOXEL(vol,z0,y0,x0) =  
           VOLVOXEL(vol,z1,y1,x1) = VOLVOXEL(vol,z2,y2,x2) = 
           VOLVOXEL(vol,z3,y3,x3) = VOLVOXEL(vol,z4,y4,x4) = 
	  (VOLVOXEL(vol,z ,y ,x ) + VOLVOXEL(vol,z0,y0,x0) +
	   VOLVOXEL(vol,z1,y1,x1) + VOLVOXEL(vol,z2,y2,x2) +
	   VOLVOXEL(vol,z3,y3,x3) + VOLVOXEL(vol,z4,y4,x4))/6.0;

	   }//for end

}				
#undef wrap_as_Crystal
#undef DEBUG
