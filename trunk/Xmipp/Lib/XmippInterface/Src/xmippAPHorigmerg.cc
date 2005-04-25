/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
 *
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

//#include <vector>

#include "../xmippAPHorigmerg.hh"
#include <XmippData/xmippArgs.hh>
#include <fstream>
#include <iomanip>

#define VERBOSE
//#define DEBUG
// APH =====================================================================
void APHFileorigmerg::read(const FileName &fn,
                                           const int mrc_label) _THROW {
   ifstream  fh_aph;
   string    line;
   // Empties current APH File
   clear();
   read_mrc_label=mrc_label;
   // Open file
   fh_aph.open(fn.c_str(), ios::in);
   if (!fh_aph)
      REPORT_ERROR(1601,"APHFileorigmerg::read: File "+fn+" not found");

   // Read first line and skip it
   getline(fh_aph,line);
   int h,k;
   int line_no=0;
   double FILM;
   spot tmp_spot;
   while (!fh_aph.eof()) {
      try {
         line_no++;
         getline(fh_aph,line);
         if(line.size()==0) continue;	 
	 h = tmp_spot.h	      = AtoI(first_token(line));
	 k = tmp_spot.k	      = AtoI(next_token());
	 tmp_spot.zstar       = AtoF(next_token());
	 tmp_spot.amp	      = AtoF(next_token());
	 tmp_spot.phase       = AtoF(next_token());
	 FILM =tmp_spot.FILM  = AtoI(next_token());
	 tmp_spot.IQ	      = AtoI(next_token());
	 tmp_spot.FLMWGT      = AtoF(next_token());
	 tmp_spot.BACK	      = AtoF(next_token());
	 tmp_spot.CTF	      = AtoF(next_token());
         if (mrc_label >= 0 && (FILM != mrc_label))
	     continue; 
	 aph_data_vector.push_back(tmp_spot);
	 max_h=MAX(max_h,h);
	 max_k=MAX(max_k,k);
	 min_h=MIN(min_h,h);
	 min_k=MIN(min_k,k);
      }
      catch (Xmipp_error) {
         cout << "aph File: Line " << line_no << " is skipped due to an error\n";
      }
   }/* while */

//   #define DEBUG_max   
   #ifdef DEBUG_max   
   cout << "max_h: " << max_h <<" max_k: " << max_k << endl;
   cout << "min_h: " << min_h <<" min_k: " << min_k << endl;
   cout << "Number Spots: " << aph_data_vector.size()<< endl;
   #endif      
   #undef DEBUG_max
   
   // remember to clear if you need to reread the file
   fh_aph.close();
   //fh_aph.clear();
   //copy file name
   fn_aph=fn;

}/*  APHFile::read */

/* ------------------------------------------------------------------------- */
void APHFileorigmerg::write(const FileName &fn) const _THROW {
   ofstream fh;
   fh.open(fn.c_str());
   char aux_char[128];
   if (!fh)
      REPORT_ERROR(1,(string)"APHFileorigmerg::write: Cannot open "+
         fn+" for output");
   
   fh << setfill('0') << setw(4) << read_mrc_label << endl;
   for(int line_no = 0; line_no < aph_data_vector.size(); line_no++){
//           1X,2I4,F8.4,F10.1,F7.1,I7,I3,F8.5,F10.1,F7.3
      sprintf(aux_char," %4d%4d%8.4f%10.1f%7.1f%7d%3d%8.5f%10.1f%7.3f\n",
			       (aph_data_vector[line_no]).h,	
			       (aph_data_vector[line_no]).k,	
			       (aph_data_vector[line_no]).zstar,
			       (aph_data_vector[line_no]).amp,	
			       (aph_data_vector[line_no]).phase,
			       (aph_data_vector[line_no]).FILM,
			       (aph_data_vector[line_no]).IQ,	
			       (aph_data_vector[line_no]).FLMWGT,
			       (aph_data_vector[line_no]).BACK,
			       (aph_data_vector[line_no]).CTF);
      fh <<aux_char;
   }
   fh.close();
}


/* ------------------------------------------------------------------------- */
void APHFileorigmerg::clear(){
   read_mrc_label=-1;
   max_h = -100000;
   max_k = -100000;
   min_h = 100000;
   min_k = 100000;
} /*clear*/

/** Move spots from assymetric unit to the plane  h>0
    those points that do not fit in the plane are ignored
    Will be used in the future. Need a,b magnitude (real space A) 
    taxa,tilt.anglefrom a to  (radians).
    
    Generate symmetrical points:
    Symmetrical points are generated through the multiplication
    \begin{verbatim}
   				  [R[0] R[2]  0    0 R[5]] 
   				  [R[1] R[3]  0    0 R[6]]
    [h' k' l' A' ph']=[h k l A ph][ 0	 0   R[4]  0  0  ]
   				  [ 0	 0    0    1  0  ]
   				  [ 0	 0    0    0 R[7]]
    \end{verbatim}
   
    These R matrices are obtained for each crystallographic group
    They can be easily obtaine from  origtiltd.for (MRC source code)
   */
   void APHFileorigmerg::unasymmetrization(const double a_mag,const  double b_mag, 
                          const double mrc_taxa,const  double mrc_tilt,
			  const double a_b_ang,const  int symmetry_group,
                          matrix2D<int> &Counter)   
  {
   int ksize=MAX(ABS(min_k),ABS(max_k));
   int hsize=MAX(ABS(min_h),ABS(max_h));
   Counter.init_zeros(2*ksize+1,2*hsize+1);
   STARTINGY(Counter)=-ksize;
   STARTINGX(Counter)=-hsize;
   //Set right symmetrization matrices and unsymmetrice vector
   if(read_mrc_label<0)
      REPORT_ERROR(1601,"APHFileorigmerg::unasymmetrization: Not implemented for several micrographs");
   switch( symmetry_group )
     {  
 	
        case sym_P1     :	//Nothing to do here
	   break;			 
        case sym_P2_122 :     // must do P22_12 but I do not think
	                      // it will be used but for debuging (ROB april
			      // 2005)
	   unsymmetrice_P222_1( a_mag,  b_mag, 
                                mrc_taxa,  mrc_tilt,
			        a_b_ang,  symmetry_group, Counter); 
	   break;
        default         : 
	   cerr << "APHFileorigmerg::unasymmetrization:" << endl;
	   cerr << "\t\tsymmetry " << symmetry_group << " not implemented" 
	        << endl;
           exit(1);
     }  
  }//unasymmetrization end
/////////////////////////////////////////////////////////////////////////

    /* Possible Rs for this group **SEEMS to be P22_12 NOT P2_122  rob april
    2005**
       h > 0
	I:[ 1 0 0  1  1 0    0  1] -> ( h, k, l)=(h,k,l)
       R2:[ 1 0 0  1 -1 0    0 -1] -> ( h, k,-l)=conj(h,k,l)
       R3:[ 1 0 0 -1  1 0  180 -1] -> ( h,-k, l)=conj(h,k,l)*(-1)^k
    R3.R2:[ 1 0 0 -1 -1 0 -180  1] -> ( h,-k,-l)=(h,k,l)*(-1)^k
       h < 0
       R1:[-1 0 0 -1 -1 0    0 -1] -> (-h,-k,-l)=conj(h,k,l)
    R1.R2:[-1 0 0 -1  1 0    0  1] -> (-h,-k, l)=(h,k,l)
    R1.R3:[-1 0 0  1 -1 0 -180  1] -> (-h, k,-l)=(h,k,l)*(-1)^k
 R1.R2.R3:[-1 0 0  1  1 0 -180 -1] -> (-h, k, l)=conj(h,k,l)*(-1)^k

       The special cases of this group are
       Real         Imag
       ----         ----
       (0,2n,l)     (0,2n+1,l)
       (h,k,0)
       (h,0,l)

       I will only use those restrictions valid for all l.
    */
   
void APHFileorigmerg::unsymmetrice_P222_1(const double a_mag,const  double b_mag, 
                               const double mrc_taxa,const  double mrc_tilt,
			       const double a_b_ang,const  int symmetry_group,
			       matrix2D<int> &Counter)
{
   // Go through all the points and compute the three alternatives
   spot first,second,third,fourth;
   spot aux_spot;
   double diff_z[4];
   // NOTE: TAXA TILT ARE IN RADIANS BUT PHASE IS IN DEG, 
   // Compute Z
   
   double h_contrib,  k_contrib, minimun_diff;
   int winner;
   int h,k;
   compute_Z( a_mag,  b_mag,mrc_taxa, mrc_tilt,
	      a_b_ang,  h_contrib,  k_contrib);

   //Search for -IQ and invert contrast
   for(int line_no = (aph_data_vector.size()-1); 
           line_no >=0; line_no--){
      if((aph_data_vector[line_no]).IQ < 0){
          (aph_data_vector[line_no]).h *= -1;
          (aph_data_vector[line_no]).k *= -1;
          (aph_data_vector[line_no]).zstar *= -1.0;
          (aph_data_vector[line_no]).phase *= -1.0;
      }	  
   }
   //h=0 is an special case even if the original values of k
   //were negative that information is lost for ever
   //fortunatelly if a and b are selectec positive in
   //the h>0 plane. k will be positive
   //so do nothing

   /////////////////////
   /////////////////////

   //Now we should loop through all the data and see if there are
   //spots with two values asigned   
   for(int line_no=0; line_no < (aph_data_vector.size()-1); 
           line_no++){
      k = (aph_data_vector[line_no]).k;
      h = (aph_data_vector[line_no]).h;
      MAT_ELEM(Counter, k, h) += 1;
      //check for errors
      if ( MAT_ELEM(Counter, k, h)>2)
           REPORT_ERROR(1234,"unsymmetrice_P222_1: Count can not be > 2.");
   }
   //#define DEBUGCounter
   #ifdef DEBUGCounter
   MAT_ELEM(Counter,0,0)=9999;
   MAT_ELEM(Counter,0,1)=1111;
   #endif
   #undef DEBUGCounter
   
//   #define DEBUGCounter1
   #ifdef DEBUGCounter1
   cout << Counter;   
   #endif
   #undef DEBUGCounter1

   //For all those positions with two points the value of Z should allow us to 
   //relocate them. OF course this is not valid for mrc_tilt=0
   
   //Special case mrc_tilt=0
   if(mrc_tilt==0){
//      for(int line_no = 0; 
//           line_no < aph_data_vector.size(); line_no++){
      for(int line_no = (aph_data_vector.size()-1); 
              line_no >=0; line_no--){
	 //if Counter=1 there is no way to decide the origin of this point
	 //h must be positive because we are in P2221
	 //so apply R3. (Note, we will never recover h<0)
	 
	 //I am assuming that the spot that should is the first in the aph file
	 //since it is more likelly that k < 0
	 //Once one spot is modified I do not want to modify its pair
         k = (aph_data_vector[line_no]).k;
         h = (aph_data_vector[line_no]).h;
	 if ( MAT_ELEM(Counter, k, h)> 1){
             (aph_data_vector[line_no]).k  *= -1;
	     (aph_data_vector[line_no]).phase   = 
	                           -1.* (aph_data_vector[line_no]).phase+
	                          180.0*k;			       
             MAT_ELEM(Counter, k, h) -= 1;	
	     //#define DEBUGTILTZERO
	     #ifdef  DEBUGTILTZERO
	     cout << (aph_data_vector[line_no]);
	     #endif
	     #undef DEBUGTILTZERO
	 }//if ( MAT_ELEM(Counter, k, h)> 1){   
      }//for(int line_no
   }//mrc_tilt==0
   //warp the phase (this should go at the end
   for(int line_no = (aph_data_vector.size()-1); 
           line_no >=0; line_no--)
      (aph_data_vector[line_no]).phase = 
               realWRAP((aph_data_vector[line_no]).phase,-180.,180.);

#ifdef NEVERDEFINED
   for(int line_no = (aph_data_vector.size()-1); 
           line_no >=0; line_no--){
      first = aph_data_vector[line_no];
      //IS h positive?
      fourth=third=second=first=aph_data_vector[line_no];
      if((aph_data_vector[line_no]).h>0){
	 second.phase *= -1.; second.zstar *= -1.;
	 third.phase   = -1.*third.phase+180.0*third.k;third.k  *= -1;
	 fourth.phase  =    fourth.phase-180.0*fourth.k;
	     fourth.k *= -1;fourth.zstar     *= -1.;
      }
      //if h=0 k must be positive and as the input k is positive
      //we allow transformation that leave k unchange so we
      //should leave the popint unchanged unless IQ < 0
      //in this case -h and -k should be aplied and phase
      //kept the saem. Do not know about l yet but I think must be
      //negated also
      else if((aph_data_vector[line_no]).h==0){
	 third.phase *= -1.; third.zstar *= -1.;
	 second.phase = second.phase-180.0*second.k;
	 second.h  = 0;second.zstar  *= -1.;
	 fourth.phase  = -1.*fourth.phase-180.0*fourth.k;
	 fourth.h = 0;
      }
      else{
	 first.h  *= -1;  first.k *= -1;  first.zstar *= -1.;
	                  first.phase *= -1.;
	 second.h *= -1;second.k *= -1;
	 third.phase = third.phase-180.0*third.k;
	 third.h  *= -1;third.zstar  *= -1.;
	 fourth.phase  = -1.*fourth.phase-180.0*fourth.k;
	 fourth.h *= -1;
      }      
//      #define DEBUG
      #ifdef DEBUG
      cout << "Theorical_Z 1 =" << first.h  *  h_contrib + 
                                 first.k  *  k_contrib << endl;
      cout << first ;
      cout << "Theorical_Z 2=" << second.h  *  h_contrib + 
                                 second.k  *  k_contrib << endl;
      cout << second ;
      cout << "Theorical_Z 3=" << third.h  *  h_contrib + 
                                 third.k  *  k_contrib << endl;
      cout << third ;
      cout << "Theorical_Z 4=" << fourth.h  *  h_contrib + 
                                 fourth.k  *  k_contrib << endl;
      cout << fourth ;
      #endif
      #undef DEBUG
      // Select  right point 
      diff_z[0]=(ABS(first.zstar- (first.h  *  h_contrib + 
                                         first.k  *  k_contrib)));
      diff_z[1]=(ABS(second.zstar-(second.h * h_contrib + 
                                         second.k * k_contrib)));
      diff_z[2]=(ABS(third.zstar- (third.h  * h_contrib + 
                                         third.k  * k_contrib)));
      diff_z[3]=(ABS(fourth.zstar-(fourth.h * h_contrib + 
                                         fourth.k * k_contrib)));
      minimun_diff = diff_z[0]; aux_spot=first;
      if (minimun_diff > diff_z[1]) {aux_spot=second;
                                     minimun_diff = diff_z[1];}
      if (minimun_diff > diff_z[2]) {aux_spot=third;
                                     minimun_diff = diff_z[2];}
      if (minimun_diff > diff_z[3]) {aux_spot=fourth;
                                     minimun_diff = diff_z[3];}

      //record how many spots are allocated in each point
      k = aux_spot.k;
      h = aux_spot.h;
      MAT_ELEM(Counter, k, h) += 1;
      //if mrc_tilt=0 them average spots. the z criterium is not valid since
      // z= -z
      if(mrc_tilt==0 && h!=0)
         {
         if ( MAT_ELEM(Counter, k, h)> 1)
	     {
	      //if IQ positive h must be positive and z=0
	      //so I must apply R3 (o R3 .R2) (h , -k, l,am,180*k-ph
	      //otherwise R1 . R3 (-h,k,-l,Am -180k +ph)
//cout << "(h,k,aux_spot.IQ) " << h << " " << k << " " << aux_spot.IQ << endl;	      
		 if(aux_spot.IQ>0){
        	    aux_spot.k  *= -1;
		    aux_spot.phase   = -1.*aux_spot.phase+
	                	      180.0*aux_spot.k;
        	 }
		 else{
        	    aux_spot.h  *= -1;
		    aux_spot.phase   = aux_spot.phase -
	                	      180.0*aux_spot.k;
        	 } 
		 //have no idea what will happend when h and k < 0
		MAT_ELEM(Counter, k, h) -= 1;
		if ( MAT_ELEM(Counter, aux_spot.k, aux_spot.h) >2)
	           REPORT_ERROR(1234,"unsymmetrice_P222_1: Count can not be =1.");
		MAT_ELEM(Counter, aux_spot.k, aux_spot.h) += 1;
		}// if ( MAT_ELEM(Counter, k, h)> 1)

         }//if(mrc_tilt==0)
	aux_spot.phase = realWRAP(aux_spot.phase,-180.,180.);
      aph_data_vector[line_no]=aux_spot;
cout << "Selected point: " << aph_data_vector[line_no];
      if ( MAT_ELEM(Counter, k, h)>2)
           REPORT_ERROR(1234,"unsymmetrice_P222_1: Count can not be > 2.");
//      #define DEBUG
      #ifdef DEBUG
      cout << "diff_z[0]" << diff_z[0] << " diff_z[1]" << diff_z[1] << endl;
      cout << "diff_z[2]" << diff_z[2] << " diff_z[3]" << diff_z[3] << endl;
      
      cout << "diff_z[0] again " <<first.zstar- (first.h  *  h_contrib + 
                                         first.k  *  k_contrib) <<endl;
      //cout << aph_data_vector[line_no] ;
      #endif
      #undef DEBUG
   }//for line_no
   //#define DEBUG
#endif//NEVERDEEFINED   
}//unsymmetrice_P222_1			       
void APHFileorigmerg::compute_Z(double a_mag, double b_mag, 
               double mrc_taxa, double mrc_tilt,
	       double a_b_ang, double & h_contrib, double & k_contrib)
	       
{   
    if(mrc_tilt==0) {h_contrib = k_contrib = 0.0;return;}
    //Not sure this will work with non square or hexagonal crystal
    // mod_a, mod_b, a_b_angle in real space
    double mod_astar_in_A =1./a_mag;//Ansgtroms ^-1
    double mod_bstar_in_A =1./b_mag;
    mod_astar_in_A /=sin(PI-a_b_ang);//angle in reciprocal space -> Pi-angle
    mod_bstar_in_A /=sin(PI-a_b_ang);//angle in reciprocal space
    
    double taxb = mrc_taxa + a_b_ang; // angle from tilt axis to bstar
    //astar and bstar component perpendicular to tilt axis
    double per_tilt_a = mod_astar_in_A * sin(mrc_taxa);
    double per_tilt_b = mod_bstar_in_A * sin(taxb);
    double tilt_tangent = tan(mrc_tilt);
    h_contrib = /* -1.  * */ tilt_tangent * per_tilt_a;
    k_contrib = /* -1.  * */ tilt_tangent * per_tilt_b;
    //#define DEBUG
    #ifdef DEBUG
    cout << "mod_astar_in_A: " << mod_astar_in_A << endl;
    cout << "mod_bstar_in_A: " << mod_bstar_in_A << endl;
    cout << "per_tilt_a staxa: " << per_tilt_a << endl;//staxa
    cout << "per_tilt_b staxb: " << per_tilt_b << endl;//staxb
    cout << "a_b_ang_star  : " << RAD2DEG(a_b_ang) << endl;
    cout << "taxb          : " << RAD2DEG(taxb) << endl;
    cout << "ttilt         : " <<  tilt_tangent;
    double z_star = h_contrib;
    cout << " z(hk 10) " <<  z_star << endl;
    z_star = k_contrib;
    cout << " z(hk 01) " <<  z_star << endl;
    #endif
    #undef DEBUG    

}//compute_Z  
