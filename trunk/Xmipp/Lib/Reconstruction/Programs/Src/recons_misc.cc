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

#include "../recons_misc.hh"

/* Fill Reconstruction info structure -------------------------------------- */
void build_recons_info(SelFile &selfile, SelFile &selctf, const FileName &fn_ctf,
   const SymList &SL,
   Recons_info * &IMG_Inf) {
   matrix2D<double>  L(4,4), R(4,4);    // A matrix from the list
   FileName          fn_proj;
   FileName          fn_ctf1; 
   Projection        read_proj;
   bool              is_there_ctf=false;
   bool              is_ctf_unique=false;

   int trueIMG=selfile.ImgNo();
   selfile.go_first_ACTIVE();
   int numIMG=trueIMG*(SL.SymsNo() + 1);
   
   // The next two ifs check whether there is a CTF file, and 
   // whether it is unique 
   if (fn_ctf!="") {
      is_there_ctf=1;
      is_ctf_unique=1;
      if (Is_FourierImageXmipp(fn_ctf)==FALSE) {
         is_ctf_unique=0;
         int truectf=selctf.ImgNo();
         selctf.go_first_ACTIVE();
         int numctf=truectf*(SL.SymsNo() + 1);
         if (numctf!=numIMG) 
		REPORT_ERROR(9696,"Number of CTF and image files differ");
      }
   }

   if (IMG_Inf!=NULL) delete [] IMG_Inf;
   if ((IMG_Inf = new Recons_info[numIMG])== NULL)
      REPORT_ERROR(3008,"Build_Recons_Info: No memory for the sorting");

   int i = 0; // It will account for the number of valid projections processed
   cerr << "Reading angle information ...\n";
   init_progress_bar(trueIMG);
   while (!selfile.eof()) {
     fn_proj=selfile.NextImg();
     if ((is_there_ctf==1)&&(is_ctf_unique!=1)) fn_ctf1=selctf.NextImg(); 
     if (fn_proj!="") {
        read_proj.read(fn_proj);

        // Filling structure
        IMG_Inf[i].fn_proj = fn_proj;
        if      (is_ctf_unique) IMG_Inf[i].fn_ctf = fn_ctf;
        else if (is_there_ctf)  IMG_Inf[i].fn_ctf = fn_ctf1;
        IMG_Inf[i].sym     = -1;
        read_proj.get_eulerAngles(IMG_Inf[i].rot,IMG_Inf[i].tilt,IMG_Inf[i].psi);
        EULER_CLIPPING(IMG_Inf[i].rot,IMG_Inf[i].tilt,IMG_Inf[i].psi);

        // Any symmetry?
        if (SL.SymsNo()>0) {
           for (int j=0; j<SL.SymsNo(); j++) {
               int sym_index=SYMINDEX(SL,j,i,trueIMG);
               IMG_Inf[sym_index].fn_proj=IMG_Inf[i].fn_proj;
               if (is_ctf_unique)     IMG_Inf[sym_index].fn_ctf = fn_ctf;
	       else if (is_there_ctf) IMG_Inf[sym_index].fn_ctf = IMG_Inf[i].fn_ctf;
               IMG_Inf[sym_index].sym=j;
               SL.get_matrices(j,L,R);
               L.resize(3,3); // Erase last row and column
               R.resize(3,3); // as only the relative orientation
                              // is useful and not the translation
               double drot, dtilt, dpsi;
      	       Euler_apply_transf(L,R,
      	          IMG_Inf[i].rot,IMG_Inf[i].tilt,IMG_Inf[i].psi,
                  drot, dtilt, dpsi);
               IMG_Inf[sym_index].rot=(float)drot;
               IMG_Inf[sym_index].tilt=(float)dtilt;
      	       IMG_Inf[sym_index].psi=(float)dpsi;
           }
        }
     }

     i ++; // I have processed one more image
     if (i%25==0) progress_bar(i);
   }
   progress_bar(trueIMG);
}
