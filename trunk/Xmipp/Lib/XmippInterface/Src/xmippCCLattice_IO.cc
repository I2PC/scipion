/***************************************************************************
 *
 * Authors:     Debora Gil
 *              Roberto Marabini
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
 *                                      <                               
 * You should have received a copy of the GNU General Public License   
 * along with this program; if not, write to the Free Software         
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA            
 * 02111-1307  USA                                                     
 *                                                                     
 *  All comments concerning this program package may be sent to the    
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/



///////////////////////////// COMMON LIBRARIES /////////////////////////

#include <Reconstruction/symmetries.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippImages.hh>
#include <fstream>
#include <iomanip>

//////////////////////////// SPECIFIC LIBRARIES ////////////////////////
#include "../xmippCCLattice_IO.hh"


#define GCC_VERSION (__GNUC__ * 10000 \
   + __GNUC_MINOR__ * 100 \
   + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
   #include <sstream>
#else
   #include <strstream.h>
#endif

#define VERBOSE
//#define DEBUG


// Constructor =============================================================
   CCLattice_IO::CCLattice_IO()
   {
    a.resize(2);
    b.resize(2);
  
    cc_max=0.;
    cc_peak_factor=0;
   }
// APH =====================================================================
void CCLattice_IO::read(const FileName &fn) {
  
   ifstream  fh_cor;
   int       line_no=0;
   string    line;
   float Th;
   
      
   // Open file
   fh_cor.open(fn.c_str(), ios::in);
   if (!fh_cor)
      REPORT_ERROR(1601,"UmbendFile::read: File "+fn+" not found");
   // Read fith line and skip them
   fh_cor.peek();
   for(int i=1;i<6;i++)
      getline(fh_cor,line);
   // Read Lattice Parameters   
   try {
       getline(fh_cor,line);
       dim[0] =AtoI(first_token(line));
       dim[1] =AtoI(next_token());
       O[0] =AtoI(next_token());
       O[1] =AtoI(next_token());
       a(0) =AtoF(next_token());
       a(1) =AtoF(next_token());
       getline(fh_cor,line);
       b(0) =AtoF(first_token(line));
       b(1) =AtoF(next_token());
       min_i =AtoI(next_token());
       max_i =AtoI(next_token());
       min_j =AtoI(next_token());
       max_j =AtoI(next_token());
       getline(fh_cor,line);
       cc_max=AtoF(first_token(line));
   } catch (...) {
      REPORT_ERROR(1601,"UmbendFile::read: Wrong  line in file "+fn);
   }
   //#define DEBUG
   #ifdef DEBUG
   cout << "astar vector              :   " << a  << endl
        << "bstar vector              :   " << b  << endl
        << "Cristal_dim               :   " << dim[0] << " " << dim[1] << endl
        << "Cristal_Search_Origen     :   " << O[0]   << " " <<  O[1]  << endl
        << "min_h, max_h, min_k, max_k:   " << min_i  << " " <<  max_i << " "
                                            << min_j  << " " <<  max_j << endl
        << "Maximum Correlation       :   " << cc_max  <<endl;  			       
     
   cout.flush();
   #endif
   #undef DEBUG
   
   
   
   // Read each line and keep it in the list of the Exp Lattice object
  
   float auxX,auxY,auxCCcoeff;
   //Threshold over CCvalue
   Th=cc_peak_factor* cc_max;
   
   while (!fh_cor.eof()) {
      try {
         getline(fh_cor,line);
	 if (line.length()!=0) {
	     auxX = AtoF(first_token(line));
	     auxY = AtoF(next_token());
	     auxCCcoeff = AtoF(next_token());
	   //  if(cc_max<auxCCcoeff) cc_max=auxCCcoeff;
	     if(auxX != 0 && auxY != 0 && auxCCcoeff>Th){
	     	   MRC_Xcoord.push_back(auxX);
	           MRC_Ycoord.push_back(auxY);
	           MRC_CCcoheficiente.push_back(auxCCcoeff);
	     }
	 }
      }
      catch (Xmipp_error) {
         cout << "aph File: Line " << line_no << " is skipped due to an error\n";
      }
      line_no++;
   }/* while */
   // Close file
   fh_cor.close();
   //#define DEBUG
   #ifdef DEBUG
   cout << "cor_data vector" << endl;
   for (int kk = 0; kk < MRC_Xcoord.size(); kk++)
        cout << MRC_Xcoord[kk] << " " 
	     << MRC_Ycoord[kk] << " "
	     << MRC_CCcoheficiente[kk] << endl;   
   #endif
   #undef DEBUG
}/*  UmbendFile::read */


