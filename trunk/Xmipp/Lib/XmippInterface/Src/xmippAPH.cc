/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
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

#include "../xmippAPH.hh"
#include <XmippData/xmippArgs.hh>
#include <fstream.h>
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
// APH =====================================================================
void APHFile2D::read(const FileName &fn) _THROW {
   ifstream  fh_aph;
   int       line_no=1;
   int       hmax=0, kmax=0, hmin=0, kmin=0;
   string    line;
   
   // Empties current APH File
   clear();
   
   // Open file
   fh_aph.open(fn.c_str(), ios::in);
   if (!fh_aph)
      REPORT_ERROR(1601,"aphFile::read: File "+fn+" not found");

   // Read first line and save as title
   fh_aph.peek();
   getline(fh_aph,line);
   astar.resize(2); bstar.resize(2);
   #if GCC_VERSION < 30300
      istrstream is(line.c_str());
   #else
      istringstream is(line.c_str());
   #endif
   try {
       double dummy;
       is >> dummy >> XX(astar) >> YY(astar) >> XX(bstar) >> YY(bstar)
          >> Xdim >> Ydim >> sampling_rate;
   } catch (...) {
      REPORT_ERROR(1601,"aphFile::read: Wrong first line in file "+fn);
   }

   // Read each line and keep it in the list of the aphFile object
   fh_aph.peek();
   while (!fh_aph.eof()) {
      try {
         getline(fh_aph,line);
	 int h = AtoI(first_token(line));
	 int k = AtoI(next_token());
	 hmax=MAX(hmax,h);
	 kmax=MAX(kmax,k);
	 hmin=MIN(hmin,h);
	 kmin=MIN(kmin,k);
      }
      catch (Xmipp_error) {
         cout << "aph File: Line " << line_no << " is skipped due to an error\n";
      }
      line_no++;
      fh_aph.peek();
   }/* while */
   #ifdef DEBUG   
   cout << "hmax: " << hmax <<" kmax: " << kmax <<endl;
   cout << "hmin: " << hmin <<" kmin: " << kmin <<endl;
   #endif      
     
   // Ask for memory
   spots_abs.init_zeros(kmax-kmin+1,hmax-hmin+1);
      STARTINGX(spots_abs)=hmin; STARTINGY(spots_abs)=kmin;
   spots_arg.init_zeros(kmax-kmin+1,hmax-hmin+1);
      STARTINGX(spots_arg)=hmin; STARTINGY(spots_arg)=kmin;
   IQ.init_zeros(kmax-kmin+1,hmax-hmin+1);
      STARTINGX(IQ)=hmin; STARTINGY(IQ)=kmin;
   background.init_zeros(kmax-kmin+1,hmax-hmin+1);
      STARTINGX(background)=hmin; STARTINGY(background)=kmin;
   CTF.init_zeros(background);

   // Read each line (again) and copy values to the matrices
   fh_aph.close();
   fh_aph.open(fn.c_str(), ios::in);
   line_no=1;

   // Read first line and skip it
   fh_aph.peek();
   getline(fh_aph,line);

   fh_aph.peek();
   while (!fh_aph.eof()) {
      try {
         getline(fh_aph,line);
	 int   h     	 = AtoI(first_token(line));
	 int   k     	 = AtoI(next_token());
         spots_abs(k,h)  = AtoF(next_token());
         spots_arg(k,h)  = AtoF(next_token());
	 IQ(k,h)     	 = AtoI(next_token());
	 background(k,h) = AtoF(next_token());
	 // CTF(k,h)        = AtoF(next_token());
      }
      catch (Xmipp_error XE) {
         cout << XE;
         cout << "aph File: Line " << line_no << " is skipped due to an error\n";
      }
      line_no++;
      fh_aph.peek();
   }/* while */
   // Close file
   fh_aph.close();
   
   fn_aph=fn;
}/*  APHFile2D::read */
/* ------------------------------------------------------------------------- */
void APHFile2D::clear(){
   astar.clear();
   bstar.clear();
   Xdim=Ydim=0;
   sampling_rate=0;
   spots_abs.clear();
   spots_arg.clear();
   IQ.clear();
   background.clear();
   CTF.clear();
} /*clear*/
/*----transform Xmipp Euler angles into MRC angles *--------------*/
void Euler_to_MRC(double rot, double tilt, double psi,
                  double * mrc_tilt, double * mrc_taxa)
{

//cout << "\nDEBUG rot: " << rot<<endl;
//cout << "\nDEBUG tilt: " << tilt<<endl;
//cout << "\nDEBUG psi: " << psi<<endl;

   EULER_CLIPPING_RAD(rot,tilt,psi);
   if(tilt==0) rot=0;//rot irrelevant if tilt = 0
   else if(rot<=PI/2.)   *mrc_taxa=PI/2-rot;
   else if(rot<PI*3./2.) *mrc_taxa=PI*3./2-rot;
   else if(rot<PI*2.)    *mrc_taxa=PI*5./2-rot;
   else
     {cerr <<"\nHORROR: (Euler_to_MRC) Can't find taxa\n)"; exit(1);}
   if((rot <PI/2+0.1   && rot >PI/2-.1)|| 
      (rot <PI*3/2+0.1 && rot >PI*3/2-.1)||
      (rot <0.1)|| 
      (rot >(PI*2)-.1))
      cerr <<"\nWARMING, rot close 0,90,270 or 360 degrees, conversion not reliable\n";
   

   if      ( (SGN(tilt)== +1) &&  (rot<=PI*3./2. && rot >PI/2.)) 
      { cout <<"one\n";*mrc_tilt = -tilt;}//nrg
   else if ( (SGN(tilt)== -1) &&  (rot>PI*3./2.)) 
      { cout <<"two\n"; *mrc_tilt = tilt;}//neg
   else if ( (SGN(tilt)== -1) &&  (rot<=PI/2.  )) 
      { cout <<"three\n";   *mrc_tilt = tilt;}//neg
   else if ( (SGN(tilt)== +1) &&  (rot>PI*3./2.)) 
      { cout <<"four\n"; *mrc_tilt = tilt;}//plus
   else if ( (SGN(tilt)== -1) &&  (rot<=PI*3./2. && rot >PI/2.)) 
      { cout <<"five\n"; *mrc_tilt = -tilt;}//plus
   else if ( (SGN(tilt)== +1) &&  (rot<=PI/2.  )) 
      { cout <<"six\n";   *mrc_tilt = +tilt;} //plu
   else
     {cerr <<"\nHORROR: (Euler_to_MRC) Can't find tilt\n)"; exit(1);}

//cout << "\nDEBUG rot: " << rot<<endl;
//cout << "\nDEBUG tilt: " << tilt<<endl;
//cout << "\nDEBUG psi: " << psi<<endl;
//cout << "\nDEBUG *mrc_tilt: " << *mrc_tilt<<endl;
        
}		  

