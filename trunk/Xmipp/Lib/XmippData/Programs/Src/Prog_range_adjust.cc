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

#include "../Prog_range_adjust.hh"
#include <XmippData/xmippArgs.hh>

/* Read parameters --------------------------------------------------------- */
void Prog_Range_adjust_Parameters::read(int argc, char **argv) _THROW {
   min_sigma=AtoF(get_param(argc,argv,"-min"));
   max_sigma=AtoF(get_param(argc,argv,"-max"));
   randomize_random_generator();
}

/* Usage ------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::usage() {
   cerr << "   -min <min_sigma %> : Variation of the minimum value\n"
        << "   -max <max_sigma %> : Variation of the maximum value\n"
   ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::show() {
   cout << "Minimum noise: " << min_sigma << endl
        << "Maximum noise: " << max_sigma << endl;
}

/* Apply ------------------------------------------------------------------- */
void Prog_Range_adjust_Parameters::apply(Image *I) {
   double amin=rnd_gaus(1,min_sigma);
   double amax=rnd_gaus(1,max_sigma);
   double minval, maxval;
   (*I)().compute_double_minmax(minval,maxval);
   minval*=amin;
   maxval*=amax;
   (*I)().range_adjust(minval,maxval);
}

void Prog_Range_adjust_Parameters::apply(Volume *V) {
   double amin=rnd_gaus(1,min_sigma);
   double amax=rnd_gaus(1,max_sigma);
   double minval, maxval;
   (*V)().compute_double_minmax(minval,maxval);
   minval*=amin;
   maxval*=amax;
   (*V)().range_adjust(minval,maxval);
}
