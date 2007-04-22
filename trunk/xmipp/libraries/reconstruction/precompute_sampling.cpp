/***************************************************************************
 *
 * Authors:
 *
 * Roberto Marabini 
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

#include "precompute_sampling.h"

/* Empty constructor ------------------------------------------------------- */
Prog_Sampling_Parameters::Prog_Sampling_Parameters() {
   /** sampling object 1 by default*/
   mysampling.SetSampling(1);
}


/* Read parameters --------------------------------------------------------- */
void Prog_Sampling_Parameters::read(int argc, char **argv) {
   int i=0;
#ifdef NEVERDEFINED
   fn_pdb=get_param(argc,argv,"-i");
   fn_out=get_param(argc,argv,"-o","");
   if (fn_out=="") fn_out=fn_pdb.without_extension();
   Ts=AtoF(get_param(argc,argv,"-sampling_rate","1"));
   highTs=AtoF(get_param(argc,argv,"-high_sampling_rate","0.1"));
   output_dim=AtoI(get_param(argc,argv,"-output_dim","-1"));
   #endif
}

/* Usage ------------------------------------------------------------------- */
void Prog_Sampling_Parameters::usage() {
   cerr << "PDBphantom\n"
        << "   -i <pdb file>                    : File to process\n"
        << "  [-o <fn_root>]                    : Root name for output\n"
        << "  [-sampling_rate <Ts=1>]           : Sampling rate (Angstroms/pixel)\n"
        << "  [-high_sampling_rate <highTs=0.1>]: Sampling rate before downsampling\n"
        << "  [-size <output_dim>]              : Final size in pixels (must be a power of 2)\n"
        << "\n"
        << "Example of use: Sample at 1.6A and limit the frequency to 10A\n"
        << "   xmipp_pdbphantom -i 1o7d.pdb -sampling_rate 1.6\n"
        << "   xmipp_fourierfilter -i 1o7d.vol -o 1o7d_filtered.vol -low_pass 10 -sampling 1.6 -fourier_mask raised_cosine 0.1\n"
   ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Sampling_Parameters::show() {
   cout
        << "Sampling rate:      " << sampling    << endl
        << "symmetry file name: " << symmetry_file << endl
   ;
}

/* Run --------------------------------------------------------------------- */
void Prog_Sampling_Parameters::run() {
   sampling=4;
   symmetry_file="ico.sym";
   mysampling.SetSampling(sampling);
   mysampling.Compute_sampling_points(false);
   cerr << " read_sym_file " << endl;
   SL.read_sym_file(symmetry_file);
   //fill vector with symmetry axis
   cerr << " remove_redundant_points " << endl;
   mysampling.remove_redundant_points(SL);
   #define DEBUG6
   #ifdef DEBUG6
       for (int i=0; i<mysampling.no_redundant_sampling_points_vector.size(); i++) {
           cout << mysampling.no_redundant_sampling_points_vector[i].transpose() << " 1 1 " << endl;
       }
   #endif
   #undef DEBUG6
}

