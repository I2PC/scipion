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
   sampling_file_root=get_param(argc,argv,"-o");
   symmetry=get_param(argc,argv,"-symmetry","cn");
   sym_order=AtoI(get_param(argc,argv,"-sym_order","1"));
   sampling=AtoF(get_param(argc,argv,"-sampling_rate","5"));
   neighborhood=AtoF(get_param(argc,argv,"-neighborhood","1"));
}

/* Usage ------------------------------------------------------------------- */
void Prog_Sampling_Parameters::usage() {
   cerr << "precompute_sampling\n"
        << "   -o root_file_name           : Root for output files\n"
        << "  [-symmetry cn]   :One of the 17 possible symmetries in\n" 
        << "                                single particle electronmicroscopy\n"
        << "                                i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv, dnh, t, td, th, o, oh, i, ih\n"
        << "  [-sym_order 1]               : For infinite groups symmetry order\n"
        << "  [-sampling_rate 5]           : Distance in degrees between sampling points\n"
        << "  [-neighborhood 1]            : A sampling point is neighbor if closer than this value in degrees\n"
        << "\n"
        << "Example of use: Sample at 2degres and compute neighboor at "
        << " 5 degrees for c6 symmetry\n"
        << "   xmipp_precompute_sampling -o out -symmetry c6 "
        << " -sampling_rate 2 -neighborhood 5\n"
   ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Sampling_Parameters::show() {
   cout
        << "Sampling rate:      " << sampling    << endl
        << "output files root:  " << sampling_file_root << endl
        << "symmetry group:     " << symmetry << endl
        << "symmetry order:     " << sym_order << endl
        << "neighborhood:       " << neighborhood << endl
  ;
}



/* Run --------------------------------------------------------------------- */
void Prog_Sampling_Parameters::run() {
   show();
   mysampling.SetSampling(sampling);
   mysampling.Compute_sampling_points(false);
   mysampling.create_sym_file(symmetry,sym_order);
   mysampling.remove_redundant_points(symmetry,sym_order);
   #define DEBUG6
   #ifdef DEBUG6
       for (int i=0; i<mysampling.no_redundant_sampling_points_vector.size(); i++){  
           cout << mysampling.no_redundant_sampling_points_vector[i].transpose() << " 1.1 2 " << endl;
           //cout << mysampling.no_redundant_sampling_points_angles[i].transpose() << " 1.21 1 " << endl;
           }
       for (int i = 0; 
            i < mysampling.sampling_points_vector.size(); 
            i++)
          cout  <<  mysampling.sampling_points_vector[i].transpose()  << " 1 1 " << endl;
   #endif
   #undef DEBUG6
}

