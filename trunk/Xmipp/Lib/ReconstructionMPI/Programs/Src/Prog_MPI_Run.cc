/***************************************************************************
 *
 * Authors:
 *
 * Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include "../Prog_MPI_Run.hh"
#include <XmippData/xmippArgs.hh>

/* Empty constructor ------------------------------------------------------- */
Prog_MPI_Run_Parameters::Prog_MPI_Run_Parameters() {
 int i=0;
}


/* Read parameters --------------------------------------------------------- */
void Prog_MPI_Run_Parameters::read(int argc, char **argv) {
   fn_commands=get_param(argc,argv,"-i");
}

/* Usage ------------------------------------------------------------------- */
void Prog_MPI_Run_Parameters::usage() {
   cerr << "MPI_Run\n"
        << "   -i <command file>    : File with commands to send to mpirun\n"
        << "\n"
        << "Example of use:\n"
        << "   xmipp_mpi_run -i commandd_file\n"
   ;
}

/* Show -------------------------------------------------------------------- */
void Prog_MPI_Run_Parameters::show() {
   cout << "Commands  file:           " << fn_commands << endl
   ;
}


/* Run --------------------------------------------------------------------- */
void Prog_MPI_Run_Parameters::run() {
   ifstream fh_in;
   fh_in.open(fn_commands.c_str());
   if (!fh_in)
      REPORT_ERROR(1,(string)"Cannot open "+fn_commands);
   
   while (!fh_in.eof()) {
      string line;
      getline(fh_in,line);
      cout << line << endl;
   }
   
   fh_in.close();
}
