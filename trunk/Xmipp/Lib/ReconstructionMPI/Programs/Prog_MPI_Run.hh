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
#ifndef _Prog_MPI_Run_HH
#define _Prog_MPI_Run_HH

#include <fstream>
#include <iostream>
#include <XmippData/xmippFuncs.hh>
using namespace std;

/**@name MPI_Run program */
//@{
/* MPI Run Program Parameters ------------------------------------------ */
/** Parameter class for the MPI run program */
class Prog_MPI_Run_Parameters {
public:
   /** PDB file */
   FileName fn_commands;
   
   /** Number of Procesors **/
   int np;
   
   /** Machinefile name (file with working node)

public:
   /** Empty constructor */
   Prog_MPI_Run_Parameters();

   /** Read from a command line.
       An exception might be thrown by any of the internal conversions,
       this would mean that there is an error in the command line and you
       might show a usage message. */
   void read(int argc, char **argv);

   /** Usage message.
       This function shows the way of introducing this parameters. */
   void usage();

   /** Show parameters. */
   void show();
   
   /** Run. */
   void run();
};
//@}
#endif
