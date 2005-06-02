/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao Castro (jrbcast@cnb.uam.es)
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
/* INCLUDES ---------------------------------------------------------------- */
#include <Reconstruction/Programs/Prog_WBP.hh>
#include <mpi.h>

int main (int argc, char **argv) {

  // For parallelization
  int           rank, size,num_img_tot;
  // For program
  VolumeXmipp   vol,aux;
  Prog_WBP_prm  prm;

  // Init Parallel interface		
  MPI_Init(&argc, &argv);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  try {
    
    // Read command line & produce side info
    prm.read(argc,argv);
    if (rank==0) prm.show();
    else  { prm.verb=0 }

    // First produce the filter and then cut selfile in smaller parts
    prm.produce_Side_info();

    // Select only relevant part of selfile for this rank
    prm.SF.mpi_select_part(rank,size,num_img_tot);

    // Actual backprojection
    prm.apply_2Dfilter_arbitrary_geometry(prm.SF, vol);

    aux().resize(vol());
    MPI_Allreduce(MULTIDIM_ARRAY(vol()),MULTIDIM_ARRAY(aux()),
		  MULTIDIM_SIZE(vol()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    vol=aux;
    
    if (rank==0) vol.write(prm.fn_out);

  } catch (Xmipp_error XE) {if (rank==0) {cout << XE; prm.usage();} MPI_Finalize(); exit(1);}

  MPI_Finalize();	
  return 0;

}

