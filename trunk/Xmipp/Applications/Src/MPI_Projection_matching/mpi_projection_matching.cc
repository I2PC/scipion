/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)   
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
#include <Reconstruction/Programs/Prog_projection_matching.hh> 
#include <mpi.h>

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv) {

  // For parallelization
  int num_img_tot, num_img_node;
  int myFirst, myLast, remaining, Npart;
  int rank, size;

  double                        aux,sumCC,sumZ;
  FileName                      fn_img,fn_tmp;
  DocFile                       DFo;
  Prog_projection_matching_prm  prm;

  // Init Parallel interface		
  MPI_Init(&argc, &argv);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Get input parameters
  try {
    // Read command line & produce side info
    prm.read(argc,argv);
    if (rank==0) prm.show();
    else  { prm.verb=0; prm.output_refs=false; }

    prm.produce_Side_info();

    // Calculate indices myFirst and myLast and adapt prm.SF
    prm.SF.clean_comments();
    prm.SF.clean();
    num_img_tot = prm.SF.ImgNo();
    Npart = (int) ceil ((float)num_img_tot / (float)size);
    remaining = num_img_tot % size;
    if ( rank < remaining ) {
      myFirst = rank * (Npart + 1);
      myLast = myFirst + Npart;
    } else {
      myFirst = rank * Npart + remaining;
      myLast = myFirst + Npart - 1;
    }
    // Now discard all images in Selfile that are outside myFirst-myLast
    prm.SF.go_beginning();
    SelFile SFpart;
    SFpart.clear();
    for (int nr=0; nr<num_img_tot; nr++) {
      if ((nr>=myFirst) && (nr<=myLast)) {
	prm.SF.go_beginning();
	prm.SF.jump_lines(nr);
	SFpart.insert(prm.SF.current());
     }
    }
    prm.SF=SFpart;

  } catch (Xmipp_error XE) {if (rank==0) {cout << XE; prm.usage();} MPI_Finalize(); exit(1);}
    
  try {

    DFo.clear();
    if (rank==0) 
      DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), maxCC (6), Z-score (7)");

    // Process all images
    prm.PM_loop_over_all_images(prm.SF,DFo,sumCC,sumZ);
 
    // Here MPI_allreduce 
    MPI_Allreduce(&sumCC,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    sumCC=aux;
    MPI_Allreduce(&sumZ,&aux,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    sumZ=aux;
    if (prm.verb>0) cerr << " Average maxCC = "<<sumCC/num_img_tot<<" average Z-score = "<<sumZ/num_img_tot <<endl;


     // All nodes write out temporary DFo
    fn_img.compose(prm.fn_out,rank,"tmpdoc");
    DFo.write(fn_img);

    MPI_Barrier(MPI_COMM_WORLD);

    // Master writes out final docfile
    if (rank==0) {
      DFo.clear();
      for (int rank2=0; rank2<size; rank2++) {
	fn_img.compose(prm.fn_out,rank2,"tmpdoc");
	int ln=DFo.LineNo();
	DFo.append(fn_img);
	DFo.locate(DFo.get_last_key());
	DFo.next();
	DFo.remove_current();
	system(((string)"rm -f "+fn_img).c_str());
      }
      DFo.write(prm.fn_out);
    }


  } catch (Xmipp_error XE) {if (rank==0) {cout << XE; prm.usage();} MPI_Finalize(); exit(1);}

  MPI_Finalize();	
  return 0;

}




