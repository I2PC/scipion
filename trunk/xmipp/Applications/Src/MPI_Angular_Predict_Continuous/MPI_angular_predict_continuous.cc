/***************************************************************************
 *
 * Authors:     Carlos Óscar Sánchez Sorzano (coss.eps@ceu.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 * Lab. de Bioingeniería, Univ. San Pablo CEU
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

#include <Reconstruction/Programs/Prog_Angular_Predict_continuous.hh>
#include <mpi.h>

int main (int argc, char **argv) {
   // Initialize MPI
   int rank, NProcessors;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &NProcessors);

   Prog_angular_predict_continuous_prm prm;
   try {
      // Read input parameters
      prm.read(argc,argv);
   } catch (Xmipp_error XE) {
      if (rank==0) {
	 cout << XE;
	 prm.usage();
      }
      exit(1);
   }

   try {
      // Prepare side information
      prm.produce_side_info();

      // Divide the selfile in chunks
      int imgNbr = prm.get_images_to_process();
      int Nchunk = (int) ((float)imgNbr/(float)(NProcessors-1));
      int myFirst=(rank-1)*Nchunk;
      int myLast=rank*Nchunk-1;
      if (rank==NProcessors-1) myLast=imgNbr-1;
      
      // Make the alignment, rank=0 receives all the assignments
      // The rest of the ranks compute the angular parameters for their
      // assigned images
      double v[7];
      if (rank==0) {	
	 int toGo=imgNbr;
	 MPI_Status status;
	 while( toGo > 0) {
	    MPI_Recv(v, 7, MPI_DOUBLE, MPI_ANY_SOURCE,
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    int i=(int)v[0];
	    prm.predicted_rot[i]=v[1];
	    prm.predicted_tilt[i]=v[2];
	    prm.predicted_psi[i]=v[3];
	    prm.predicted_shiftX[i]=v[4];
	    prm.predicted_shiftY[i]=v[5];
	    prm.predicted_cost[i]=v[6];
	    toGo--;	
	 }
      } else {
         SelFile SF_in(prm.fn_in);
	 SF_in.go_beginning();
	 SF_in.jump(myFirst);
	 init_progress_bar(myLast-myFirst+1);
	 for (int i=myFirst; i<=myLast; i++) {
	    // Read image and estimate angular parameters
            ImageXmipp I;
	    I.read(SF_in.NextImg(),false,false,false);
	    I().set_Xmipp_origin();
	    double shiftX, shiftY, psi, rot, tilt;
	    prm.predict_angles(I, shiftX, shiftY, rot, tilt, psi);
	    double cost=prm.predicted_cost[i-myFirst];
	    
	    // Rewrite the image header if necessary
	    if (!prm.dont_modify_header) {
	       I.read(I.name());
	       I.set_eulerAngles(rot,tilt,psi);
	       I.set_originOffsets(shiftX,shiftY);
	       I.write();
	    }

      	    // Send the alignment parameters to the master
	    v[0]=i;
	    v[1]=rot;
	    v[2]=tilt;
	    v[3]=psi;
	    v[4]=shiftX;
	    v[5]=shiftY;
	    v[6]=cost;
	    MPI_Send(v, 7, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	    progress_bar(i-myFirst+1);
	 }
      }
      progress_bar(myLast-myFirst+1);

      MPI_Finalize();	
      if (rank==0) prm.finish_processing();
      return 0 ;
   } catch (Xmipp_error XE) {
      cout << XE << endl;
   }
}
