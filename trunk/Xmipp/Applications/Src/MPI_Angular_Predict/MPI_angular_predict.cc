/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao Castro (coss@cnb.uam.es)
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

#include <Reconstruction/Programs/Prog_Angular_Predict.hh>
#include <mpi.h>

int main (int argc, char **argv) {
SF ImagesFile;
FileName fn_sym;
int MyFirst, MyLast;	// First and last file to process each node
Prog_angular_predict_prm prm;
double shiftX, shiftY, psi, rot, tilt;
int remaining;
int NPart;
   	
	MPI_Init(&argc, &argv);
   	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
   	
	prm.produce_side_info();
	ImagesFile.read( prm.fn_in );
	Img_Nbr = prm.get_images_to_process();
	
	Npart = (int) ceil ((float)Img_Nbr / (float)size);
	remaining = num_img_tot % size;
	
	if( !remaining || rank < remaining)
        	myFirst = rank * Npart;
	else                  // for rank >= remaining > 0
        	myFirst = rank * Npart - (rank - remaining) - 1;

     	myLast = myFirst + Npart - 1;
	 
	for( int i = MyFirst; i <= MyLast; i++)
	{
		XmippImage image;
		image.read(ImagesFile[i].fn_proj);
     		double corr = prm->predict_angles( image, shiftX, shiftY, rot, tilt, psi);
	} 
	
	if( rank == 0){
		// Receives results from others and stores values in a file.
		waiting_for = Img_Nbr - MyLast; 
		MPI_DOUBLE buff[8];
		while( waiting_for > 0){
			MPI_Recv( buff, 8, MPI_DOUBLE, MPI_ANY, MPI_ANY, MPI_WORLD_COMM, &status );
			waiting_for--;
		}
	   	// Save predicted angles
   		int p = (prm.predicted_rot).size();
   		DocFile DF;
   		DF.reserve(p+1);
   		DF.append_comment("Predicted_Rot Predicted_Tilt Predicted_Psi Predicted_ShiftX Predicted_ShiftY Corr");
   		matrix1D<double> v(7);
   		for (int i=0; i<p; i++) {
      			v(0) = prm.predicted_rot[i];
      			v(1) = prm.predicted_tilt[i];
      			v(2) = prm.predicted_psi[i];
      			v(3) = prm.predicted_shiftX[i];
      			v(4) = prm.predicted_shiftY[i];
      			v(5) = prm.predicted_corr[i];
      			DF.append_data_line(v);
   		}
	   	DF.write(prm.fn_out_ang);
	}
   	else{
	   	// Sends predicted angles
   		MPI_DOUBLE buff[8];
   		for (int i=myFirst; i<=MyLast; i++) {
      			buff[0] = prm.predicted_rot[i];
      			buff[1] = prm.predicted_tilt[i];
      			buff[2] = prm.predicted_psi[i];
      			buff[3] = prm.predicted_shiftX[i];
      			buff[4] = prm.predicted_shiftY[i];
      			buff[5] = prm.predicted_corr[i];
			buff[6] = i;
			MPI_Send(buff, 8, MPI_DOUBLE, 0, 0, MPI_WORLD_COMM);
   		}
	}
	
        MPI_Finalize();	
        return 0 ;
}
