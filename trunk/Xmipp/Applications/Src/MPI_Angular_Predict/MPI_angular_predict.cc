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

#include <Reconstruction/Programs/Prog_Angular_Predict.hh>
#include <XmippData/xmippImages.hh>
#include <mpi.h>

int main (int argc, char **argv) {
SelFile ImagesFile;
FileName fn_sym;
int MyFirst, MyLast;	// First and last file to process each node
Prog_angular_predict_prm prm;
double shiftX, shiftY, psi, rot, tilt;
int remaining;
int Npart;
int rank, size, Img_Nbr, myFirst, myLast;
MPI_Status  status;            	// Stores MPI directives status
double buff[8];
   		   	
	MPI_Init(&argc, &argv);
   	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
   	
	prm.produce_side_info();
	ImagesFile.read( prm.fn_in );
	Img_Nbr = prm.get_images_to_process();
	
//	Npart = (int) ceil ((float)Img_Nbr / (float)size);
	Npart = (int) ((float)Img_Nbr / (float)size);

	remaining = Img_Nbr % size;
	
	if( !remaining || rank < remaining)
        	myFirst = rank * Npart;
	else                  // for rank >= remaining > 0
        	myFirst = rank * Npart - (rank - remaining) - 1;

     	myLast = myFirst + Npart - 1;
	
	int p = (prm.predicted_rot).size();	
   	 
	if(rank == 0)
	{	
  		matrix1D<double> v(7);
 		DocFile DF;
   		DF.reserve(p+1);
   		DF.append_comment("Predicted_Rot Predicted_Tilt Predicted_Psi Predicted_ShiftX Predicted_ShiftY Corr");		
		
		int to_save = Img_Nbr;
		while( to_save > 0)
		{
			MPI_Recv( buff, 8, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			   	
      			v(0) = buff[ 0 ];
      			v(1) = buff[ 1 ];
      			v(2) = buff[ 2 ];
      			v(3) = buff[ 3 ];
      			v(4) = buff[ 4 ];
      			v(5) = buff[ 5 ];
		
			// position in the file
			int position = (int)buff[ 6 ];
			DF.go_first_data_line();
			DF.jump( position );
	      		DF.insert_data_line(v);	
			to_save --;	
		}	
		DF.write(prm.fn_out_ang);
	}
	else
	{
		ImageXmipp image;
		ImagesFile.go_beginning();
		ImagesFile.jump(MyFirst); 
		for( int i = MyFirst; i <= MyLast; i++)
		{
			image.read(ImagesFile.get_current_file());
			double corr = prm.predict_angles( image, shiftX, shiftY, rot, tilt, psi);
			buff[0] = prm.predicted_rot[i];
      			buff[1] = prm.predicted_tilt[i];
      			buff[2] = prm.predicted_psi[i];
      			buff[3] = prm.predicted_shiftX[i];
      			buff[4] = prm.predicted_shiftY[i];
      			buff[5] = prm.predicted_corr[i];
			buff[6] = i;
			MPI_Send(buff, 8, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
   			ImagesFile.next();
		}
	}
	
        MPI_Finalize();	
        return 0 ;
}
