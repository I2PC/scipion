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
#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFFT.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippProjection.hh>
#include <Reconstruction/projection.hh>
#include <mpi.h>

int main (int argc, char **argv) 
{

   ImageXmipp img,orig_img,aux;
   VolumeXmipp vol;
   FourierImageXmipp IMG;
   FourierVolumeXmipp VOL;
   SelFile SF_in;
   FileName fn_sel, fn_out, fn_in;
   Projection curr_proj;
   int rank, size; // MPI identifiers
   int Npart, remaining, myFirst, myLast;
   MPI_Status		  status;            	// Stores MPI directives status

   // Read command line .................................................... 

   // Only the master process stores the result

   // Init Parallel interface		
   MPI_Init(&argc, &argv);

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if( rank == 0 ) fn_out = get_param(argc,argv,"-o"); 

   fn_sel =  get_param(argc,argv,"-i"); 
	
   if (!exists( fn_sel ))
   {
	MPI_Finalize();
	return 0;
   }
    
   SF_in.read( fn_sel );

   int xdim, ydim; 
   SF_in.ImgSize(xdim, ydim); 
   int Imagedim = max( xdim, ydim );
   VolumeXmipp Vol(Imagedim, Imagedim, Imagedim);
   Vol().set_Xmipp_origin();

   FileName fn_proj;
   // Discard unuseful images  
   SF_in.clean();

   Npart = (int) ceil ((float)SF_in.ImgNo() / (float)size);
   remaining = SF_in.ImgNo() % size;
	
   // each process will only see the images it is iterested in.            
   if( !remaining || rank < remaining)
       	myFirst = rank * Npart;
   else                  // for rank >= remaining > 0
       	myFirst = rank * Npart - (rank - remaining) - 1;

   myLast = myFirst + Npart - 1;
	
   if ( rank == 0 )
   {
	VolumeXmipp Volbuff(Imagedim, Imagedim, Imagedim);
	Volbuff().set_Xmipp_origin();

	SF_in.go_first_ACTIVE();
        for( int i=0; i<Npart; i++ )
	{
		fn_proj = SF_in.get_current_file();
		SF_in.next();		
		curr_proj.read( fn_proj );
		curr_proj().set_Xmipp_origin();

		FourierTransform( curr_proj(), IMG());
		int Nx = IMG().ColNo();
     		int Ny = IMG().RowNo();

		FOR_ALL_ELEMENTS_IN_MATRIX2D(IMG())
		{ 
		     int i2 = i+Ny/2;
		     int j2 = j+Nx/2;
	       	     if( i2 > Ny/2 ) i2 -= Ny;
	             if( j2 > Nx/2 ) j2 -= Nx;
	
	             double sqcomp = sqrt((double)(i2*i2)+(double)(j2*j2));
	             IMG( i,j ) *= sqcomp;
                }
                InverseFourierTransform(IMG(), curr_proj());
     
                singleWBP( Vol(), curr_proj );
	}
	for(int i=0;i<size-1;i++)
	{
		MPI_Recv( MULTIDIM_ARRAY(Volbuff()), 
			MULTIDIM_SIZE(Volbuff()), 
			MPI_DOUBLE, 
			MPI_ANY_SOURCE,
			MPI_ANY_TAG, 
			MPI_COMM_WORLD, 
			&status );
		Vol()+=Volbuff();	
	}	
   }
   else
   {
   	SF_in.go_first_ACTIVE();
   	SF_in.jump( myFirst );
   	for( int i=0; i<Npart; i++ )
	{
     		fn_proj = SF_in.get_current_file();
		SF_in.next();
		curr_proj.read( fn_proj );
     
     		curr_proj().set_Xmipp_origin();
     		FourierTransform(curr_proj(), IMG()); 
     		IMG().set_Xmipp_origin();
     
     		int Nx = IMG().ColNo();
     		int Ny = IMG().RowNo();

		FOR_ALL_ELEMENTS_IN_MATRIX2D(IMG())
		{ 
		    int i2 = i+Ny/2;
	   	    int j2 = j+Nx/2;
 		    if( i2 > Ny/2 ) i2 -= Ny;
		    if( j2 > Nx/2 ) j2 -= Nx;
	
		    double sqcomp = sqrt((double)(i2*i2)+(double)(j2*j2));
		    IMG( i,j ) *= sqcomp;
	        }
		InverseFourierTransform(IMG(), curr_proj());
     
	        singleWBP( Vol(), curr_proj );
        }
 	MPI_Send(MULTIDIM_ARRAY(Vol()), 
		MULTIDIM_SIZE(Vol()), 
		MPI_DOUBLE, 
		0, 
		0, 
		MPI_COMM_WORLD );
    }
    
    if( rank > 0 )
    {
	MPI_Finalize();	  // Must exist for each proccess on MPI evironment
       	return 0;
    }
    
    Vol() = Vol()/SF_in.ImgNo();
    Vol.write( fn_out );
    MPI_Finalize();
    return 0;
}

 
