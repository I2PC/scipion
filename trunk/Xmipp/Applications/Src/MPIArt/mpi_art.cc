/***************************************************************************
 *
 * Authors:     José Román Bilbao Castro (jrbcast@teleline.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 1999 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#include <Reconstruction/Programs/Prog_art.hh>
#include <Reconstruction/Programs/Prog_art_crystal.hh>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <mpi.h>

/* ------------------------------------------------------------------------- */
/* Time managing stuff                                                       */
/* ------------------------------------------------------------------------- */

typedef struct {
	double user;  /* User time. */
	double sys;   /* System time. */
	double cpu;   /* CPU time = User + System. */
	double wall;  /* Wall time. */
} USWtime_t;

void uswtime(USWtime_t *tm)
{ 
	struct rusage buffer; 

	tm->wall = MPI_Wtime();
	getrusage(RUSAGE_SELF, &buffer);
	tm->user = (double) buffer.ru_utime.tv_sec + 1.0e-6 * buffer.ru_utime.tv_usec; 
	tm->sys  = (double) buffer.ru_stime.tv_sec + 1.0e-6 * buffer.ru_stime.tv_usec;  
}

/* ------------------------------------------------------------------------- */
/* Main code                                                                 */
/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[]) {
   Basic_ART_Parameters   art_prm;
   Plain_ART_Parameters   eprm;
   Crystal_ART_Parameters crystal_art_prm;
   VolumeXmipp            vol_voxels,vol_voxels_aux;
   GridVolume             vol_blobs;
   int                    crystal_mode;
   MPI_Status		  status;
   int 			  num_img_tot;		// The total amount of images 
   int			  num_img_node;		// Stores the total amount of images each node deals with
   int 			  remaining;
   int 			  Npart;
   int 			  myFirst, myLast;	// Initial and finishing indexes for each node (absolute value)

   double		  comms_t,aux_comm_t;		// Comunications time
   double		  it_t;			// iteration time
   double 		  cav_t;		// time for CAV weights calculation
   double 		  cavk_total_t,cavk_it_t;
   USWtime_t		  recons_t;		// Reconstruction time
   double		  total_t;		// Program execution time
   double 		  comms_t_it,aux_t;		// communicattions time in one iteration
	
   int 			  rank, size;	// MPI number of proccess and number of proccesses
   
	
	// Init Parallel interface		
  	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if( rank == 0 ) total_t = MPI_Wtime();		// to perform overall execution time
		
	MPI_Comm_size(MPI_COMM_WORLD, &size);
   	
	// Read Art Parameters
	try {
		art_prm.read(argc, argv);
      		// Crystal
      		crystal_mode = check_param(argc,argv,"-crystal");
      		if (crystal_mode) crystal_art_prm.read(argc, argv, art_prm);
   	} catch (Xmipp_error &XE) {
  	 	cout << XE;
  		art_prm.usage();
      		crystal_art_prm.usage();
      		exit(1);
	}
	
	/*************************** PARAMETERS INITIALIZATION ***************************/
	int num_iter = art_prm.no_it;
	cavk_total_t = 0.0;
	if( rank == 0 )
	{  				// Master code
		comms_t = 0.0;
	
      	 	art_prm.produce_Side_Info(vol_blobs,FULL,rank);
		eprm.produce_Side_Info(art_prm,vol_blobs);
	
		Basic_ART_Init_history(art_prm, eprm, vol_blobs);

		// ordered list must be the same in all nodes
		aux_comm_t = MPI_Wtime();
		MPI_Bcast( MULTIDIM_ARRAY( art_prm.ordered_list ), MULTIDIM_SIZE( art_prm.ordered_list ), MPI_INT, 0, MPI_COMM_WORLD);
		comms_t += MPI_Wtime() - aux_comm_t;
		
		cout << "Filename root: " << art_prm.fn_root << endl;
	}
	else
	{
		// each proccess can handle its own history file
		// so we add the id number to the root filename
		FileName aux = art_prm.fn_root;
		
		art_prm.fn_root = art_prm.fn_root+ItoA( rank );
		art_prm.produce_Side_Info(vol_blobs,FULL,rank);
		
		// Restore original filename.
		art_prm.fn_root = aux;
		eprm.produce_Side_Info(art_prm,vol_blobs);
		
		// ordered list must be the same in all nodes
		MPI_Bcast( MULTIDIM_ARRAY( art_prm.ordered_list ), MULTIDIM_SIZE( art_prm.ordered_list ), MPI_INT, 0, MPI_COMM_WORLD);
	}

	// How many images processes each node. It is done in such a way that every node receives the
	// same amount of them
	
	num_img_tot = art_prm.numIMG;

	Npart = (int) ceil ((float)num_img_tot / (float)size);
	remaining = num_img_tot % size;
	
	art_prm.numIMG = Npart;
	
	// each process will only see the images it is iterested in.            
	if( !remaining || rank < remaining)
        	myFirst = rank * Npart;
	else                  // for rank >= remaining > 0
        	myFirst = rank * Npart - (rank - remaining) - 1;

     	myLast = myFirst + art_prm.numIMG - 1;
	
	num_img_node = art_prm.numIMG;

	art_prm.ordered_list.window( myFirst, myLast );
	art_prm.ordered_list.startingX()=0;		
	art_prm.no_it = 1;
	
	GridVolume   vol_blobs_aux = vol_blobs;
	GridVolume   vol_aux2 = vol_blobs;
   	
	// Print some data
	if( rank == 0 ){
		if      ( art_prm.parallel_mode == Basic_ART_Parameters::SART ){ 
	             if (art_prm.block_size < size ) art_prm.block_size = size; // Assures that each processor will have at least one projection
		     cout << "SART " << "TB = " << art_prm.block_size << endl;
		}
		else if ( art_prm.eq_mode == CAV )
		     cout << "CAV" << endl;
		else if ( art_prm.parallel_mode == Basic_ART_Parameters::BiCAV ){
		     if (art_prm.block_size < size ) art_prm.block_size = size; // Assures that each processor will have at least one projection
		     cout << "BiCAV " << "TB = " << art_prm.block_size << endl;
	        }
		else if ( art_prm.parallel_mode == Basic_ART_Parameters::SIRT )
		     cout << "SIRT" << endl;
		else cout << "AVSP"<< endl;
				
		cout << "Number of processors: "<< size <<endl;
		cout << "Lambda: " << art_prm.lambda_list(0) << endl;
		cout << "Iterations: " << num_iter << endl;
		cout << "Projections: " << num_img_tot << endl;
		cout << "________________________________________________________________________________\n\n" << endl;
	}
	
	/*************************** CAV weights precalculation *************************/
	if ( art_prm.eq_mode == CAV ){
		
		// Creates and initializes special variables needed to CAV weights computation.
		GridVolumeT<int> GVNeq_aux; // This is a buffer for the communication
		
		/*
		EXTRA CALCULATIONS FOR CAV WEIGHTS: Each node computes its own part related to its images
		and after that send each others their results and sum up them. This part of code has been taken and
		modified from Basic_art.cc produce_side_info().
		*/
		cav_t = MPI_Wtime();
		art_prm.compute_CAV_weights(vol_blobs,num_img_node);
		GVNeq_aux = *(art_prm.GVNeq);
		for ( int n = 0 ; n < (art_prm.GVNeq)->VolumesNo(); n++ ){
		 	aux_comm_t = MPI_Wtime();
			MPI_Allreduce( MULTIDIM_ARRAY((*(art_prm.GVNeq))(n)()), MULTIDIM_ARRAY(GVNeq_aux(n)()), MULTIDIM_SIZE( (*(art_prm.GVNeq))(n)() ), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			comms_t += MPI_Wtime() - aux_comm_t;
			(*(art_prm.GVNeq))(n)() = GVNeq_aux(n)();	
		}
		if( rank==0 ) cout <<"Elapsed time for CAV weights computation: "<< MPI_Wtime()-cav_t<<endl;
	} 
	
	/*************************** PARALLEL ART ALGORITHM ***************************/
	
	for( int i = 0; i< num_iter; i++ )
	{
		comms_t_it = 0.0;
		cavk_it_t  = 0.0;
		
		it_t = MPI_Wtime();
		// A bit tricky. Allows us to use the sequential Basic_ART_iterations as
		// in parallel it runs internally only one iteration.
		art_prm.lambda_list(0) = art_prm.lambda(i);
		
      	        if( art_prm.parallel_mode==Basic_ART_Parameters::SART ){
		
			int blocksize;
			
			int numsteps = Npart / art_prm.block_size;
			
			// could be necessary another step for remaining projections
			if (( Npart % art_prm.block_size ) != 0) numsteps++;
			int processed = 0;
			
			art_prm.numIMG = art_prm.block_size;
			
			for(int k = 0; k < numsteps ; k++){
			
				if( k == numsteps-1 )	art_prm.numIMG = num_img_node - processed;
				
				art_prm.ordered_list.startingX( ) = -processed; 
				
				processed += art_prm.numIMG;
				
				for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
					vol_aux2(j)() = vol_blobs(j)();
				
				Basic_ART_iterations(art_prm, eprm, vol_blobs, rank);
			
				blocksize = art_prm.numIMG * size;
				
				for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
				{
					vol_blobs(j)() = vol_blobs(j)()-vol_aux2(j)(); // Adapt result to parallel ennvironment from sequential routine
					vol_blobs(j)() *= ((double)art_prm.numIMG /(double) blocksize);
					aux_comm_t = MPI_Wtime();
					MPI_Allreduce( MULTIDIM_ARRAY(vol_blobs(j)()), MULTIDIM_ARRAY(vol_blobs_aux(j)()), MULTIDIM_SIZE( vol_blobs(j)()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					aux_t = MPI_Wtime() - aux_comm_t;
					comms_t += aux_t;
					comms_t_it += aux_t;
					vol_blobs(j)() = vol_aux2(j)() + vol_blobs_aux(j)();
				}	
			}
	        }
		else if( art_prm.parallel_mode==Basic_ART_Parameters::BiCAV ){
			// Creates and initializes special variables needed to CAVK weights computation.
			GridVolumeT<int> GVNeq_aux; // This is a buffer for communications
			
			int numsteps = Npart / art_prm.block_size;
	
			if(( Npart % art_prm.block_size) != 0) numsteps++;
			int processed = 0;
				
			art_prm.numIMG = art_prm.block_size;
			
			art_prm.eq_mode = CAV; // Another trick
			for(int ns = 0; ns < numsteps ; ns++){
				
				if( ns == numsteps-1 )	art_prm.numIMG = num_img_node - processed;
				
				art_prm.ordered_list.startingX()= -processed;	
				
			        /*
				EXTRA CALCULATIONS FOR CAVK WEIGHTS: Each node computes its own part related to its images
				and after that send each others their results and sum up them. This part of code has been taken and
				modified from Basic_art.cc produce_side_info().
				*/
				
				cav_t = MPI_Wtime();
                                art_prm.compute_CAV_weights(vol_blobs,art_prm.numIMG);
				GVNeq_aux = *(art_prm.GVNeq);
				 
				for ( int n = 0 ; n < (art_prm.GVNeq)->VolumesNo(); n++ )
				{
					aux_comm_t = MPI_Wtime();
					MPI_Allreduce( MULTIDIM_ARRAY((*(art_prm.GVNeq))(n)()), MULTIDIM_ARRAY(GVNeq_aux(n)()), MULTIDIM_SIZE( (*(art_prm.GVNeq))(n)() ), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
					aux_t = MPI_Wtime() - aux_comm_t;
					comms_t += aux_t;
					comms_t_it += aux_t;
					(*(art_prm.GVNeq))(n)() = GVNeq_aux(n)();
				}	
				cavk_it_t += MPI_Wtime()-cav_t;
				
				for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
					vol_aux2(j)() = vol_blobs(j)();
					
				Basic_ART_iterations(art_prm, eprm, vol_blobs, rank);
				
				processed += art_prm.numIMG;
				
				for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
				{
					vol_blobs(j)() = vol_blobs(j)()-vol_aux2(j)(); // Adapt result to parallel ennvironment from sequential routine
					aux_comm_t = MPI_Wtime();
					MPI_Allreduce( MULTIDIM_ARRAY(vol_blobs(j)()), MULTIDIM_ARRAY(vol_blobs_aux(j)()), MULTIDIM_SIZE( vol_blobs(j)()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
					aux_t = MPI_Wtime() - aux_comm_t;
					comms_t += aux_t;
					comms_t_it += aux_t;
					vol_blobs(j)() = vol_aux2(j)() + vol_blobs_aux(j)();
				}	
			}
			art_prm.eq_mode = CAVK; // trick undone
		}
		else if( art_prm.eq_mode == CAV){
		
			//CAV weights calculations are done before iterations begin in order to avoid recalculate them
			for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
				vol_aux2(j)() = vol_blobs(j)();
				
			Basic_ART_iterations(art_prm, eprm, vol_blobs, rank);
				
			for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
			{
				vol_blobs(j)() = vol_blobs(j)()-vol_aux2(j)(); // Adapt result to parallel ennvironment from sequential routine
				aux_comm_t = MPI_Wtime();
				MPI_Allreduce( MULTIDIM_ARRAY(vol_blobs(j)()), MULTIDIM_ARRAY(vol_blobs_aux(j)()), MULTIDIM_SIZE( vol_blobs(j)()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				aux_t = MPI_Wtime() - aux_comm_t;
				comms_t += aux_t;
				comms_t_it += aux_t;
				vol_blobs(j)() = vol_aux2(j)() + vol_blobs_aux(j)(); 
			}	
		}
		else{   // SIRT AND AVSP
			if(art_prm.SIRT)
				for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
					vol_aux2(j)() = vol_blobs(j)();
			
			Basic_ART_iterations(art_prm, eprm, vol_blobs, rank);
		
			for ( int j = 0 ; j < vol_blobs.VolumesNo() ; j++)
			{
				if(art_prm.SIRT)
				{
					vol_blobs(j)() = vol_blobs(j)()-vol_aux2(j)(); // Adapt result to parallel ennvironment from sequential routine
					vol_blobs(j)() *= ((double)art_prm.numIMG /(double) num_img_tot);
				}
				aux_comm_t = MPI_Wtime();
				MPI_Allreduce( MULTIDIM_ARRAY(vol_blobs(j)()), MULTIDIM_ARRAY(vol_blobs_aux(j)()), MULTIDIM_SIZE( vol_blobs(j)()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				aux_t = MPI_Wtime() - aux_comm_t;
				comms_t += aux_t;
				comms_t_it += aux_t;
				if( !art_prm.SIRT )
				{
					vol_blobs(j)() = vol_blobs_aux(j)();
					vol_blobs(j)() /= size; // Non-SIRT Normalization
				}
				else vol_blobs(j)() = vol_aux2 (j)() + vol_blobs_aux(j)();
			}
		}
		cavk_total_t += cavk_it_t;
		if ( rank == 0 ){
			cout << "\nIteration " << i << endl;
			cout << "Time: " << MPI_Wtime() - it_t << " secs." << endl; 
			cout << "Comms. time: " << comms_t_it << " secs." << endl;
			if(  art_prm.parallel_mode==Basic_ART_Parameters::BiCAV  )
				cout << "BiCAV weighting time: " << cavk_it_t << endl;
		
			if( i < num_iter-1 ){
				int Xoutput_volume_size=(art_prm.Xoutput_volume_size==0) ?
				    art_prm.projXdim:art_prm.Xoutput_volume_size;
				int Youtput_volume_size=(art_prm.Youtput_volume_size==0) ?
				    art_prm.projYdim:art_prm.Youtput_volume_size;
 				int Zoutput_volume_size=(art_prm.Zoutput_volume_size==0) ?
				    art_prm.projXdim:art_prm.Zoutput_volume_size;
		 		blobs2voxels(vol_blobs, art_prm.blob, &vol_voxels, art_prm.D, 
        			Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
        			vol_voxels.write(art_prm.fn_root+"it"+ItoA(i+1)+".spi");
			}
			if ((art_prm.tell&TELL_SAVE_BLOBS) && (i < num_iter-1 ))
 			{
				vol_blobs.write(art_prm.fn_root+"it"+ItoA(i+1)+".blob");
			}
		}
	}	

	/*************************** FINISHING AND STORING VALUES ***************************/
		
	if( rank > 0 ){
		art_prm.fh_hist.close();			
		return 0;
	}
		
 	int Xoutput_volume_size=(art_prm.Xoutput_volume_size==0) ?
	  art_prm.projXdim:art_prm.Xoutput_volume_size;
	int Youtput_volume_size=(art_prm.Youtput_volume_size==0) ?
	  art_prm.projYdim:art_prm.Youtput_volume_size;
 	int Zoutput_volume_size=(art_prm.Zoutput_volume_size==0) ?
	  art_prm.projXdim:art_prm.Zoutput_volume_size;
			 		
	blobs2voxels(vol_blobs, art_prm.blob, &vol_voxels, art_prm.D,Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
 		
	vol_voxels.write(art_prm.fn_root+".spi");
		
	if (art_prm.tell&TELL_SAVE_BLOBS) 
	vol_blobs.write(art_prm.fn_root+".blob");

	art_prm.fh_hist.close();			
   	uswtime( &recons_t );
	cout << "\n\n------ FINAL STATISTICS ------" << endl;
	cout << "\nTOTAL EXECUTION TIME: " << recons_t.wall - total_t << endl;
	cout << "Communications time: " << comms_t << " secs." << endl;
	cout << "CPU time: " << recons_t.user + recons_t.sys<< " secs." << endl;
	cout << "USER: " << recons_t.user << " SYSTEM: " << recons_t.sys << "\n\n" << endl;  
	if(  art_prm.parallel_mode==Basic_ART_Parameters::BiCAV )
		cout << "total BiCAV Weighting time: "<< cavk_total_t << endl;								
        MPI_Finalize();	
        return 0 ;
}






