/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "mpi_angular_projection_matching.h"

/*Some constast to message passing tags */
#define TAG_START_JOB	 	1
#define TAG_STOP_JOB 		2

/*Constructor */
MpiProgAngularProjectionMatching::MpiProgAngularProjectionMatching()
{
    imagesBuffer = NULL;
    last_chunk = NULL;
}
/* Destructor */
MpiProgAngularProjectionMatching::~MpiProgAngularProjectionMatching()
{
    if (node->isMaster())
    {
        delete[] imagesBuffer;
        delete[] last_chunk;
    }
    delete node; //this calls MPI_Finalize
}

void MpiProgAngularProjectionMatching::read(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
    // Master should read first
    if (node->isMaster())
        ProgAngularProjectionMatching::read(argc, (const char **)argv);
    node->barrierWait();
    if (!node->isMaster())
    {

        verbose = 0;//disable verbose for slaves
        ProgAngularProjectionMatching::read(argc, (const char **)argv);
    }
}

/* Define accepted params ------------------------------------------------------------------- */
void MpiProgAngularProjectionMatching::defineParams()
{
    ProgAngularProjectionMatching::defineParams();
    addParamsLine(
        "  [--mpi_job_size <size=10>]   : Number of images sent to a cpu in a single job ");
    addParamsLine("                                : 10 may be a good value");
    addParamsLine(
        "                                : if  -1 the computer will fill the value for you");
    addParamsLine(
        "  [--chunk_angular_distance <dist=-1>]  : sample the projection sphere with this ");
    addParamsLine(
        "                                 :using the voronoi regions");
    addParamsLine(
        "  [--sym <cn=\"c1\">]            : One of the 17 possible symmetries in");
    addParamsLine(
        "                                 :single particle electronmicroscopy");
    addParamsLine(
        "                                 :i.e.  ci, cs, cn, cnv, cnh, sn, dn, dnv,");
    addParamsLine(
        "                                 :dnh, t, td, th, o, oh, i1 (default MDB), i2, i3, i4, ih");
    addParamsLine(
        "                                 :i1h (default MDB), i2h, i3h, i4h");
    addParamsLine(
        "                                : where n may change from 1 to 99");
}

/* Read parameters --------------------------------------------------------- */
void MpiProgAngularProjectionMatching::readParams()
{
    ProgAngularProjectionMatching::readParams();
    mpi_job_size = getIntParam("--mpi_job_size");
    imagesBuffer = new size_t[mpi_job_size + 1];
    chunk_angular_distance = getDoubleParam("--chunk_angular_distance");
    fn_sym = getParam("--sym");
}

void MpiProgAngularProjectionMatching::processAllImages()
{
	int		i=0, nodeID=0;					// Loop counters.
	bool	pending_Images=true;			// More images to process flag.
    std::vector<size_t> imagesToProcess;	// Data to be processed.
    size_t numberOfImages = 0;				// # images to process in TAG_START_JOB message.
    MPI_Status 	status[node->size-1];		// MPI status.

	// Master distribute jobs.
    if (node->isMaster())
    {
    	MPI_Request	requests[node->size-1];							// Non blocking requests.
        size_t 		async_Buffers[node->size-1][mpi_job_size + 1];	// Transmission buffer for slave nodes.
        size_t 		finishedImages = 0, totalImages = DFexp.size();

        if (verbose)
        {
            progress_bar_step = XMIPP_MAX(1, totalImages / 80);
            init_progress_bar(totalImages);
        }

        // Send first asynchronous messages to all nodes.
        for (nodeID=1; nodeID<node->size; nodeID++)
        {
        	// Check if there are more images to process.
        	if (pending_Images)
        	{
				// Check if there are more images to send to nodes.
				if (!distributeJobs(imagesBuffer, nodeID))
				{
					pending_Images = false;
				}
				else
				{
					// Copy data to be transmitted.
					for (i=0; i<(imagesBuffer[0] + 1) ;i++)
					{
						async_Buffers[nodeID-1][i] = imagesBuffer[i];
					}

					// Send asynchronous message to ith-node.
					MPI_Isend( async_Buffers[nodeID-1], imagesBuffer[0] + 1, XMIPP_MPI_SIZE_T, nodeID,
												TAG_START_JOB, MPI_COMM_WORLD, &requests[nodeID-1]);
				}
       		}
       	}

        int	nFinished=0;
        int	finishedNodes[node->size-1];

        // Loop that sends remaining images asynchronously.
        while (pending_Images)
        {
        	// Test how many operations have finished.
        	MPI_Testsome( node->size-1, requests, &nFinished, finishedNodes, status);

        	// Send jobs to idle nodes.
        	i = 0;
        	while ((i<nFinished) && (pending_Images))
        	{
        		// Get operation identifier and add 1 to set node value.
        		nodeID = finishedNodes[i] + 1;

	            if (verbose)
	            {
            		finishedImages += async_Buffers[nodeID-1][0];
	                progress_bar(finishedImages);
	            }

				// Check if there are more images to send to nodes.
				if (!distributeJobs(imagesBuffer, nodeID))
				{
					pending_Images = false;
				}
				else
				{
					// Copy data to be transmitted.
					for (int j=0; j<(imagesBuffer[0] + 1) ;j++)
					{
						async_Buffers[nodeID-1][j] = imagesBuffer[j];
					}

					// Send asynchronous message to ith-node.
					MPI_Isend( async_Buffers[nodeID-1], imagesBuffer[0] + 1, XMIPP_MPI_SIZE_T, nodeID,
												TAG_START_JOB, MPI_COMM_WORLD, &requests[nodeID-1]);
				}

				i++;
        	}

        	// Check if there are pending images for master.
            if (pending_Images)
            {
				// Check if there are more images to send to nodes.
				if (!distributeJobs(imagesBuffer, nodeID))
				{
					pending_Images = false;
				}
				else
				{
					numberOfImages = imagesBuffer[0];

			        // Update output parameter with images to process.
			        imagesToProcess.clear();
			        for (size_t i = 1; i <= numberOfImages; ++i)
			        {
			            imagesToProcess.push_back(imagesBuffer[i]);
			        }

		        	// Process new images.
		            processSomeImages(imagesToProcess);

		            if (verbose)
		            {
			            finishedImages += numberOfImages;
		                progress_bar(finishedImages);
		            }
				}
            }
        }

        // Wait for all nodes until they finish their last job.
        MPI_Waitall ( node->size-1, requests, status);

        // Send asynchronous messages to all nodes.
        for (nodeID=1; nodeID<node->size; nodeID++)
        {
			// Send asynchronous STOP message to ith-node (does not matter buffer content).
			MPI_Isend( async_Buffers[nodeID-1], imagesBuffer[0] + 1, XMIPP_MPI_SIZE_T, nodeID,
											TAG_STOP_JOB, MPI_COMM_WORLD, &requests[nodeID-1]);
        }

        // Wait for all nodes.
        MPI_Waitall ( node->size-1, requests, status);
    }
    // Slaves really work.
    else
    {
    	bool 	finished = false;				// Loop control flag.
        MPI_Request	requests;					// Asynchronous request.

        // Get first message from master node.
        MPI_Recv(imagesBuffer, mpi_job_size + 1, XMIPP_MPI_SIZE_T, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status[0]);

        finished = false;
        while(!finished)
        {
			// TAG_JOB_REPLY -> more images to process.
			if (status[0].MPI_TAG == TAG_START_JOB)
			{
				// Copy data from buffer.
				numberOfImages = imagesBuffer[0];
				imagesToProcess.clear();
				for (size_t i = 1; i <= numberOfImages; ++i)
				{
					imagesToProcess.push_back(imagesBuffer[i]);
				}

				// Start new reception.
				MPI_Irecv( imagesBuffer, mpi_job_size + 1, XMIPP_MPI_SIZE_T, 0, MPI_ANY_TAG,
																			MPI_COMM_WORLD, &requests);

				// Process new images.
				processSomeImages(imagesToProcess);

	    		// Wait for reception.
	    		MPI_Wait( &requests, &status[0]);
			}
			// Slave finishes.
			else
			{
				finished = true;
			}
        }
        // At this point slave node finishes its execution.
    }
}

bool MpiProgAngularProjectionMatching::distributeJobs(size_t * imagesToSent,
        int node)
{
    int node_index = last_chunk[node];
    static bool reached_last_chunk = false;

    if (node_index == -1
        || chunk_mysampling.my_exp_img_per_sampling_point[node_index].empty())
    {

        if (!reached_last_chunk)
        {
            node_index = chunk_index;
            chunk_index = (chunk_index + 1) % chunk_number;
            reached_last_chunk = (chunk_index == 0);
        }
        else
        {
            int i = 0;
            for (;
                 i < chunk_number
                 && chunk_mysampling.my_exp_img_per_sampling_point[chunk_index].empty();
                 ++i, chunk_index = (chunk_index + 1) % chunk_number)
                ;

            //            std::cerr << formatString("DEBUG: master: assigned %d chuck to node %d, i = %d", chunk_index, node, i);
            if (i == chunk_number)
            {
                imagesToSent[0] = 0;
                return false;
            }
            last_chunk[node] = node_index = chunk_index;
            chunk_index = (chunk_index + 1) % chunk_number;
        }
    }

    size_t assigned_images = std::min((size_t)mpi_job_size, chunk_mysampling.my_exp_img_per_sampling_point[node_index].size());
    imagesToSent[0] = assigned_images;

    //std::cerr << formatString("DEBUG: master: assigned %lu images to node %d\n", assigned_images, node);

    for (size_t i = 1; i <= assigned_images; ++i)
    {
        imagesToSent[i] =
            chunk_mysampling.my_exp_img_per_sampling_point[node_index].back() + FIRST_IMAGE;
        //        std::cerr << " " << imagesToSent[i];
        chunk_mysampling.my_exp_img_per_sampling_point[node_index].pop_back();
    }

    //    std::cerr << std::endl;

    return true;
}


void MpiProgAngularProjectionMatching::writeOutputFiles()
{
    node->gatherMetadatas(DFo, fn_out);
    if (node->isMaster())
    {
    	MetaData mdAux;
    	mdAux.sort(DFo,MDL_IMAGE);
    	mdAux.write(fn_out, do_overwrite);
    }
}

void MpiProgAngularProjectionMatching::produceSideInfo()
{
    ProgAngularProjectionMatching::produceSideInfo();

    if (node->isMaster())
        computeChunks();
    node->barrierWait();
}

void MpiProgAngularProjectionMatching::computeChunks()
{
	size_t max_number_of_images_in_around_a_sampling_point = 0;
    //process the symmetry file
    if (!chunk_mysampling.SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(ERR_NUMERICAL,
                     (String)"mpi_angular_proj_match::prerun Invalid symmetry: " + fn_sym);
    chunk_mysampling.SL.readSymmetryFile(fn_sym);
    // find a value for chunk_angular_distance if != -1
    if (chunk_angular_distance == -1)
        computeChunkAngularDistance(symmetry, sym_order);
    //store symmetry matrices, this is faster than computing them
    //each time we need them
    chunk_mysampling.fillLRRepository();
    int remaining_points = 0;

    while (remaining_points == 0)
    {
        //first set sampling rate
        chunk_mysampling.setSampling(chunk_angular_distance);
        //create sampling points in the whole sphere
        chunk_mysampling.computeSamplingPoints(false, 0, 180);
        //precompute product between symmetry matrices
        //and experimental data
        chunk_mysampling.fillExpDataProjectionDirectionByLR(fn_exp);
        //remove redundant sampling points: symmetry
        chunk_mysampling.removeRedundantPoints(symmetry, sym_order);
        remaining_points =
            chunk_mysampling.no_redundant_sampling_points_angles.size();
        if (chunk_angular_distance > 2)
            chunk_angular_distance -= 1;
        else
            chunk_angular_distance /= 2;
        //if(remaining_points ==0)
        if (verbose)
            std::cout << "New chunk_angular_distance " << chunk_angular_distance
            << std::endl;
        if (chunk_angular_distance < 0)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "Can't compute chunk_angular_distance");
    }
    //remove sampling points too far away from experimental data
    chunk_mysampling.removePointsFarAwayFromExperimentalData();
    //for each sampling point find the experimental images
    //closer to that point than to any other
    chunk_mysampling.findClosestExperimentalPoint();
    //print number of points per node
    chunk_number = chunk_mysampling.my_exp_img_per_sampling_point.size();
    for (int j = 0; j < chunk_number; j++)
    {
        if (max_number_of_images_in_around_a_sampling_point
            < chunk_mysampling.my_exp_img_per_sampling_point[j].size())
            max_number_of_images_in_around_a_sampling_point =
                chunk_mysampling.my_exp_img_per_sampling_point[j].size();
    }
    if (verbose)
    {
        std::cout << "number of subsets: " << chunk_number << std::endl
        << "biggest subset (EXPERIMENTAL images per chunk): "
        << max_number_of_images_in_around_a_sampling_point << std::endl
        << "maximun number of references in memory: "
        << max_nr_refs_in_memory << std::endl;
    }
    //alloc memory for buffer
    if (mpi_job_size == -1)
        mpi_job_size = (int)ceil((double)DFexp.size()/(node->size - 1));

    //Distribution related variables
    chunk_index = 0;
    last_chunk = new int[node->size];
    for (size_t i = 1; i < node->size; ++i)
        last_chunk[i] = -1;//by default no chunk assigned yet
}

void MpiProgAngularProjectionMatching::computeChunkAngularDistance(int symmetry,
        int sym_order)
{
    double non_reduntant_area_of_sphere =
        chunk_mysampling.SL.nonRedundantProjectionSphere(symmetry,
                sym_order);
    double number_cpus = (double) node->size - 1;
    //NEXT ONE IS SAMPLING NOT ANOTHERSAMPLING
    double neighborhood_radius = fabs(acos(mysampling.cos_neighborhood_radius));
    //NEXT ONE IS SAMPLING NOT ANOTHERSAMPLING
    if (mysampling.cos_neighborhood_radius < -1.001)
        neighborhood_radius = 0;
    int counter = 0;
    while (1)
    {
        if (counter++ > 1000)
        {
            chunk_angular_distance = 0.001;
            std::cerr << "****************************************************"
            << std::endl;
            std::cerr << "* WARNING: The neighbourhood does not fit in memory "
            << std::endl;
            std::cerr << "****************************************************"
            << std::endl;
            break;
        }
        double area_chunk = non_reduntant_area_of_sphere / number_cpus;
        //area chunk is area of spheric casket=2 PI h
        chunk_angular_distance = acos(1 - area_chunk / (2 * PI));
        double area_chunck_neigh = 2 * PI
                                   * (1 - cos(chunk_angular_distance + neighborhood_radius));
        //double area_chunck= 2 * PI *( 1 - cos(chunk_angular_distance));
        //let us see how many references from the reference library fit
        //in area_chunk, that is divide area_chunk between the voronoi
        //region of the sampling points of the reference library
        double areaVoronoiRegionReferenceLibrary = 2 *( 3 *(  acos(
                    //NEXT ONE IS SAMPLING NOT ANOTHERSAMPLING
                    cos(mysampling.sampling_rate_rad)/(1+cos(mysampling.sampling_rate_rad)) )  ) - PI);
        int number_of_images_that_fit_in_a_chunck_neigh =(int)
            ceil(area_chunck_neigh / areaVoronoiRegionReferenceLibrary);
        //#define DEBUG
#ifdef DEBUG

        std::cerr << "\n\ncounter " << counter << std::endl;
        std::cerr << "area_chunk " << area_chunk << std::endl;
        std::cerr << "2*chunk_angular_distance " << 2*chunk_angular_distance << std::endl;
        //NEXT ONE IS SAMPLING NOT ANOTHERSAMPLING
        std::cerr << "sampling_rate_rad " << mysampling.sampling_rate_rad
        << " " << mysampling.sampling_rate_rad*180/PI
        << std::endl;
        std::cerr << "neighborhood_radius " << neighborhood_radius
        << std::endl;
        std::cerr << "areaVoronoiRegionReferenceLibrary " << areaVoronoiRegionReferenceLibrary << std::endl;
        std::cerr << "number_of_images_that_fit_in_a_chunck_neigh " << number_of_images_that_fit_in_a_chunck_neigh << std::endl;
        std::cerr << "number_cpus " << number_cpus << std::endl;
        std::cerr << "max_nr_imgs_in_memory " << max_nr_imgs_in_memory << std::endl;
#endif
#undef DEBUG

        if (number_of_images_that_fit_in_a_chunck_neigh > max_nr_imgs_in_memory)
            number_cpus = 1.2 * number_cpus;
        else
            break;
    }
    //chunk_angular_distance -= neighborhood_radius;
    chunk_angular_distance *= 2.0;

    //#define DEBUG
#ifdef DEBUG

    std::cerr << "chunk_angular_distance " << chunk_angular_distance
    << std::endl
    << "neighborhood_radius " << neighborhood_radius
    << std::endl;
#endif
#undef DEBUG
    //chuck should not be bigger than a triangle in the icosahedra
    if (chunk_angular_distance >= 0.5 * cte_w)
        chunk_angular_distance = 0.5 * cte_w;
    chunk_angular_distance *= (180. / PI);
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "chunk_angular_distance_degrees " << chunk_angular_distance
    << std::endl;
#endif
#undef DEBUG

}
