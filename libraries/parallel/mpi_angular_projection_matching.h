/***************************************************************************
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              J.M de la Rosa   (jmdelarosa@cnb.csic.es)
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

#ifndef MPI_ANGULAR_PROJECTION_MATCHING_H_
#define MPI_ANGULAR_PROJECTION_MATCHING_H_

#include "parallel/xmipp_mpi.h"
#include "reconstruction/angular_projection_matching.h"

/**@defgroup MpiProgAngularProjectionMatching  MpiProgAngularProjectionMatching
   @ingroup Programs */
//@{
/** Program to parallelize the ML 2D alignment program */
class MpiProgAngularProjectionMatching: public ProgAngularProjectionMatching
{
private:
    MpiNode *node;

    /** Dvide the job in this number block with this number of images */
    int mpi_job_size;
    /** Buffer to send and receive images ids */
    size_t * imagesBuffer;
    /** Indexes to distribute jobs to slaves */
    int chunk_index, chunk_number; //will keep the next chunk to assign and the total number of chunks
    int *last_chunk; //store the last chunk assigned to each node


    /** classify the experimental data making voronoi regions
        with an hexagonal grid mapped onto a sphere surface */
    double chunk_angular_distance;

    /** symmetry file */
    FileName        fn_sym;

    /** sampling object */
    Sampling chunk_mysampling;

    /** Symmetry. One of the 17 possible symmetries in
        single particle electron microscopy.
         */
    int symmetry;

    /** For infinite groups symmetry order*/
    int sym_order;

public:
    /** Redefine read */
    void read(int argc, char** argv);
    /** Constructor */
    MpiProgAngularProjectionMatching();
    /** Destructor */
    ~MpiProgAngularProjectionMatching();

    /** Override virtual function implementations */
    void processAllImages();
    void writeOutputFiles();
    /** Function to distribute jobs to slaves.
     * will return false if no more images to process.
     * Will try to sent neighbours images to same node.
     */
    bool distributeJobs(size_t * imagesToSent, int node);
    /** Function of slaves nodes to ask jobs to master.
     * will return false if no more images to process.
     */
    bool requestJobs(std::vector<size_t> &imagesToProcess);

    /* Define accepted params ------------------------------------------------------------------- */
    void defineParams();
        /* Read parameters --------------------------------------------------------- */
    void readParams();

    /** Redefine produceSideInfo */
    void produceSideInfo();
    /** These two function will be executed only by master */
    void computeChunks();
    void computeChunkAngularDistance(int symmetry, int sym_order);
}
;//end of class MpiProgML2D
/** @} */
#endif /* MPI_ANGULAR_PROJECTION_MATCHING_H_ */
