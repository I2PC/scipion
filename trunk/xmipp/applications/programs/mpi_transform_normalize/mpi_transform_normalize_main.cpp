/***************************************************************************
 *
 * Authors:  Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 * Lab. de Bioingenieria, Univ. San Pablo CEU
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

#include <data/normalize.h>
#include <parallel/mpi.h>

/** Class to normalize with  MPI parallelization */
class MpiProgTransformNormalize: public ProgNormalize,
            public MpiMetadataProgram
{
public:
    /** Redefine read to initialize MPI environment */
    void read(int argc, char **argv)
    {
    	MpiMetadataProgram::read(argc,argv);
    	if (!node->isMaster())
    		verbose=0;
        ProgNormalize::read(argc, argv);
    }

    /** Preprocess */
    void preProcess()
    {
    	ProgNormalize::preProcess();
        createTaskDistributor(mdIn);
    }

    //Only master do starting progress bar stuff
    void startProcessing()
    {
        if (node->isMaster())
        	ProgNormalize::startProcessing();
    }
    //Only master show progress
    void showProgress()
    {
        if (node->isMaster())
        {
        	time_bar_done=first+1;
        	ProgNormalize::showProgress();
        }
    }

    //Now use the distributor to grasp images
    bool getImageToProcess(size_t &objId, size_t &objIndex)
    {
        return getTaskToProcess(objId, objIndex);
    }

    void finishProcessing()
    {
        node->gatherMetadatas(mdOut, fn_out);
        if (node->isMaster())
        	ProgNormalize::finishProcessing();
    }
}
;

int main(int argc, char **argv)
{
	MpiProgTransformNormalize program;
    program.read(argc, argv);
    return program.tryRun();
}
