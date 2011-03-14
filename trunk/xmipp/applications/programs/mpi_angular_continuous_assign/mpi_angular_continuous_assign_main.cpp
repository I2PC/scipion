/***************************************************************************
 *
 * Authors:  Slavica Jonic slavica.jonic@impmc.jussieu.fr
 *           Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
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
//#include <mpi.h>
#include <reconstruction/angular_continuous_assign.h>
#include <parallel/mpi.h>

/** Class to perfom the NMA Alignment with  MPI parallelization */
class MpiProgAngularContinuousAssign: public ProgAngularContinuousAssign,
            public MpiMetadataProgram
{
public:
    /** Redefine read to initialize MPI environment */
    void read(int argc, char **argv)
    {
        ProgAngularContinuousAssign::read(argc, argv);
        MpiMetadataProgram::read(argc,argv);
    }

    /** Preprocess */
    void preProcess()
    {
        ProgAngularContinuousAssign::preProcess();
        createTaskDistributor(mdIn);
    }

    /** Postprocess */
    void postProcess()
    {
        node->gatherMetadatas(mdIn, fn_out);
        mdIn.write(fn_out);
    }

    //Only master do starting progress bar stuff
    void startProcessing()
    {
        if (node->isMaster())
            ProgAngularContinuousAssign::startProcessing();
    }
    //Only master show progress
    void showProgress()
    {
        if (node->isMaster())
        {
        	time_bar_done=first+1;
            ProgAngularContinuousAssign::showProgress();
        }
    }

    //Now use the distributor to grasp images
    size_t getImageToProcess()
    {
        return getTaskToProcess();
    }

    void finishProcessing()
    {
        node->barrierWait();
        if (node->isMaster())
            ProgAngularContinuousAssign::finishProcessing();
    }
}
;

int main(int argc, char **argv)
{
    MpiProgAngularContinuousAssign program;
    program.read(argc, argv);
    return program.tryRun();
}
