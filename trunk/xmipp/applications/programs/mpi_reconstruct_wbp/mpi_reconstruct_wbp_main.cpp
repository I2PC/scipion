/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao Castro (jrbcast@cnb.csic.es)
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

#include <data/mpi.h>

#include <reconstruction/reconstruct_wbp.h>

class ProgMPIRecWbp: public ProgRecWbp, public MpiMetadataProgram
{
public:
    void read(int argc, char **argv)
    {
        MpiMetadataProgram::read(argc,argv);
        if (!node->isMaster())
            verbose=0;
        ProgRecWbp::read(argc, argv);
    }
    void preProcess()
    {
        ProgRecWbp::preProcess();
        createTaskDistributor(mdIn);
    }
    void startProcessing()
    {
        if (node->isMaster())
            ProgRecWbp::startProcessing();
    }
    void showProgress()
    {
        if (node->isMaster())
        {
            time_bar_done=first+1;
            ProgRecWbp::showProgress();
        }
    }
    bool getImageToProcess(size_t &objId, size_t &objIndex)
    {
        return getTaskToProcess(objId, objIndex);
    }
    void finishProcessing()
    {
        node->gatherMetadatas(mdOut, fn_out);
        if (node->isMaster())
            ProgRecWbp::finishProcessing();
    }
}
#ifdef NEVER
{
public:
    Image<double> vol, aux;
    int iaux;

public:
    /*  constructor ------------------------------------------------------- */
    ProgMPIRecWbp(int argc, char **argv)
    {
        //parent class constructor will be called by deault without parameters
        MPI_Comm_size(MPI_COMM_WORLD, &(size));
        MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
        if (size < 2)
            error_exit("This program cannot be executed in a single working node");

        MPI_Barrier(MPI_COMM_WORLD);
    }

    /* a short function to print a message and exit */
    void error_exit(char * msg)
    {
        fprintf(stderr, "%s", msg);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    void run()
    {
        // First produce the filter and then cut selfile in smaller parts
        produceSideInfo();

        // Select only relevant part of selfile for this rank
        mpiSelectPart(SF, rank, size, num_img_tot);

        // Actual backprojection
        apply_2Dfilter_arbitrary_geometry(SF, vol());

        aux().resize(vol());
        MPI_Allreduce(MULTIDIM_ARRAY(vol()), MULTIDIM_ARRAY(aux()),
                      MULTIDIM_SIZE(vol()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        vol = aux;
        MPI_Allreduce(&count_thr, &iaux, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0)
        {
            std::cerr << "Fourier pixels for which the threshold was not reached: "
            << (float)(iaux*100.) / (num_img_tot*dim*dim) << " %" << std::endl;
            vol.write(fn_out);
        }
    }
};


int main(int argc, char **argv)
{

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        fprintf(stderr, "MPI initialization errorn");
        exit(EXIT_FAILURE);
    }

    ProgMPIRecWbp program;

    if (program.rank == 0)
    {
        try
        {
            program.read(argc, argv);
        }
        catch (XmippError XE)
        {
            std::cerr << XE;
            MPI_Finalize();
            exit(1);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (program.rank != 0)
    {
        try
        {
            program.read(argc, argv);
        }
        catch (XmippError XE)
        {
            std::cerr << XE;
            MPI_Finalize();
            exit(1);
        }
    }

    try
    {
        program.run();
        MPI_Finalize();
    }
    catch (XmippError XE)
    {
        std::cerr << XE;
        exit(1);
    }

    exit(0);

}
#endif
