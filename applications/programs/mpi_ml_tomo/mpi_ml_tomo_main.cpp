/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.csic.es)
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

//#include <mpi.h>
#include <parallel/mpi_ml_align2d.h>

RUN_XMIPP_PROGRAM(MpiProgMLTomo)

/*
#define TAG_DOCFILE 12
#define TAG_DOCFILESIZE 13

int main(int argc, char **argv)
{
    Prog_ml_tomo_prm prm;

    int c, nn, imgno, opt_refno;
    double LL, sumw_allrefs, convv, sumcorr;
    bool converged;
    std::vector<double> conv;
    double aux, wsum_sigma_noise, wsum_sigma_offset;
    std::vector<Matrix3D<double > > wsumimgs;
    std::vector<Matrix3D<double > > wsumweds;
    std::vector<Matrix1D<double > > fsc;
    Matrix1D<double> sumw,  Vaux;
    Matrix3D<double> P_phi, Mr2, Maux, Maux2;
    FileName fn_img, fn_tmp;
    Matrix1D<double> oneline(0);
    DocFile DFo;
    // For parallelization
    int rank, size, num_img_tot;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    // Get input parameters
    try
    {

        // Read subsequently to avoid problems in restart procedure
        for (int proc = 0; proc < size; proc++)
        {
            if (proc == rank) prm.read(argc, argv);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if (rank != 0) prm.verb = 0;

        // All nodes produce general side-info
        prm.produceSideInfo();

        // Some output to screen
        if (rank == 0) prm.show();

        // Create references from random subset averages, or read them from selfile
        if (prm.fn_ref == "")
        {
            if (prm.nr_ref != 0)
            {
                if (rank == 0)
                {
                    prm.generateInitialReferences();
                }
                else
                {
                    prm.fn_ref = prm.fn_root + "_it";
                    prm.fn_ref.compose(prm.fn_ref, 0, "sel");
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            else
            {
                REPORT_ERROR(ERR_ARG_MISSING, "Please provide -ref or -nref");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Select only relevant part of selfile for this rank
        mpiSelectPart(prm.SF, rank, size, num_img_tot);

        // And produce selfile-specific side-info
        prm.produceSideInfo2();
        MPI_Barrier(MPI_COMM_WORLD);

    }
    catch (XmippError XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    try
    {

        // For MPI_Allreduce statements
        Maux.resize(prm.dim, prm.dim, prm.dim);
        Maux.setXmippOrigin();
        Maux2.resize(prm.dim, prm.dim, prm.hdim+1);

        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {
            if (prm.verb > 0) std::cerr << "  Multi-reference refinement:  iteration " << iter << " of " << prm.Niter << std::endl;

            // Save old reference images
            for (int refno = 0;refno < prm.nr_ref; refno++)
                prm.Iold[refno]() = prm.Iref[refno]();

            // Integrate over all images
            prm.expectation(prm.SF, prm.Iref, iter,
                            LL, sumcorr, DFo, wsumimgs, wsumweds,
                            wsum_sigma_noise, wsum_sigma_offset,
                            sumw);

            // Here MPI_allreduce of all wsums,LL and sumcorr !!!
            MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            LL = aux;
            MPI_Allreduce(&sumcorr, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumcorr = aux;
            MPI_Allreduce(&wsum_sigma_noise, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_noise = aux;
            MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_sigma_offset = aux;
            Vaux.resize(prm.nr_ref);
            MPI_Allreduce(MULTIDIM_ARRAY(sumw), MULTIDIM_ARRAY(Vaux),
                          MULTIDIM_SIZE(sumw), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumw = Vaux;
            for (int refno = 0; refno < 2*prm.nr_ref; refno++)
            {
                MPI_Allreduce(MULTIDIM_ARRAY(wsumimgs[refno]), MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(wsumimgs[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                wsumimgs[refno] = Maux;
                if (prm.do_missing)
                {
                    MPI_Allreduce(MULTIDIM_ARRAY(wsumweds[refno]), MULTIDIM_ARRAY(Maux2),
                                  MULTIDIM_SIZE(wsumweds[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    wsumweds[refno] = Maux2;
                }
            }

            // Update model parameters
            prm.maximization(wsumimgs, wsumweds,
                             wsum_sigma_noise, wsum_sigma_offset,
                             sumw, sumcorr, sumw_allrefs, fsc, iter);

            // Check convergence
            converged = prm.checkConvergence(conv);

            // Write intermediate files
            if (rank != 0)
            {
                // All slaves send docfile to the master
                std::ostringstream doc;
                doc << DFo;
                int s_size=  doc.str().size();
                char results[s_size];
                strncpy(results,doc.str().c_str(),s_size);
                results[s_size]='\0';
                MPI_Send(&s_size, 1, MPI_INT, 0, TAG_DOCFILESIZE, MPI_COMM_WORLD);
                MPI_Send(results, s_size, MPI_CHAR, 0, TAG_DOCFILE, MPI_COMM_WORLD);
            }
            else
            {
                // Master fills docfile
                std::ofstream myDocFile;
                FileName fn_tmp;
                fn_tmp.compose(prm.fn_root + "_it",iter,"doc");
                myDocFile.open (fn_tmp.c_str());
                myDocFile << " ; Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Zoff (6), Ref (7), Wedge (8), Pmax/sumP (9), LL (10)\n";

                // Master's own contribution
                myDocFile << DFo;
                int docCounter=1;
                while (docCounter < size)
                {
                    // receive in order
                    int iNumber, s_size;
                    MPI_Recv(&s_size, 1, MPI_INT, docCounter, TAG_DOCFILESIZE, MPI_COMM_WORLD, &status);
                    char results[s_size];
                    MPI_Recv(results, s_size, MPI_CHAR, docCounter, TAG_DOCFILE, MPI_COMM_WORLD, &status);
                    results[s_size]='\0';
                    myDocFile<<results ;
                    docCounter++;
                }

                //save doc_file and renumber it
                myDocFile.close();
                DFo.clear();
                DFo.read(fn_tmp);
                DFo.renum();

                // Output all intermediate files
                prm.writeOutputFiles(iter, DFo, wsumweds, sumw_allrefs, LL, sumcorr, conv, fsc);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            if (converged)
            {
                if (prm.verb > 0) std::cerr << " Optimization converged!" << std::endl;
                break;
            }
            else
            {
                // reset DFo
                DFo.clear();
            }
            MPI_Barrier(MPI_COMM_WORLD);

        } // end loop iterations
	if (rank == 0)
	    prm.writeOutputFiles(-1, DFo, wsumweds, sumw_allrefs, LL, sumcorr, conv, fsc);

    }
    catch (XmippError XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    MPI_Finalize();
    return 0;
}
*/
