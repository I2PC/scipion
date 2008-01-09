/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)
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

#include <reconstruction/angular_projection_matching.h>

#include <mpi.h>

int main(int argc, char **argv)
{

    // For parallelization
    int rank, size, num_img_tot;

    double                        aux, sumCC;
    Matrix2D<double>              Maux;
    FileName                      fn_img, fn_tmp, fn_base;
    DocFile                       DFo;
    SelFile                       SFo;
    Prog_projection_matching_prm  prm;

    // Init Parallel interface
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get input parameters
    try
    {
        // Read command line & produce side info
        prm.read(argc, argv);
        if (rank == 0) prm.show();
        else
        {
            prm.verb = 0;
            prm.output_refs = false;
        }

        prm.produce_Side_info();
        // Select only relevant part of selfile for this rank
        prm.SF.mpi_select_part(rank, size, num_img_tot);


    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    try
    {

        DFo.clear();
        if (rank == 0)
            DFo.append_comment("Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Refno (6), maxCC (7), Z-score (8)");

        // Process all images
        prm.PM_loop_over_all_images(prm.SF, DFo, sumCC);

        // Here MPI_allreduce
        MPI_Allreduce(&sumCC, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumCC = aux;
        if (prm.verb > 0) cerr << " Average maxCC = " << sumCC / num_img_tot << endl;

        if (prm.output_classes)
        {
            // Do MPI_Allreduce of the averages (and their weight)
            // And write out temporary selfiles
            Maux.resize(prm.dim, prm.dim);
            Maux.setXmippOrigin();
            fn_base.compose(prm.fn_root, rank, "");
            for (int dirno = 0; dirno < prm.nr_dir; dirno++)
            {
                MPI_Allreduce(MULTIDIM_ARRAY(prm.class_avgs[dirno]()),
                              MULTIDIM_ARRAY(Maux),
                              MULTIDIM_SIZE(prm.class_avgs[dirno]()),
                              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                prm.class_avgs[dirno]() = Maux;
                double weight = prm.class_avgs[dirno].weight();
                MPI_Allreduce(&weight, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                prm.class_avgs[dirno].weight() = aux;
                fn_img.compose(fn_base, dirno, "tmpsel");
                prm.class_selfiles[dirno].write(fn_img);
            }
        }

        // All nodes write out temporary DFo
        fn_img.compose(prm.fn_root, rank, "tmpdoc");
        DFo.write(fn_img);
        MPI_Barrier(MPI_COMM_WORLD);

        // Master writes out final docfile
        if (rank == 0)
        {
            DFo.clear();
            for (int rank2 = 0; rank2 < size; rank2++)
            {
                fn_img.compose(prm.fn_root, rank2, "tmpdoc");
                int ln = DFo.LineNo();
                DFo.append(fn_img);
                DFo.locate(DFo.get_last_key());
                DFo.next();
                DFo.remove_current();
                system(((string)"rm -f " + fn_img).c_str());
            }
            fn_tmp = prm.fn_root + ".doc";
            DFo.write(fn_tmp);
        }

        // Master writes out all class averages and combines all selfiles
        if (prm.output_classes && rank == 0)
        {
            for (int dirno = 0; dirno < prm.nr_dir; dirno++)
            {
                SFo.clear();
                for (int rank2 = 0; rank2 < size; rank2++)
                {
                    fn_base.compose(prm.fn_root, rank2, "");
                    fn_img.compose(fn_base, dirno, "tmpsel");
                    SFo.append(fn_img);
                    system(((string)"rm -f " + fn_img).c_str());
                }
                prm.class_selfiles[dirno] = SFo;
            }
            prm.write_classes();
        }

    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    MPI_Finalize();
    return 0;

}




