/***************************************************************************
 *
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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


#include "parallel/mpi_project_XR.h"


int main(int argc, char *argv[])
{

    MpiNode node(argc, argv);

//    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
//    {
//        fprintf(stderr, "MPI initialization error\n");
//        exit(EXIT_FAILURE);
//    }
    ProgProjectXR      prog_prm;
    Projection                           proj;
    MetaData                              SF;

    // Check the command line
    try
    {
        prog_prm.read(argc, argv);

    }
    catch (XmippError &XE)
    {
        std::cout << XE;
        prog_prm.usage();
        exit(1);
    }

    try
    {
        // Really project
        ROUT_mpi_XR_project(prog_prm, proj, SF, node);
    }
    catch (XmippError XE)
    {
        std::cout << XE;
    }
    return 0;
}

