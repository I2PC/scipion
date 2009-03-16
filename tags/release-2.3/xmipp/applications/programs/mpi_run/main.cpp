/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.uam.es)
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

//#include "mpi_run.h"

#include <mpi.h>
#include <data/args.h>

#include <cstring>
#include <cstdlib>
#include <data/funcs.h>

#define TAG_WORK   0
#define TAG_STOP   1
#define TAG_WAIT   2
class Prog_MPI_Run_Parameters
{
    public:
        /** PDB file */
        FileName fn_commands;

        /** Number of Procesors **/
        int nprocs;

        /** computing node number. Master=0 */
        int rank;

        /** status after am MPI call */
        MPI_Status status;
    /** Empty constructor */
    //Prog_MPI_Run_Parameters(int argc, char **argv);

    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    //void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introducing this parameters. */
    //void usage();

    /** Show parameters. */
    //void show();

    /** Run. */
    //void run();

    /* Empty constructor ------------------------------------------------------- */
    Prog_MPI_Run_Parameters(int argc, char **argv)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &(nprocs));
        MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
        if (nprocs < 2)
            error_exit("This program cannot be executed in a single working node");
        //Blocks until all process have reached this routine.
        MPI_Barrier(MPI_COMM_WORLD);
    }


    /* Read parameters --------------------------------------------------------- */
    void read(int argc, char **argv)
    {
        fn_commands = getParameter(argc, argv, "-i");
    }

    /* Usage ------------------------------------------------------------------- */
    void usage()
    {
        std::cerr << "MPI_Run\n"
        << "   -i <command file>    : File with commands to send to mpirun\n"
        << "\n"
        << "Example of use:\n"
        << "   xmipp_mpi_run -i commandd_file\n"
        ;
    }

    /* Show -------------------------------------------------------------------- */
    void show()
    {
        std::cout << "Commands  file:           " << fn_commands << std::endl
        ;
    }


    /* Run --------------------------------------------------------------------- */
    void run()
    {


        if (rank == 0)
        {
            std::ifstream fh_in;
            fh_in.open(fn_commands.c_str());
            if (!fh_in)
                REPORT_ERROR(1, (std::string)"Cannot open " + fn_commands);

    #define MAX_LINE 1024
            std::string line;
            char szline[MAX_LINE];
            int number_of_node_waiting = 0; // max is nprocs -1
            while (!fh_in.eof())
            {
                //wait untill a server is free
                MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0,
                         MPI_COMM_WORLD, &status);
                number_of_node_waiting++;
                getline(fh_in, line);
                strcpy(szline, line.c_str());
                std::string::size_type loc = line.find("MPI_Barrier", 0);
                if (loc != std::string::npos)
                {
                    while (number_of_node_waiting < (nprocs - 1))
                    {
                        MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0,
                                 MPI_COMM_WORLD, &status);
                        number_of_node_waiting++;
                    }
                    while (number_of_node_waiting > 0)
                    {
                        MPI_Send(&szline,
                                 1,
                                 MPI_CHAR,
                                 number_of_node_waiting,
                                 TAG_WAIT,
                                 MPI_COMM_WORLD);
                        number_of_node_waiting--;
                    }
                    continue;
                }
                //send work
                MPI_Send(&szline,
                         MAX_LINE,
                         MPI_CHAR,
                         status.MPI_SOURCE,
                         TAG_WORK,
                         MPI_COMM_WORLD);
                number_of_node_waiting--;
                //std::cout << line << std::endl;
            }

            fh_in.close();
            for (int i = 1; i < nprocs; i++)
            {
                MPI_Send(0, 0, MPI_INT, i, TAG_STOP, MPI_COMM_WORLD);
            }

        }
        else
        {

            while (1)
            {
                char szline[1024];
                //any message from de master, is tag is TAG_STOP then stop
                MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
                //get yor next task
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == TAG_STOP)//I am free
                {
                    break;
                }
                if (status.MPI_TAG == TAG_WAIT)//I am free
                {
                    MPI_Recv(&szline, 1, MPI_CHAR, 0, TAG_WAIT, MPI_COMM_WORLD, &status);
                    continue;
                }
                MPI_Recv(&szline, MAX_LINE, MPI_CHAR, 0, TAG_WORK, MPI_COMM_WORLD, &status);
                std::cout << szline << std::endl;
                //do the job
                system(szline);


            }
        }

        MPI_Finalize();
    }

    /* a short function to print a message and exit */
    void error_exit(char * msg)
    {
        fprintf(stderr, "%s", msg);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

};

int main(int argc, char *argv[])
{
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        fprintf(stderr, "MPI initialization error\n");
        exit(EXIT_FAILURE);
    }

    Prog_MPI_Run_Parameters prm(argc, argv);
    try
    {
        prm.read(argc, argv);
    }

    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(1);
    }

    try
    {
        if (prm.rank == 0) prm.show();
        prm.run();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    exit(0);
}


