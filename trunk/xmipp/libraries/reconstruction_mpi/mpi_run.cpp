/***************************************************************************
 *
 * Authors:
 *
 * Roberto Marabini
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

#include "mpi_run.h"

#include <data/args.h>

#include <cstring>
#include <cstdlib>

#define TAG_WORK   0
#define TAG_STOP   1
#define TAG_WAIT   2

/* Empty constructor ------------------------------------------------------- */
Prog_MPI_Run_Parameters::Prog_MPI_Run_Parameters(int argc, char **argv)
{
    MPI_Comm_size(MPI_COMM_WORLD, &(nprocs));
    MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
    if (nprocs < 2)
        error_exit("This program cannot be executed in a single working node");
    //Blocks until all process have reached this routine.
    MPI_Barrier(MPI_COMM_WORLD);
}


/* Read parameters --------------------------------------------------------- */
void Prog_MPI_Run_Parameters::read(int argc, char **argv)
{
    fn_commands = get_param(argc, argv, "-i");
}

/* Usage ------------------------------------------------------------------- */
void Prog_MPI_Run_Parameters::usage()
{
    cerr << "MPI_Run\n"
    << "   -i <command file>    : File with commands to send to mpirun\n"
    << "\n"
    << "Example of use:\n"
    << "   xmipp_mpi_run -i commandd_file\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_MPI_Run_Parameters::show()
{
    cout << "Commands  file:           " << fn_commands << endl
    ;
}


/* Run --------------------------------------------------------------------- */
void Prog_MPI_Run_Parameters::run()
{


    if (rank == 0)
    {
        ifstream fh_in;
        fh_in.open(fn_commands.c_str());
        if (!fh_in)
            REPORT_ERROR(1, (string)"Cannot open " + fn_commands);

#define MAX_LINE 1024
        string line;
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
            string::size_type loc = line.find("MPI_Barrier", 0);
            if (loc != string::npos)
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
            //cout << line << endl;
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
            cout << szline << endl;
            //do the job
            system(szline);


        }
    }

    MPI_Finalize();
}

/* a short function to print a message and exit */
void Prog_MPI_Run_Parameters::error_exit(char * msg)
{
    fprintf(stderr, "%s", msg);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}
