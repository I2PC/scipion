/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

#include <parallel/xmipp_mpi.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_program.h>

#define TAG_WORK   0
#define TAG_STOP   1
#define TAG_WAIT   2

class ProgMPIRun: public XmippProgram
{
public:
    /** Command file */
    FileName fn_commands;

    /** status after am MPI call */
    MPI_Status status;

    // Mpi node
    MpiNode *node;

    ProgMPIRun(int argc, char **argv)
    {
        node=new MpiNode(argc,argv);
        if (!node->isMaster())
            verbose=0;
    }

    void readParams()
    {
        fn_commands = getParam("-i");
        if (node->size < 2)
            REPORT_ERROR(ERR_ARG_INCORRECT,
                         "This program cannot be executed in a single working node");
    }

    void defineParams()
    {
    	addUsageLine("Run commands in a text file in a parallel environment");
    	addUsageLine("+You may use the tag MPI_BARRIER to separate execution blocks.");
    	addUsageLine("+You may use the tag MPI_NEWLINE to force a line break (this is useful for programs accepting parameters directly from stdin).");
    	addParamsLine("-i <commandFile>    : File with commands in different lines");
    }

    void show()
    {
    	if (!verbose)
    		return;
        std::cout << "Commands  file: " << fn_commands << std::endl;
    }

    /* Run --------------------------------------------------------------------- */
#define MAX_LINE 2048
    char szline[MAX_LINE+1];
    void run()
    {
        if (node->rank == 0)
        {
            std::ifstream fh_in;
            fh_in.open(fn_commands.c_str());
            if (!fh_in)
                REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"Cannot open " + fn_commands);
            std::string line;
            int number_of_node_waiting = 0; // max is nprocs -1
            while (!fh_in.eof())
            {
                //wait until a server is free
                MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0,
                         MPI_COMM_WORLD, &status);
                number_of_node_waiting++;
                getline(fh_in, line);
                line=findAndReplace(line,"MPI_NEWLINE","\n");
                strcpy(szline, line.c_str());

                std::string::size_type loc = line.find("MPI_BARRIER", 0);
                if (loc != std::string::npos)
                {
                    while (number_of_node_waiting < (node->size - 1))
                    {
                        MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0,
                                 MPI_COMM_WORLD, &status);
                        number_of_node_waiting++;
                    }
                    while (number_of_node_waiting > 0)
                    {
                        MPI_Send(&szline, 1, MPI_CHAR, number_of_node_waiting,
                                 TAG_WAIT, MPI_COMM_WORLD);
                        number_of_node_waiting--;
                    }
                    continue;
                }

                //send work
                MPI_Send(&szline, MAX_LINE, MPI_CHAR, status.MPI_SOURCE,
                         TAG_WORK, MPI_COMM_WORLD);
                number_of_node_waiting--;
            }

            fh_in.close();
            for (size_t i = 1; i < node->size; i++)
                MPI_Send(0, 0, MPI_INT, i, TAG_STOP, MPI_COMM_WORLD);
        }
        else
        {
            while (1)
            {
                //I am free
                MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
                //get your next task
                MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == TAG_STOP)//I am free
                    break;
                else if (status.MPI_TAG == TAG_WAIT)//wait
                {
                    MPI_Recv(&szline, 1, MPI_CHAR, 0, TAG_WAIT, MPI_COMM_WORLD, &status);
                    continue;
                }
                else if (status.MPI_TAG == TAG_WORK)//work to do
                {
                    MPI_Recv(&szline, MAX_LINE, MPI_CHAR, 0, TAG_WORK, MPI_COMM_WORLD, &status);
                    //do the job
                    if(strlen(szline)<1)
                        continue;
                    else
                    {
                        if (system(szline)==-1)
                        	REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
                    }
                }
                else
                    std::cerr << "WRONG TAG RECEIVED" << std::endl;

            }
        }

    }
};

int main(int argc, char *argv[])
{
	ProgMPIRun prm(argc, argv);
    prm.read(argc, argv);
    int result = prm.tryRun();
    MPI_Finalize();
    return result;
}
