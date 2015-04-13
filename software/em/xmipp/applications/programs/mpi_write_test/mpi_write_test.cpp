/***************************************************************************
 *
 * Authors: J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
 *
 ****************************************************************************/

/****************************************************************************
 * This program is a simple example of using the MPI parallel library
 * in wich some needed MPI calls are encapsulated in the MpiNode class.
 * In the constructor the MPI initialization are done and the finalize
 * in destructor.
 *
 * In this example PI is calculated by sampling discrete
 * points of a quarter of a circle with radius R and aproximating
 * the area by the number of points inside the circle.
 *
 * Also the points to be process(parallel tasks) are dynamically distrubuted
 * between the MPI nodes. Two distributors are available, one base on
 * filesystem locking and another in wich one of the MPI nodes will be
 * distributing the work.
 ***************************************************************************/

#include "parallel/xmipp_mpi.h"

//Some useful macros
#define CREATE_LOG() FILE * _logML = fopen(formatString("nodo%02d.log", node->rank).c_str(), "w+")
#define LOG(msg) do{fprintf(_logML, "%s\t%s\n", getCurrentTimeString(), msg); fflush(_logML); }while(0)
#define CLOSE_LOG() fclose(_logML)

#define IS_MASTER (node->rank == 0)

ParallelTaskDistributor * distributor;
MpiNode * node;
//#define stackSize 1024
#define stackSize 65536
#define Xdim 51
#define Ydim 77

int main(int argc, char **argv)
{
	FileName fnIN="test_delete_me.mrcs";
    //First argument should be the md filename
    //FileName fn(argv[1]);
    MetaData md;
    //MPI Initialization
    node = new MpiNode(argc, argv);
    int rank = node->rank;
    int size = node->size;
    ////CREATE_LOG();
    //create blank file
    if(rank==0)
    {
        unlink(fnIN.c_str());
        createEmptyFile(fnIN, Xdim, Ydim, 1, stackSize);
    }
    Image<double> Iaux(Xdim,Ydim);
    //Be sure all sync here
    ////LOG("waiting on barrier...");
    node->barrierWait();

    //    if (IS_MASTER)
    for (int var = 1; var <= stackSize; var++)
    {
        if(var%size==rank)
        {
            String ss = formatString("%03d@%s", var,fnIN.c_str());
            std::cerr << "ssIN: value" << ss << " " << (double)rank << std::endl;
            Iaux().initConstant((double)rank);
            node->barrierWait();
            Iaux.write(ss);
        }
    }
    node->barrierWait();
    //check results:
    if(rank==0)
    {
    	int errors=0;
        for (int var = 1; var <= stackSize; var++)
        {
            double value = (double) (var%size);
            String ss = formatString("%03d@%s", var,fnIN.c_str());
            Iaux.read(ss.c_str());
            double min,max,std, avg=0.;
            Iaux().computeStats(avg,std,min,max);
            if (
                ABS(avg-value)> 0.00001 ||
                ABS(std-0)> 0.00001 ||
                ABS(min-value)> 0.00001 ||
                ABS(max-value)> 0.00001
            ){
            	errors++;
                std::cerr << "Error in image :" <<  var
                << " avg " << avg
                << " std " << std
                << " min " << min
                << " max " << max
                << std::endl;}
        }
        std::cerr << "errors:" << errors << std::endl;
    }

    ////LOG("FINISHED!!!");
    ////CLOSE_LOG();
    delete node;
    return 0;
}
