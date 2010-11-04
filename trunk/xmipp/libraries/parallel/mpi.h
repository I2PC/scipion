/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#ifndef XMIPP_MPI_H_
#define XMIPP_MPI_H_

#include <mpi.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "data/threads.h"

/** @addtogroup Parallel
 * @{
 */

/** Class to wrapp some MPI common calls in an work node.
 *
 */
class MpiNode
{

public:
    int rank, size;
    MpiNode(int &argc, char ** argv);
    ~MpiNode();

    bool isMaster() const;
    /** Wait on a barrier for the other MPI nodes */
    void barrierWait();
};

//mpi macros
#define TAG_WORK   0
#define TAG_STOP   1
#define TAG_WAIT   2

/** This class is another implementation of ParallelTaskDistributor with MPI workers.
 * It extends from ThreadTaskDistributor and add the MPI call
 * for making the distribution and extra locking mechanims between
 * MPI nodes.
 */
class MpiTaskDistributor: public ThreadTaskDistributor
{
protected:
    MpiNode * node;
    ThreadManager * manager;

    virtual bool distribute(longint &first, longint &last);

public:
    MpiTaskDistributor(longint nTasks, longint bSize, MpiNode *node);
    ~MpiTaskDistributor();

    friend void __threadMpiMasterDistributor(ThreadArgument &arg);
};//end of class MpiTaskDistributor

/** Another implementation of ParallelTaskDistributor using a file as lock mechanism.
 * It will extends from ThreadTaskDistributor for also be compatible with several
 * threads running in the same process and also syncronization between different
 * process since only one will get the lock on the file.
 */
class FileTaskDistributor: public ThreadTaskDistributor
{
private:
    void createLockFile();
    void loadLockFile();
    void readVars();
    void writeVars();

protected:
    MpiNode * node;
    int lockFile;
    char lockFilename[L_tmpnam];
    bool fileCreator;

    virtual void lock();
    virtual void unlock();

public:
    FileTaskDistributor(longint nTasks, longint bSize, MpiNode * node  = NULL);
    virtual ~FileTaskDistributor();
};//end of class FileTaskDistributor

/** @} */
#endif /* XMIPP_MPI_H_ */
