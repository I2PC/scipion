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

#include "mpi.h"

/** Function for a thread waiting on master MPI node distributing tasks.
 * This function will be called for one thread from the constructor
 * of the MpiTaskDistributor in the master.
 */
void __threadMpiMasterDistributor(ThreadArgument &arg)
{
    MpiTaskDistributor * distributor = (MpiTaskDistributor*) arg.workClass;
    int size = distributor->node->size;
    size_t workBuffer[3];
    MPI_Status status;
    int finalizedWorkers = 0;

    while (finalizedWorkers < size - 1)
    {
        //wait for request form workers
        MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        workBuffer[0] = distributor->ThreadTaskDistributor::distribute(
                            workBuffer[1], workBuffer[2]) ? 1 : 0;
        if (workBuffer[0] == 0) //no more jobs
            finalizedWorkers++;
        //send work
        MPI_Send(workBuffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, TAG_WORK,
                 MPI_COMM_WORLD);
    }
}

MpiTaskDistributor::MpiTaskDistributor(size_t nTasks, size_t bSize,
                                       MpiNode *node) :
        ThreadTaskDistributor(nTasks, bSize)
{
    this->node = node;
    //if master create distribution thread
    if (node->isMaster())
    {
        manager = new ThreadManager(1, (void*) this);
        manager->runAsync(__threadMpiMasterDistributor);
    }
}

MpiTaskDistributor::~MpiTaskDistributor()
{
    if (node->isMaster())
    {
        manager->wait();
        delete manager;
    }
}

bool MpiTaskDistributor::distribute(size_t &first, size_t &last)
{
    if (node->isMaster())
        return ThreadTaskDistributor::distribute(first, last);

    //If not master comunicate with master thread
    //to get tasks
    //workBuffer[0] = 0 if no more jobs, 1 otherwise
    //workBuffer[1] = first
    //workBuffer[2] = last
    size_t workBuffer[3];
    MPI_Status status;
    //any message from the master, is tag is TAG_STOP then stop
    MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Recv(workBuffer, 3, MPI_LONG_LONG_INT, 0, TAG_WORK, MPI_COMM_WORLD,
             &status);

    first = workBuffer[1];
    last = workBuffer[2];

    return (workBuffer[0] == 1);
}

// ================= FILE MUTEX ==========================
MpiFileMutex::MpiFileMutex(MpiNode * node)
{
    fileCreator = false;

    if (node == NULL || node->isMaster())
    {
        fileCreator = true;
        strcpy(lockFilename, "pijol_XXXXXX");
        if ((lockFile = mkstemp(lockFilename)) == -1)
        {
            perror("MpiFileMutex::Error generating tmp lock file");
            exit(1);
        }
        close(lockFile);
    }
    //if using mpi broadcast the filename from master to slaves
    if (node != NULL)
        MPI_Bcast(lockFilename, L_tmpnam, MPI_CHAR, 0, MPI_COMM_WORLD);

    if ((lockFile = open(lockFilename, O_RDWR)) == -1)
    {
        perror("MpiFileMutex: Error opening lock file");
        exit(1);
    }
}

void MpiFileMutex::lock()
{
    Mutex::lock();
    lseek(lockFile, 0, SEEK_SET);
    lockf(lockFile, F_LOCK, 0);
}

void MpiFileMutex::unlock()
{
    lseek(lockFile, 0, SEEK_SET);
    lockf(lockFile, F_ULOCK, 0);
    Mutex::unlock();
}

MpiFileMutex::~MpiFileMutex()
{
    close(lockFile);
    if (fileCreator && remove(lockFilename) == -1)
    {
        perror("~MpiFileMutex: error deleting lock file");
        exit(1);
    }
}

FileTaskDistributor::FileTaskDistributor(size_t nTasks, size_t bSize,
        MpiNode * node):
        ThreadTaskDistributor(nTasks, bSize)
{
    fileMutex=new MpiFileMutex(node);
    if (node == NULL || node->isMaster())
        createLockFile();
    loadLockFile();
}

FileTaskDistributor::~FileTaskDistributor()
{
    delete fileMutex;
}

void FileTaskDistributor::createLockFile()
{
    int buffer[] = { numberOfTasks, assignedTasks, blockSize };

    if ((lockFile = open(fileMutex->lockFilename, O_CREAT | O_RDWR | O_TRUNC, S_IRUSR
                         | S_IWUSR | S_IRGRP | S_IROTH)) == -1)
    {
        perror("FileTaskDistributor::createLockFile: Error opening lock file");
        exit(1);
    }

    if (write(lockFile, buffer, 3 * sizeof(int)) == -1)
    {
        perror("FileTaskDistributor::createLockFile: Error writing to lock file");
        exit(1);
    }

    writeVars();
}//function createLockFile

void FileTaskDistributor::loadLockFile()
{
    if ((lockFile = open(fileMutex->lockFilename, O_RDWR)) == -1)
    {
        perror("FileTaskDistributor::loadLockFile: Error opening lock file");
        exit(1);
    }
    readVars();
}

void FileTaskDistributor::readVars()
{
    lseek(lockFile, 0, SEEK_SET);
    read(lockFile, &numberOfTasks, sizeof(size_t));
    read(lockFile, &assignedTasks, sizeof(size_t));
    read(lockFile, &blockSize, sizeof(size_t));
}

void FileTaskDistributor::writeVars()
{
    lseek(lockFile, 0, SEEK_SET);
    write(lockFile, &numberOfTasks, sizeof(size_t));
    write(lockFile, &assignedTasks, sizeof(size_t));
    write(lockFile, &blockSize, sizeof(size_t));
}

void FileTaskDistributor::lock()
{
    fileMutex->lock();
    readVars();
}

void FileTaskDistributor::unlock()
{
    writeVars();
    fileMutex->unlock();
}


//------------ MPI ---------------------------
MpiNode::MpiNode(int &argc, char ** argv)
{
    //MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

MpiNode::~MpiNode()
{
    MPI_Finalize();
}

bool MpiNode::isMaster() const
{
    return rank == 0;
}

void MpiNode::barrierWait()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

void MpiNode::gatherMetadatas(MetaData &MD, const FileName &rootname,
                              MDLabel sortLabel)
{
    FileName fn;

    if (!isMaster())//workers just write down partial results
    {
        fn = formatString("%s_node%d.xmd", rootname.c_str(), rank);
        MD.write(fn);
    }
    ///Wait for all workers write results
    barrierWait();
    if (isMaster()) //master should collect and join workers results
    {
        MetaData mdAll, mdSlave;
        fn = formatString("%s_node%d.xmd", rootname.c_str(), nodeRank);
        mdAll.read(fn);
        for (int nodeRank = 2; nodeRank < size; nodeRank++)
        {
            fn = formatString("%s_node%d.xmd", rootname.c_str(), nodeRank);
            mdSlave.read(fn);
            //make sure file is not empty
            if (!mdSlave.isEmpty())
                mdAll.unionAll(mdSlave);
            //Remove blockname
            fn=fn.removeBlockName();
            remove(fn.c_str());
        }
        MD.sort(mdAll, MDL_IMAGE);
    }
}

/* -------------------- MpiMetadataProgram ------------------- */
MpiMetadataProgram::~MpiMetadataProgram()
{
    delete node;
    delete fileMutex;
    delete distributor;
}

void MpiMetadataProgram::read(int argc, char **argv)
{
    node = new MpiNode(argc, argv);
    fileMutex = new MpiFileMutex(node);
    last=0;
    first=1;
}

void MpiMetadataProgram::createTaskDistributor(const MetaData &mdIn)
{
    int blockSize=mdIn.size()/(node->size*5);
    if (blockSize<1)
        blockSize=1;
    mdIn.findObjects(imgsId);
    distributor = new FileTaskDistributor(mdIn.size(), blockSize, node);
}

//Now use the distributor to grasp images
size_t MpiMetadataProgram::getTaskToProcess()
{
    bool moreTasks=true;
    if (first>last)
        moreTasks = distributor->getTasks(first, last);
    if (moreTasks)
        return imgsId[first++];
    first = distributor->numberOfTasks-1;
    return BAD_OBJID;
}

