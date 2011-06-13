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
#include "data/program.h"

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
    /** Gather metadatas */
    void gatherMetadatas(MetaData &MD, const FileName &rootName,
    		MDLabel sortLabel=MDL_IMAGE);
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

    virtual bool distribute(size_t &first, size_t &last);

public:
    MpiTaskDistributor(size_t nTasks, size_t bSize, MpiNode *node);
    ~MpiTaskDistributor();

    friend void __threadMpiMasterDistributor(ThreadArgument &arg);
}
;//end of class MpiTaskDistributor

/** Mutex on files.
 * This class extends threads mutex to also provide file locking.
 */
class MpiFileMutex: public Mutex
{
protected:
	MpiNode * node;
	    int lockFile;
	    bool fileCreator;
public:
    /** Default constructor. */
	MpiFileMutex(MpiNode * node);

    /** Destructor. */
    ~MpiFileMutex();

    /** Function to get the access to the mutex.
     * If the some thread has the mutex and other
     * ask to lock will be waiting until the first one
     * release the mutex
     */
    void lock();

    /** Function to release the mutex.
     * This allow the access to the mutex to other
     * threads that are waiting for it.
     */
    void unlock();

    char lockFilename[L_tmpnam];
}
;//end of class MpiFileMutex
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
    MpiFileMutex *fileMutex;
    int           lockFile;

    virtual void lock();
    virtual void unlock();

public:
    FileTaskDistributor(size_t nTasks, size_t bSize, MpiNode * node  = NULL);
    virtual ~FileTaskDistributor();
}
;//end of class FileTaskDistributor

class MpiMetadataProgram
{
protected:
    MpiNode *node;
    FileTaskDistributor *distributor;
    std::vector<size_t> imgsId;
    MpiFileMutex *fileMutex;
    size_t first, last;

public:
    /** Constructor */
    MpiMetadataProgram();
    /** Destructor */
    ~MpiMetadataProgram();
    /** Read arguments */
    void read(int argc, char **argv);
    /** Create task distributor */
    void createTaskDistributor(const MetaData &mdIn, int blockSize = 0);
    /** Get task to process */
    bool getTaskToProcess(size_t &objId, size_t &objIndex);
};

/** Macro to define a simple MPI paralelization
 * of a program based on XmippMetaDataProgram */
#define CREATE_MPI_METADATA_PROGRAM(baseClassName, mpiClassName) \
class mpiClassName: public baseClassName, public MpiMetadataProgram\
{\
  int blockSize;\
public:\
    void defineParams()\
    {\
        baseClassName::defineParams();\
        addParamsLine("== MPI ==");\
        addParamsLine(" [--mpi_job_size <size=0>]     : Number of images sent simultaneously to a mpi node");\
    }\
    void readParams()\
    {\
        blockSize = getIntParam("--mpi_job_size");\
    }\
    void read(int argc, char **argv)\
    {\
        MpiMetadataProgram::read(argc,argv);\
        if (!node->isMaster())\
            verbose=0;\
        baseClassName::read(argc, argv);\
    }\
    void preProcess()\
    {\
        baseClassName::preProcess();\
        createTaskDistributor(mdIn, blockSize);\
    }\
    void startProcessing()\
    {\
        if (node->isMaster())\
            baseClassName::startProcessing();\
    }\
    void showProgress()\
    {\
        if (node->isMaster())\
        {\
            time_bar_done=first+1;\
            baseClassName::showProgress();\
        }\
    }\
    bool getImageToProcess(size_t &objId, size_t &objIndex)\
    {\
        return getTaskToProcess(objId, objIndex);\
    }\
    void finishProcessing()\
    {\
        node->gatherMetadatas(mdOut, fn_out);\
        if (node->isMaster())\
            baseClassName::finishProcessing();\
    }\
}\

/** @} */
#endif /* XMIPP_MPI_H_ */
