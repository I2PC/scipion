 /*
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
#include <data/xmipp_threads.h>
#include <data/xmipp_program.h>
#define XMIPP_MPI_SIZE_T MPI_UNSIGNED_LONG
/** @defgroup MPI MPI
 *  @ingroup ParallelLibrary
 * @{
 */

/** Class to wrapp some MPI common calls in an work node.
*
*/
class MpiNode
{
public:

    //MPI_Comm *comm;
    size_t rank, size, active;//, activeNodes;
    MpiNode(int &argc, char ** argv);
    ~MpiNode();

    /** Check if the node is master */
    bool isMaster() const;

    /** Wait on a barrier for the other MPI nodes */
    void barrierWait();

    /** Gather metadatas */
    void gatherMetadatas(MetaData &MD, const FileName &rootName);

    /** Update the MPI communicator to connect the currently active nodes */
//    void updateComm();

protected:
    /** Calculate the number of still active nodes */
    size_t getActiveNodes();

};

//mpi macros
#define TAG_WORK   0
#define TAG_STOP   1
#define TAG_WAIT   2

/** This class is another implementation of ParallelTaskDistributor with MPI workers.
 * It extends from ThreadTaskDistributor and adds the MPI call
 * for making the distribution and extra locking mechanisms among
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
    virtual ~MpiTaskDistributor();

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

/** This class represent an Xmipp MPI Program.
 *  It includes the basic MPI functionalities to the programs,
 *  like an mpinode, a mutex...
 *
 *  To be compatible with inheritance multiple, the  BaseXmippProgram
 *  must be declared with XmippProgramm as virtual.
 *
 * @code
 * class BaseProgram: public virtual XmippProgram {};
 * class BaseMpiProgram: public BaseProgram, public XmippMpiProgram {};
 * @endcode
 */
class XmippMpiProgram: public virtual XmippProgram
{
protected:
    /** Mpi node */
    MpiNode * node;
    bool created_node;
    /** Number of Processors **/
    size_t nProcs;
    /** Number of independent MPI jobs **/
    size_t numberOfJobs;
    /** status after an MPI call */
    MPI_Status status;

    XmippMpiProgram();
    ~XmippMpiProgram();

    /** Provide a node when calling from another MPI program  */
    void setNode(MpiNode * node);

public:
    /** Read MPI params from command line */
    void read(int argc, char **argv);
    /** Call the run function inside a try/catch block
    * sending an abort signal to the rest of mpi nodes.
    * */
    int tryRun();
};

class MpiMetadataProgram: public XmippMpiProgram
{
protected:
    /** Divide the job in this number block with this number of images */
    int blockSize;
    MpiTaskDistributor *distributor;
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

    void defineParams();
    void readParams();
    /** Create task distributor */
    void createTaskDistributor(MetaData &mdIn, size_t blockSize = 0);
    /** Preprocess */
    void preProcess();
    /** finishProcessing */
    void finishProcessing();
    /** Get task to process */
    bool getTaskToProcess(size_t &objId, size_t &objIndex);
};

/** Macro to define a simple MPI parallelization
 * of a program based on XmippMetaDataProgram */
#define CREATE_MPI_METADATA_PROGRAM(baseClassName, mpiClassName) \
class mpiClassName: public baseClassName, public MpiMetadataProgram\
{\
public:\
    void defineParams()\
    {\
        baseClassName::defineParams();\
        MpiMetadataProgram::defineParams();\
    }\
    void readParams()\
    {\
        MpiMetadataProgram::readParams();\
        baseClassName::readParams();\
    }\
    void read(int argc, char **argv, bool reportErrors = true)\
    {\
        MpiMetadataProgram::read(argc,argv);\
        baseClassName::read(argc, argv, reportErrors);\
    }\
    void preProcess()\
    {\
        baseClassName::preProcess();\
        MetaData &mdIn = *getInputMd();\
        mdIn.addItemId(); \
        createTaskDistributor(mdIn, blockSize);\
    }\
    void startProcessing()\
    {\
        if (node->isMaster())\
            baseClassName::startProcessing();\
        node->barrierWait();\
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
        node->gatherMetadatas(*getOutputMd(), fn_out);\
    	MetaData MDaux; \
    	MDaux.sort(*getOutputMd(), MDL_ITEM_ID); \
        MDaux.removeItemId(); \
        *getOutputMd()=MDaux; \
        if (node->isMaster())\
            baseClassName::finishProcessing();\
    }\
};\

/** @} */
#endif /* XMIPP_MPI_H_ */
