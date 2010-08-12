/***************************************************************************
 *
 * Authors:     Joaquin Oton (jotons@cnb.csic.es)
 *              J. M. de la Rosa (jmdelarosa@cnb.csic.es)
 *              R. Marabini (roberto@cnb.csic.es)
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

#ifndef PARALLELJOBHANDLER_H_
#define PARALLELJOBHANDLER_H_

#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include "funcs.h"

/** @defgroup ParallelLibrary Parallel Library
 *  @ingroup DataLibrary
 */

/** This class distributes dynamically N tasks between parallel workers.
 * @ingroup ParallelLibrary
 * This class is a generalization of a common task in a parallel
 * environment of dynamically distribute N tasks between workers(threads or mpi proccess).
 * Each worker will ask for a group of tasks, proccess it and ask for more tasks
 * until there is not more task to process.
 *
 * A lock is required for access the number of unprocessed tasks each time
 * a worker ask for job, this is implemented now by a filesystem lock in a file.(the lockfile)
 */
class ParallelJobHandler
{
private:
    //The lock file handler and temporaly file name
    int lockFile;
    char lockFilename[L_tmpnam];
    //Flag to know who created the file and should deleted
    bool fileCreator;
    //The total number of tasks to be distributed
    long long int numberOfJobs;
    //How many tasks give in each request
    long long int blockSize;
    //The number of tasks that have been assigned
    long long int assignedJobs;

public:
    /** @name Constructors Constructors
     *
     * The are two constructors, one that should be called by the master
     * in wich the lock file is created and parameters should be supplied
     * and the constructor for the slaves, which read read the parameters
     * from the lockfile
     * @{
     */

    /** Constructor for Master node.
     */
    ParallelJobHandler(long long int nJobs, long long int bSize, char *fName = NULL);

    /** Constructor for Slaves nodes.
     */
    ParallelJobHandler(const char *fName);

    /** Destructor.
     */
    ~ParallelJobHandler();
    /** @} */

    /** Set the number of tasks assigned in each request */
    bool setBlockSize(long long int blockSize);
    /** Return the number of tasks assigned in each request */
    int getBlockSize() const;

    /** Gets parallel tasks.
     *  @ingroup ParallelJobHandler
     *  This function will be called by workers for asking tasks
     *  until there are not more tasks to process.
     *  Example:
     *  @code
     *  //...
     *  ParallelJobHandler jobHand;
     *  //...
     *  //function to perform some operation
     *  //to N images executed in parellel     *
     *  void processImages()
     *  {
     *      long long int firstImage, lastImage;
     *      while (jobHand->getJobs(firstImage, lastImage))
     *          for (int image = firstImage; image <= lastImage; ++image)
     *          {
     *              //...
     *              processImage(image);
     *              //...
     *          }
     *  }     *
     *  @endcode
     */
    bool getJobs(long long int &first, long long int &last); // False = no more jobs, true = more jobs


private:
    //This function should be only called in the master
    //for create the lock file
    void createLockFile();
    //Load the file if its already created
    //usually this constructor should be called by slaves
    void loadLockFile();

    //Utils functions to read and write from file
    void readVars();
    void writeVars();

    //Functions for lock and unlock using the lockFile
    void lock();
    void unlock();

}
;//class ParallelJobHandler

/* Measure time differences, the function predcision is miliseconds
 *
 */
#include <sys/time.h>
class elapsedTime
{
    struct timeval start_time, end_time;
public:
    /* Return diference in seconds between two times set by setStart and setEndTime respectively
     *
     */

    double getElapsedTime() const
    {
        return  (double)  (end_time.tv_sec  - start_time.tv_sec ) +
                ((double) (end_time.tv_usec - start_time.tv_usec)/1000000.);
    }
    void setStartTime()
    {
        gettimeofday(&start_time, 0);
    }
    void setEndTime()
    {
        gettimeofday(&end_time, 0);
    }
    bool saneInterval() const
    {
        if(getElapsedTime() < 1.)//less than one second
            return false;
        return true;
    }

}
;

#endif /* PARALLELJOBHANDLER_H_ */
