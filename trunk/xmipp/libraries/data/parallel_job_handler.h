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


/** @defgroup ParallelJobHandler
 *
 * When processing in parallel, this class distributes the work load.
 *
 * It creates a database with all the jobs and assigns them in small blocks
 * under request from the working nodes.
 *
 *The database may be created in memory (threads) or in disk (mpi).
 */


class ParallelJobHandler
{
private:

    char lockFilename[L_tmpnam];
    int lockFile;
    bool fileCreator;

    int numberOfJobs;
    int blockSize;
    int assignedJobs;

public:


    //TODO: todos los casos: inicializar variable blocksize, crear base de datos,
    ParallelJobHandler(int nJobs, int bSize, char *fName);
    ParallelJobHandler(const char *fName);
    ~ParallelJobHandler();


    bool setBlockSize(int blockSize);
    int getBlockSize();

    bool getJobs(int &first, int &last); // False = no more jobs, true = more jobs



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


};//class ParallelJobHandler


#endif /* PARALLELJOBHANDLER_H_ */
