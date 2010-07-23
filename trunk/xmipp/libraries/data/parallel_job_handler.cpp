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


#include "parallel_job_handler.h"

ParallelJobHandler::ParallelJobHandler(int nJobs, int bSize, char *fName)
{
    numberOfJobs = nJobs;
    blockSize = bSize;
    assignedJobs = 0;

    if (fName[0] == '\0')
    {
        strcpy(lockFilename, "pijol_XXXXXX");
        if ((lockFile = mkstemp(lockFilename)) == -1)
        {
            perror("ParallelJobHandler::Error generating tmp lock file");
            exit(1);
        }
        else
            close(lockFile);
        strcpy(fName, lockFilename);
        std::cerr << "lockFilename: " << lockFilename << std::endl;
    }
    else
        strcpy(lockFilename, fName);
    createLockFile();
}

ParallelJobHandler::ParallelJobHandler(const char *fName)
{
    if (fName[0] == '\0')
        REPORT_ERROR(-33, "ParallelJobHandler: lock filename couldn't be empty");

    strcpy(lockFilename, fName);
    loadLockFile();
}

ParallelJobHandler::~ParallelJobHandler()
{
    close(lockFile);
    if (fileCreator && remove(lockFilename) == -1)
        perror("ParalleJobHandler: error deleting lock file");
}

void ParallelJobHandler::createLockFile()
{

    int bytes;
    int buffer[] = {numberOfJobs, assignedJobs, blockSize};

    if (( lockFile = open(lockFilename, O_CREAT | O_RDWR | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) ) == -1)
        //if ((lockFile = mkstemp(lockFilename)) == -1)
    {
        perror("ParallelJobHandler::createLockFile: Error opening lock file");
        exit(1);
    }

    if (write(lockFile, buffer, 3*sizeof(int)) == -1)
    {
        perror("ParallelJobHandler::createLockFile: Error writing to lock file");
        exit(1);
    }

    writeVars();
    fileCreator = true;
}//function createLockFile

void ParallelJobHandler::loadLockFile()
{
    if (( lockFile = open(lockFilename, O_RDWR) ) == -1)
    {
        perror("ParallelJobHandler::loadLockFile: Error opening lock file");
        exit(1);
    }
    readVars();
    fileCreator = false;
}

void ParallelJobHandler::readVars()
{
    lseek(lockFile, 0, SEEK_SET);
    read(lockFile, &numberOfJobs, sizeof(int));
    read(lockFile, &assignedJobs, sizeof(int));
    read(lockFile, &blockSize, sizeof(int));
}

void ParallelJobHandler::writeVars()
{
    lseek(lockFile, 0, SEEK_SET);
    write(lockFile, &numberOfJobs, sizeof(int));
    write(lockFile, &assignedJobs, sizeof(int));
    write(lockFile, &blockSize, sizeof(int));
}

void ParallelJobHandler::lock()
{
    lseek(lockFile, 0, SEEK_SET);
    lockf(lockFile, F_LOCK, 0);
}

void ParallelJobHandler::unlock()
{
    lseek(lockFile, 0, SEEK_SET);
    lockf(lockFile, F_ULOCK, 0);
}

bool ParallelJobHandler::getJobs(int &first, int &last)
{
    bool ret = true;
    first = last = -1;

    lock();
    readVars();

    if (assignedJobs >= numberOfJobs)
    {
        ret = false;
    }
    else
    {
        first = assignedJobs;
        assignedJobs = XMIPP_MIN(assignedJobs + blockSize, numberOfJobs);
        last = assignedJobs - 1;
        writeVars();
    }

    unlock();

    return ret;
}
