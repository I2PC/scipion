/***************************************************************************
 *
 * Authors:     Jose Roman Bilbao (jrbcast@ace.ual.es)
 *    			Roberto Marabini (roberto@cnb.csic.es)
 *    			Vahid Abrishami (vabrishami@cnb.csic.es)
 *    			David Strelak (davidstrelak@gmail.com)
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
#ifndef MPI_RECONSTRUCT_FOURIER_ACCEL_H_
#define MPI_RECONSTRUCT_FOURIER_ACCEL_H_

#include "xmipp_mpi.h"
#include <data/args.h>
#include <reconstruction/reconstruct_fourier_accel.h>
#include <data/projection.h>
#include <cstring>
#include <cstdlib>
#include <data/xmipp_funcs.h>
#include <data/matrix2d.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#define TAG_WORKFORWORKER    0
#define TAG_STOP     1
#define TAG_TRANSFER   2
#define TAG_FREEWORKER    3
#define TAG_SETVERBOSE  5

//TODO (MARIANA) Please give more documentation and in a good structure e.g. @name

/**@defgroup ProgMPIRecFourier ProgMPIRecFourier
   @ingroup Programs */
//@{

class ProgMPIRecFourierAccel: public ProgRecFourierAccel, public XmippMpiProgram
{
public:
    /** Empty constructor */
	ProgMPIRecFourierAccel(){};

    /*  constructor ------------------------------------------------------- */
	ProgMPIRecFourierAccel(int argc, char *argv[]);

    /* constructor providing an MpiNode
     * this is useful for using this programs from others
     */
	ProgMPIRecFourierAccel(MpiNode * node);

    /** Special way of reading to sync all nodes */
    void read(int argc, char** argv);

private:
    /** Divide the job in this number block with this number of images */
	int mpi_job_size;

	/** Read parameters from command line */
	void readParams();

	/** Specify supported command line arguments */
	void defineParams();

	/** Pre Run PreRun for all nodes but not for all works */
	void preRun();

	/** Main processing method */
	void run();


};
//@}
//end of class MPI reconstruct fourier

#endif /* MPI_RECONSTRUCT_FOURIER_ACCEL_H_ */
