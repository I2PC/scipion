/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 *     Joaquin Oton         (joton@cnb.csic.es)
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

#ifndef MPI_RECONSTRUCT_WBP_H_
#define MPI_RECONSTRUCT_WBP_H_

#include "xmipp_mpi.h"
#include <reconstruction/reconstruct_wbp.h>

class ProgMPIRecWbp: public ProgRecWbp, public MpiMetadataProgram
{
public:
    //Empty constructor
    ProgMPIRecWbp()
    {
    }
    ProgMPIRecWbp(int argc, char **argv);
    ProgMPIRecWbp(MpiNode *node);
    void defineParams();
    void readParams();
    void read(int argc, char **argv);
    void produceSideInfo();
    void showProgress();
    bool getImageToProcess(size_t &objId, size_t &objIndex);
    void finishProcessing();
};


#endif /* MPI_RECONSTRUCT_WBP_H_ */
