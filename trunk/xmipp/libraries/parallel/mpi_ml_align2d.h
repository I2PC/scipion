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

#ifndef MPI_ML_ALIGN2D_H_
#define MPI_ML_ALIGN2D_H_

#include "parallel/mpi.h"
#include "reconstruction/ml_align2d.h"


/**@defgroup MPI_Programs Programs that parallelize using MPI library
   @ingroup ParallelLibrary */
//@{
/** Program to parallelize the ML 2D alignment program */
class MpiProgML2D: public ProgML2D
{
private:
    MpiNode *node;

public:

    /** Destructor */
    ~MpiProgML2D();
    /** Redefine XmippProgram read to perform some syncronization */
    void read(int argc, char** argv);
    /** Only take a part of images for process */
    void setNumberOfLocalImages();
    /** All mpi nodes should syncronize at this point
     * that's why the need of override the implementation.
     */
    void produceSideInfo2();
    /// Add docfiledata to docfile
    void addPartialDocfileData(const MultidimArray<double> &data, int first, int last);
    /// Write model parameters
    void writeOutputFiles(const ModelML2D &model, int outputType = OUT_FINAL);
    /// After normal ML2D expectation, data must be collected from nodes
    void expectation();

};//end of class MpiProgML2D

#endif /* MPI_ML_ALIGN2D_H_ */
