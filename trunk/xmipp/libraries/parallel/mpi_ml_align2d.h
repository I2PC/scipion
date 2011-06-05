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
#include "reconstruction/ml_refine3d.h"

/**@defgroup MPI_Programs Programs that parallelize using MPI library
   @ingroup ParallelLibrary */
//@{

/** Class to organize some useful MPI-functions for ML programs
 * It will also serve as base for those programs*/
class MpiML2DBase
{
protected:
    MpiNode *node;
    bool created_node;
    //Reference to the program to be parallelized
    XmippProgram * program;

public:
    /** Read arguments sequentially to avoid concurrency problems */
    void readMpi(int argc, char **argv);
    /** Default constructor */
    MpiML2DBase(XmippProgram * prm);
    /** Constructor passing the MpiNode */
    MpiML2DBase(XmippProgram * prm, MpiNode * node);
    /** Destructor */
    ~MpiML2DBase();
    /** This function is only valid for 2D ML programs*/
    void sendDocfile(const MultidimArray<double> &data);
}
;//end of class MpiML

/** Class to parallelize the ML 2D alignment program */
class MpiProgML2D: public ProgML2D, public MpiML2DBase
{

public:
    /** Default constructor */
    MpiProgML2D();
    /** Constructor passing the MpiNode */
    MpiProgML2D(MpiNode * node);
    /** Only take a part of images for process */
    void setNumberOfLocalImages();
    /** All mpi nodes should syncronize at this point
     * that's why the need of override the implementation.
     */
    void produceSideInfo2();
    /// Write model parameters
    void writeOutputFiles(const ModelML2D &model, OutputType outputType = OUT_FINAL);
    /// After normal ML2D expectation, data must be collected from nodes
    void expectation();
    /// Redefine endIteration for some syncronization
    void endIteration();
    //Just for debugging
    void printModel(const String &msg, const ModelML2D & model);
    //Redefine usage, only master should print
    virtual void usage(int verb = 0) const;

}
;//end of class MpiProgML2D

/** Class to parallelize the ML 3D refinement program */
class MpiProgMLRefine3D: public ProgMLRefine3D, public MpiML2DBase
{
public:
    /** Constructor */
    MpiProgMLRefine3D(int argc, char ** argv, bool fourier = false);
    /** Destructor */
    virtual ~MpiProgMLRefine3D();
    /** Only master copy reference volumes before start processing */
    void copyVolumes();
    /** Only master postprocess volumes */
    void postProcessVolumes();
    /** Project volumes, sync after projection */
    void projectVolumes(MetaData &mdProj);
    /** Make noise images, only master */
    void makeNoiseImages();
    /// Calculate 3D SSNR, only master and broadcast result
    void calculate3DSSNR(MultidimArray<double> &spectral_signal);
    /// Convergency check, only master and broadcast result
    bool checkConvergence() ;


}
;//end of class  MpiProgMLRefine3D

/** Class to paralleliza the MLF 2D alignment program */
class MpiProgMLF2D: public ProgMLF2D, public MpiML2DBase
{
public:
  /** Default constructor */
  MpiProgMLF2D();
  /** Constructor passing the MpiNode */
  MpiProgMLF2D(MpiNode * node);
  /** All mpi nodes should syncronize at this point
   * that's why the need of override the implementation.
   */
  void produceSideInfo();
  void produceSideInfo2();
  /// Write model parameters
  void writeOutputFiles(const ModelML2D &model, OutputType outputType = OUT_FINAL);
  /// After normal ML2D expectation, data must be collected from nodes
  void expectation();
  /// Redefine endIteration for some syncronization
  void endIteration();
}
;//end of class MpiProgMLF2D


#endif /* MPI_ML_ALIGN2D_H_ */
