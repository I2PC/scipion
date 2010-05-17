/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

#ifndef _VIRUS_VERTEX_HH
#define _VIRUS_VERTEX_HH

#include <data/matrix1d.h>
#include <data/matrix3d.h>
#include <data/docfile.h>
#include <data/selfile.h>
#include <data/args.h>
#include <data/volume.h>
#include <data/image.h>
#include <data/symmetries.h>
#include <data/projection.h>

//@{

class VirusVertex
{
public:
    /** Read parameters from command line.
        */
    void read(int argc, char **argv);

    /** Show. */
    void show();

    /** Usage. */
    void usage();
    /** Main program */
    void run();
    /** load icosahedron vertex and apply symmetry */
    void loadIcosahedronVertex();
    /** get angles either from a document file or from an image */
    void get_angles_for_image(const FileName &fn, double &rot,
    double &tilt, double &psi, double &xoff, double &yoff, double &flip,
    double &weight);
    /** Process one set of angles */
    void processAngles();
    /** Relate vertex and projection matrices   */
    void assignSymmetryMatricesToVertex();
public:
    /** Filenames */
    FileName fn_sel, fn_doc, fn_root, fn_sym;
    /** Docfile with experimental images */
    DocFile DFangles;
    /** Selfile with experimental images */
    SelFile SF;
    /** Virus radius, distance from center to a vertex in pixels */
    double virusRadius;
    /** extract vertex located at a radius greater than min_radius */
    double minVirusRadius;
    /** set verbose mode on */
    bool verbose;
    /* remove vertex closer than  dim*/
    bool removeCloseVertex;
    /** Output file dimensions */
    int dim;
    /** vector with icosahedron vertex */
    std::vector <Matrix1D<double> > vertices_vectors;
    /** Symmetry information */
    SymList  SL;
    /** symmetry group */
    int symmetry;
    // Column numbers in the docfile
    int col_rot, col_tilt, col_psi, col_xoff, col_yoff, col_flip, col_weight;
    // Symmetry Matrices
    std::vector <Matrix2D<double> > R_repository;
    //relation between symmetry matrices and vertex
    Matrix2D<int> symmetryMatrixVertex;//12 vertex, 60 symmetry matrices, 5-fold symmetry
};


//@}
#endif
