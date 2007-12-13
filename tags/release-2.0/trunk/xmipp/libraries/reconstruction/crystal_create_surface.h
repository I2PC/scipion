/***************************************************************************
 *
 * Authors:     Javier Rodriguez Falces (jrodriguez@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef CRYSTAL_CREATE_SURFACE_H
#define CRYSTAL_CREATE_SURFACE_H

#include <data/docfile.h>
#include <data/matrix1d.h>
#include <data/integration.h>

#include <iostream>

/**@defgroup CrystalSurface crystal_create_surface (Create a surface from a crystal volume)
   @ingroup ReconsLibraryPrograms */
//@{
/* Create_Surface Program Parameters ------------------------------------------ */
/** Please go to this URL
 http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/CrystalSurface
for further details
*/

class Prog_create_surface
{
public:
    /** Struct to keep the data of the plane lattice and the
    mapped lattice */
    struct surface_coordinates
    {
        double cox_ideal;
        double coy_ideal;
        double coz_ideal;

        double cox_real;
        double coy_real;
        double coz_real;

        double Nx;
        double Ny;
        double Nz;
    };

    /** Input file */
    FileName fn_in;

    /** Output file */
    FileName fn_out;

    /** Paramter that determines the shape of the surface to
    be mapped */
    string option;

    /** parameters included in the input file */
    int hmax, hmin, kmin, kmax;
    Matrix1D<double > a, b;
    double x_des, y_des, z_des, x_nor_des, y_nor_des, z_nor_des;

public:
    /** Constructor */
    Prog_create_surface();

    /** Read from a command line.  The program is expecting an
    input file (fn_in), an output file (fn_out) and the type of
    surface to map (option).  An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Read lattice parameters from an input file (fn_in) */
    void read_input_file();

    /** Usage message.
        This function shows the way of introducing this parameters. */
    void usage();

    /** The maping function maps each (h,k) point into a parabole or cosine surface */
    void maping_function(Matrix1D<double > a, Matrix1D<double > b, int h, int k,  surface_coordinates &result);

    /** Run function. This function carries out the maping
    of the points of the plane lattice into other shape
    lattice   */
    void run();
};
//@}
#endif
