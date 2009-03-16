/***************************************************************************
 *
 * Authors:     Javier Rodr�guez Falces
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

#include <reconstruction/crystal_create_surface.h>
#include <data/error.h>

#include <iostream>

/* ------------------------------------------------------------------------- */
/* Program                                                                   */
/* ------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    Prog_create_surface prm;
    // FileName filein;
    try
    {
        prm.read(argc, argv);
        prm.read_input_file();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(1);
    }

    try
    {
        prm.run();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }

    exit(0);
}
