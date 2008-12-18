/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include <reconstruction/phantom_create_micrograph.h>

int main(int argc, char **argv)
{
    Prog_Phantom_Create_Micrograph_Parameters prm;
    
    try {
        prm.read(argc,argv);
    } catch (Xmipp_error XE) {
        std::cout << XE << std::endl;
        prm.usage();
        exit(1);
    }
    
    try {
        prm.show();
        prm.produce_side_info();
        prm.run();
    } catch (Xmipp_error XE) {
        std::cout << XE << std::endl;
        exit(1);
    }
    
    return 0;
}
