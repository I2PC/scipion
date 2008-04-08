/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <reconstruction/ctf_wiener_filter.h>

int main(int argc, char **argv)
{
    WienerFilterParams prm;
    try
    {
       prm.read(argc,argv);
       prm.produceSideInfo();
    }
    catch (Xmipp_error XE)
    {
       std::cerr << XE << std::endl;
       prm.usage();
       exit(1);
    }
    
    try
    {
       prm.run();
       std::cerr << " done!" <<std::endl;
    }
    catch (Xmipp_error XE)
    {
       std::cerr << XE << std::endl;
       prm.usage();
       return 1;
    }
    return 0;
}
