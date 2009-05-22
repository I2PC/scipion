/***************************************************************************
 *
 * Authors:     Slavica JONIC (slavica.jonic@impmc.jussieu.fr, slavica.jonic@a3.epfl.ch)
 *		Jean-Noel PIOCHE (jnp95@hotmail.com) 
 *		
 * Biomedical Imaging Group, EPFL (Lausanne, Suisse).
 * Structures des Assemblages Macromoléculaires, IMPMC UMR 7590 (Paris, France).
 * IUT de Reims-Châlons-Charleville (Reims, France).
 *
 * Last modifications by JNP the 12/05/2009 13:22:12 
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

#include <reconstruction/projection_real_shears.h>

int main(int argc, char *argv[])
{
    Prog_Project_Parameters_2      prog_prm;
    Projection                     proj;
    SelFile                        SF;

// Check the command line
    try
    {
        prog_prm.read(argc, argv);
    }
    catch (Xmipp_error &XE)
    {
        std::cout << XE;
        prog_prm.usage();
        exit(1);
    }

// Really project
    try
    {
        ROUT_project_real_shears(prog_prm, proj, SF);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    return 0;
}

