/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.csic.es)
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

#include <reconstruction/break_symmetry.h>

int main(int argc, char **argv)
{

    int                         iter;
    std::vector<SelFile>             SFout;
    double                      avecorr;
    FileName                    fn_tmp;
    Prog_Break_Sym_prm          prm;

    // Get input parameters
    try
    {

        // Read command line
        prm.read(argc, argv);
        prm.show();

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

    try
    {

        // Loop over all iterations
        iter = prm.istart;
        while (iter <= prm.Niter)
        {

            if (prm.verb > 0)
            {
                std::cerr        << "--> Break symmetry:  iteration " << iter << " of " << prm.Niter << std::endl;
                prm.fh_hist << "--> Break symmetry:  iteration " << iter << " of " << prm.Niter << std::endl;
            }

            prm.process_selfile(prm.SF, SFout, avecorr);

            if (prm.verb > 0)
            {
                std::cerr        << "--> Average ccf: " << avecorr << std::endl;
                prm.fh_hist << "--> Average ccf: " << avecorr << std::endl;
            }

            for (int volno = 0; volno < prm.Nvols; volno++)
                prm.reconstruction(argc, argv, SFout[volno], iter, volno);

            iter++;
        } // end loop iterations

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}




