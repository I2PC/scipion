/***************************************************************************
 *
 * Authors:   Sjors H.W. Scheres
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

#include <data/args.h>
#include <data/progs.h>
#include <data/image.h>
#include <cstdio>

class Image_convert_parameters: public Prog_parameters
{
public:
	bool verbose;
	bool output_is_stack, input_is_stack;
    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        verbose = checkParameter(argc, argv, "-verbose");
        // FIXME not yet implemented...
        output_is_stack = checkParameter(argc, argv, "-output_stack");
        input_is_stack = checkParameter(argc, argv, "-input_stack");
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-verbose]                : Verbose output\n";
    }


};


bool process_img(Image<double> &img, const Prog_parameters *prm)
{
	Image_convert_parameters *eprm = (Image_convert_parameters *) prm;
	if (eprm->verbose)
	{
		img().printShape();
	}
	return true;
}

int main(int argc, char **argv)
{
    Image_convert_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img);
}
