/***************************************************************************
 *
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#include <data/psf_xr.h>
#include <data/args.h>

int main(int argc, char **argv)
{
	FileName fnPSF;
	XmippXRPSF psf;

    try
    {
    	fnPSF=getParameter(argc,argv,"-psf");
    }
    catch (Xmipp_error &XE)
    {
    	std::cerr << XE << std::endl;
    	psf.usage();
    	// usage();
    	return 1;
    }

    try
    {
//    	XmippXRPSF psf;
    	psf.read(fnPSF);
    	std::cout << psf << std::endl;
    }
    catch (Xmipp_error &XE)
    {
    	std::cerr << XE << std::endl;
    	return 1;
    }
    return 0;
}
