/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.uam.es)
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

#include "ctf_wiener_filter.h"

/* Read parameters from command line. -------------------------------------- */
void WienerFilterParams::read(int argc, char **argv)
{
    fn_sel = getParameter(argc, argv, "-i");
    fn_ctfdat = getParameter(argc, argv, "-ctfdat");
    fn_root   = getParameter(argc, argv, "-o","wien");
    phase_flipped = checkParameter(argc, argv, "-phase_flipped");
    wiener_constant = textToFloat(getParameter(argc, argv, "-wiener_constant","10."));

}

/* Show -------------------------------------------------------------------- */
void WienerFilterParams::show()
{

    std::cerr << "  Input selection file    : "<< fn_sel << std::endl;
    std::cerr << "  Input ctfdat file       : "<< fn_ctfdat << std::endl;
    std::cerr << "  Output rootname         : "<< fn_root << std::endl;
    if (phase_flipped)
    {
	std::cerr << " -> Output Wiener filters for data that are PHASE FLIPPED"<<std::endl;
    }
    else
    {
	std::cerr << " -> Output Wiener filters for data that are NOT PHASE FLIPPED"<<std::endl;
    }
}

/* Usage ------------------------------------------------------------------- */
void WienerFilterParams::usage()
{
    std::cerr << "   -i <selfile>                        : Input selection file\n"
              << "   -ctfdat <ctfdat file>               : Input CTFdat file for all data\n"
              << "  [-o <oext=\"wien\">]                   : Output root name\n"
	      << "  [-phase_flipped]                     : Output filters for phase-flipped data\n"
              << "  [-wiener_constant <float=10>]        : Wiener constant to be added to denominator\n"
    ;
}

/* Produce Side information ------------------------------------------------ */
void WienerFilterParams::produceSideInfo()
{

    FileName fnt_img, fnt_ctf, fnt;
    ImageXmipp img;
    XmippCTF ctf;
    Matrix2D<double> Mctf;
    Matrix2D<std::complex<double> >  ctfmask;
    SelFile SF;
    CTFDat ctfdat;
    int dim, ydim, iifocus;
    bool is_unique, found;

    // Read input selfile and get image dimensions
    SF.read(fn_sel);
    SF.ImgSize(dim,ydim);
    if ( dim != ydim )
	REPORT_ERROR(1,"ctf_wiener_filter ERROR%% Only squared images are allowed!");
    Mctf.resize(dim,dim);

    // Read ctfdat and determine the number of CTF groups
    all_fn_ctfs.clear();
    ctfdat.read(fn_ctfdat);
    ctfdat.goFirstLine();
    while (!ctfdat.eof())
    {
	ctfdat.getCurrentLine(fnt_img,fnt_ctf);
	is_unique=true;
	for (int i = 0; i< all_fn_ctfs.size(); i++)
	{
	    if (fnt_ctf == all_fn_ctfs[i])
	    {
		is_unique=false;
		break;
	    }
	}
	if (is_unique)
	{
	    all_fn_ctfs.push_back(fnt_ctf);
	    count_defocus.push_back(0);
		
	    // Read CTF in memory
	    ctf.read(fnt_ctf);
	    ctf.enable_CTF = true;
	    ctf.Produce_Side_Info();
	    ctf.Generate_CTF(dim, dim, ctfmask);
	    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mctf)
	    {
		if (phase_flipped) dMij(Mctf, i, j) = fabs(dMij(ctfmask, i, j).real());
		else dMij(Mctf, i, j) = dMij(ctfmask, i, j).real();
	    }
	    all_Mctfs.push_back(Mctf);
	}
	ctfdat.nextLine();
    }
	
    // Fill count_defocus vectors
    SF.go_beginning();
    while (!SF.eof())
    {
	fnt=SF.NextImg();
	// Find which CTF group it belongs to
	found = false;
	ctfdat.goFirstLine();
	while (!ctfdat.eof())
	{
	    ctfdat.getCurrentLine(fnt_img,fnt_ctf);
	    if (fnt_img == fnt)
	    {
		found = true;
		for (iifocus = 0; iifocus< all_fn_ctfs.size(); iifocus++)
		    if (fnt_ctf == all_fn_ctfs[iifocus])
			break;
		break;
	    }
	    ctfdat.nextLine();
	}
	if (!found)
	    REPORT_ERROR(1, "ctf_wiener_filter ERROR%% Did not find image "+fnt+" in the CTFdat file");
	count_defocus[iifocus]++;
    }

    // Precalculate denominator term of Wiener filter
    //
    //                               C_{ictf}
    //  W_{ictf} =  --------------------------------------------
    //                 wiener_constant + Sum_{i=1}^N C^2_{ictf} 
    //
    Msum.resize(dim,dim);
    Msum.init_constant(wiener_constant);
    for (int ictf = 0; ictf< all_fn_ctfs.size(); ictf++)
    {
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Msum)
	{
	    dMij(Msum,i,j) += count_defocus[ictf] * dMij(all_Mctfs[ictf],i,j) * dMij(all_Mctfs[ictf],i,j); 
	}
    }

}

// Do the actual work
void WienerFilterParams::run()
{

    FourierImageXmipp wien;
    FileName fn_out;
    double aux;

    wien().resize(dim,dim);

    // Loop over all groups
    for (int ictf = 0; ictf< all_fn_ctfs.size(); ictf++)
    {
	// Calculate the actual Wiener filter
	FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(wien())
	{
	    aux = dMij(all_Mctfs[ictf],i,j) / dMij(Msum,i,j);
	    dMij(wien(),i,j) = (std::complex<double>)(aux,0.);
	}

	// Write to disc
	fn_out.compose(fn_root, ictf+1,".fft"); 
	wien.write(fn_out);

    }

}




