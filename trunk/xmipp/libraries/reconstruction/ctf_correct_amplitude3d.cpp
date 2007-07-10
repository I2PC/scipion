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

#include "ctf_correct_amplitude3d.h"

#include <data/args.h>

/* Read parameters from command line. -------------------------------------- */
void CorrectAmplitude3DParams::read(int argc, char **argv)
{
    fnCtfdat  = getParameter(argc, argv, "-ctfdat");
    fnOut     = getParameter(argc, argv, "-o","wiener");
    fnNrImgs  = getParameter(argc, argv, "-nr_imgs");
    minResol  = textToFloat(getParameter(argc, argv, "-minres", "-1"));
    fnEnvs    = getParameter(argc, argv, "-env","");
    wienConst = textToFloat(getParameter(argc, argv, "-wc", "0.01"));
    isFlipped = checkParameter(argc, argv, "-phase_flipped");

}

/* Show -------------------------------------------------------------------- */
void CorrectAmplitude3DParams::show()
{
    std::cout << "CTF datfile:                       " << fnCtfdat << std::endl
              << "Wiener constant:                   " << wienConst << std::endl
              << "Output rootname:                   " << fnOut << std::endl
              << "Phase flipped:                     " << isFlipped << std::endl
	      << "Number of images in:               " << fnNrImgs << std::endl
	;
    if (minResol>0)
    {
	std::cout << "Apply Wiener filter only beyond:   " << minResol << " Angstroms" << std::endl;
    }
    if (fnEnvs!="")
    {
	std::cout << "Envelopes in file:                 " << fnEnvs << std::endl;
    }
    
}

/* Usage ------------------------------------------------------------------- */
void CorrectAmplitude3DParams::usage()
{
    std::cerr << "   -ctfdat <CTF datfile>               : 2-column ASCII file with name of volume and CTF param for each defocus group\n"
	      << "   -nr_imgs <docfile>                  : Docfile with number of images in each defocus group\n"
              << "  [-o \"wiener\"]                     : Output rootname \n"
              << "  [-minres <Ang>]                     : Apply Wiener filter only beyond this resolution (in Angstrom)\n"
              << "  [-phase_flipped]                    : Use this if the maps were reconstructed from phase corrected images \n"
              << "  [-wc <0.05>]                         : Wiener constant to prevent boosting of noise \n"
//              << "  [-env <selfile>]                     : selfile to ASCII files with envelope for each defocus group\n"
    ;
}

/* Produce Side information ------------------------------------------------ */
void CorrectAmplitude3DParams::produceSideInfo()
{
    // Read the CTFdat
    ctfdat.read(fnCtfdat);

    // Get dimensions of the volumes
    ctfdat.goFirstLine();
    FileName fnVol, fnCTF;
    ctfdat.getCurrentLine(fnVol,fnCTF);
    VolumeXmipp V;
    V.read(fnVol);
    V().getDimension(Zdim,Ydim,Xdim);
    V.clear();

    // Read the envelopes
    if (fnEnvs!="")
    {
	SFenv.read(fnEnvs);
	if (SFenv.LineNo()!=ctfdat.lineNo())
	{
	    REPORT_ERROR(1, "CTFdat and envelope selfile have unequal number of entries.");
	}
    }

    // Read the docfile with the number of images per group
    DFimgs.read(fnNrImgs);
    if (DFimgs.LineNo()!=ctfdat.lineNo())
    {
	REPORT_ERROR(1, "CTFdat and docfile with number of images per group have unequal number of entries.");
    }

}


/* Make 3D CTF ------------------------------------------------------------- */
void CorrectAmplitude3DParams::generateCTF1D(const FileName &fnCTF, const double nr_steps,
				       Matrix1D<double> &CTF1D)
{
    // Read the CTF 
    ctf.FilterBand = CTF;
    ctf.ctf.enable_CTFnoise = false;
    ctf.ctf.read(fnCTF);
    ctf.ctf.Produce_Side_Info();

    double maxres = ( 0.5 * sqrt(3.) ) / ctf.ctf.Tm;
    double stepsize = maxres / nr_steps;
    CTF1D.resize(nr_steps);
    double res = 0.;

    res=0.;
    for (int step=0; step < nr_steps; step++) 
    {
	if ( (minResol < 0) || (1./res < minResol) )
	{
	    CTF1D(step)=ctf.ctf.CTF_at(res, 0);
	    if (isFlipped) CTF1D(step)=ABS(CTF1D(step));
	}
	else
	{
	    CTF1D(step)=1.;
	}
	res+=stepsize;
    }

}

/* Make Wiener filters ------------------------------------------------------------- */
void CorrectAmplitude3DParams::generateWienerFilters()
{
    Matrix1D<double> CTF1D, sumterm;
    FileName fn_tmp;
    int nrimgs;
    ofstream  fh;
    double res;
    double tot_nr_imgs = 0;
    double minsum=99.e99;
    double maxsum=0.;
    double minwien=99.e99;
    double maxwien=0.;
    // Oversample the 1D CTF and Wiener filter vectors OVERSAMPLE times
    double nr_steps= CEIL(OVERSAMPLE * 0.5 * sqrt((double)(Zdim*Zdim + Ydim*Ydim + Xdim*Xdim)));

    Vctfs1D.clear();
    Vwien1D.clear();
    int ii = 0;
    ctfdat.goFirstLine();
    while (!ctfdat.eof())
    {
	// Calculate 1D CTF
        FileName fnVol, fnCTF;
	ctfdat.getCurrentLine(fnVol,fnCTF);
	generateCTF1D(fnCTF,nr_steps,CTF1D);
	Vctfs1D.push_back(CTF1D);

	// Get the number of images contributing to this group
	DFimgs.search(ii+1);
	tot_nr_imgs += DFimgs(0);

	// Calculate denominator of the Wiener filter
	if (ii==0) 
	{
	    sumterm.resize(CTF1D);
	}
	FOR_ALL_ELEMENTS_IN_MATRIX1D(sumterm) 
	{
	    sumterm(i) += DFimgs(0) * CTF1D(i) * CTF1D(i);
	}
	ctfdat.nextLine();
	ii++;
    }

    FOR_ALL_ELEMENTS_IN_MATRIX1D(sumterm) 
    {
	// Find min and max values of the sumterm
	if (sumterm(i)>maxsum) maxsum=sumterm(i);
	if (sumterm(i)<minsum) minsum=sumterm(i);
	// Add (normalized) Wiener filter constant
	sumterm(i) += tot_nr_imgs*wienConst;
    }

    int iimax = ii;
    // Fill the Wiener filter vector
    for (ii=0; ii<iimax; ii++) 
    {
	// Get the number of images contributing to this group
	DFimgs.search(ii+1);

	FOR_ALL_ELEMENTS_IN_MATRIX1D(CTF1D) 
	{
	    CTF1D(i) = DFimgs(0) * Vctfs1D[ii](i) / sumterm(i);
	    if (CTF1D(i)>maxwien) maxwien=CTF1D(i);
	    if (CTF1D(i)<minwien) minwien=CTF1D(i);

	}
	Vwien1D.push_back(CTF1D);

	// Write CTF and Wiener filter curves to disc
	fn_tmp = fnOut + "_wien";
	fn_tmp.compose(fn_tmp, ii+1, "txt");
	fh.open((fn_tmp).c_str(), ios::out);
	if (!fh) REPORT_ERROR(1, (string)"Error: Cannot write file: " + fn_tmp);
	for (int step = 0; step < nr_steps; step++) 
	{
	    res = (step * sqrt(3.) ) / 
		(OVERSAMPLE * sqrt( (double) (Zdim*Zdim + Ydim*Ydim + Xdim*Xdim) ) );
	    fh << res << " " << Vwien1D[ii](step) << " " << Vctfs1D[ii](step) << "\n";
	}
	fh.close();
    }

    // Some output to screen
    cerr <<" ---------------------------------------------------"<<endl;
    cerr <<" + Number of defocus groups      = "<<ctfdat.lineNo()<<endl;
    cerr <<" + Total number of images        = "<<tot_nr_imgs<<endl;
    cerr <<" + Normalized Wiener constant    = "<<tot_nr_imgs*wienConst<<endl;
    cerr <<" + Minimum of sum in denominator = "<<minsum<<endl;
    cerr <<" + Maximum of sum in denominator = "<<maxsum<<endl;
    cerr <<" + Minimum Wiener filter value   = "<<minwien<<endl;
    cerr <<" + Maximum Wiener filter value   = "<<maxwien<<endl;
    cerr <<" ---------------------------------------------------"<<endl;

}
	
/* Make deconvolved volume -------------------------------------------------- */
void CorrectAmplitude3DParams::generateVolumes()
{

    VolumeXmipp V;
    FileName fnVol, fnCTF;
    Matrix3D<complex<double> > fft,fft_out;
    Matrix1D<int>    idx(3);
    Matrix1D<double> freq(3);
    int ires;

    int ii = 0;
    ctfdat.goFirstLine();
    while (!ctfdat.eof())
    {
	ctfdat.getCurrentLine(fnVol,fnCTF);
	V.read(fnVol);
	FourierTransform(V(),fft);
	if (ii == 0)
	{
	    fft_out.resize(fft);
	}
	FOR_ALL_ELEMENTS_IN_MATRIX3D(fft)
	{
	    XX(idx) = j;
	    YY(idx) = i;
	    ZZ(idx) = k;
	    FFT_idx2digfreq(fft, idx, freq);
	    ires= ROUND(OVERSAMPLE*sqrt(XX(freq)*XX(freq)*Xdim*Xdim+YY(freq)*YY(freq)*Ydim*Ydim+ZZ(freq)*ZZ(freq*Zdim*Zdim)));
	    fft_out(k,i,j)+=(Vwien1D[ii])(ires)*fft(k,i,j);
	}
	ctfdat.nextLine();
	ii++;
    }
    
    // Inverse Fourier Transform
    InverseFourierTransform(fft_out,V());
    fnVol=fnOut+"_deconvolved.vol";
    V.write(fnVol);

    // Calculate CTF-affected volumes
    ii = 0;
    ctfdat.goFirstLine();
    while (!ctfdat.eof())
    {
	ctfdat.getCurrentLine(fnVol,fnCTF);
	FOR_ALL_ELEMENTS_IN_MATRIX3D(fft)
	{
	    XX(idx) = j;
	    YY(idx) = i;
	    ZZ(idx) = k;
	    FFT_idx2digfreq(fft, idx, freq);
	    ires= ROUND(OVERSAMPLE*sqrt(XX(freq)*XX(freq)*Xdim*Xdim+YY(freq)*YY(freq)*Ydim*Ydim+ZZ(freq)*ZZ(freq*Zdim*Zdim)));
	    fft(k,i,j)=Vctfs1D[ii](ires)*fft_out(k,i,j);
	}
	ctfdat.nextLine();
	ii++;

	fnVol=fnOut+"_ctffiltered_group";
	fnVol.compose(fnVol,ii,"vol");
	InverseFourierTransform(fft,V());
	V.write(fnVol);
    }

}

/* Correct a set of images ------------------------------------------------- */
void CorrectAmplitude3DParams::run()
{

    // Calculate the filters 
    generateWienerFilters();

    // Calculate the deconvolved and CTF-filtered volumes
    generateVolumes();
    
}
