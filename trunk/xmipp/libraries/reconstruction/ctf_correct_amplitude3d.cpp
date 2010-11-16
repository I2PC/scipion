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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "ctf_correct_amplitude3d.h"

#include <data/args.h>
#include <data/fft.h>

/* Read parameters from command line. -------------------------------------- */
void CorrectAmplitude3DParams::read(int argc, char **argv)
{
    fnCtfdat  = getParameter(argc, argv, "-ctfdat");
    fnOut     = getParameter(argc, argv, "-o","wiener");
    minResol  = textToFloat(getParameter(argc, argv, "-minres", "-1"));
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
    ;
    if (minResol>0)
        std::cout << "Apply Wiener filter only beyond:   " << minResol << " Angstroms" << std::endl;
}

/* Usage ------------------------------------------------------------------- */
void CorrectAmplitude3DParams::usage()
{
    std::cerr << "-ctfdat <CTF datfile>               : Metadata with the volumes, ctfs, number of images in that group,\n"
    << "                                         and optionally the envelopes\n"
    << "  [-o \"wiener\"]                        : Output rootname \n"
    << "  [-minres <Ang>]                      : Apply Wiener filter only beyond this resolution (in Angstrom)\n"
    << "  [-phase_flipped]                     : Use this if the maps were reconstructed from phase corrected images \n"
    << "  [-wc <0.05>]                         : Wiener constant (to be multiplied by the total number of images) \n"
    ;
}

/* Produce Side information ------------------------------------------------ */
void CorrectAmplitude3DParams::produceSideInfo()
{
    // Read the CTFdat
    ctfdat.read(fnCtfdat);

    // Get dimensions of the volumes
    ctfdat.firstObject();
    FileName fnVol, fnCTF;
    ctfdat.getValue(MDL_IMAGE,fnVol);
    ctfdat.getValue(MDL_CTFMODEL,fnCTF);
    Image<double> V;
    V.read(fnVol);
    unsigned long Ndim;
    V.getDimensions(Xdim,Ydim,Zdim,Ndim);
}

/* Make 3D CTF ------------------------------------------------------------- */
void CorrectAmplitude3DParams::generateCTF1D(const FileName &fnCTF, const double nr_steps,
        MultidimArray<double> &CTF1D)
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
        	ctf.ctf.precomputeValues(res, 0.0);
            CTF1D(step)=ctf.ctf.CTF_at();
            if (isFlipped)
                CTF1D(step)=ABS(CTF1D(step));
        }
        else
            CTF1D(step)=1.;
        res+=stepsize;
    }

}

/* Make Wiener filters ------------------------------------------------------------- */
void CorrectAmplitude3DParams::generateWienerFilters()
{
    MultidimArray<double> CTF1D, sumterm;
    FileName fn_tmp;
    int nrimgs;
    std::ofstream  fh;
    double res;
    double tot_nr_imgs = 0;
    double minsum=99.e99;
    double maxsum=0.;
    double minwien=99.e99;
    double maxwien=0.;
    // Oversample the 1D CTF and Wiener filter vectors OVERSAMPLE times
    // Use 0.55*sqrt(3) to make sure all pixels fit in...
    double nr_steps= CEIL(OVERSAMPLE * 0.55 * sqrt((double)(Zdim*Zdim + Ydim*Ydim + Xdim*Xdim)));

    Vctfs1D.clear();
    Vwien1D.clear();
    int ii = 0;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        // Calculate 1D CTF
        FileName fnVol, fnCTF;
        ctfdat.getValue(MDL_IMAGE,fnVol);
        ctfdat.getValue(MDL_CTFMODEL,fnCTF);
        generateCTF1D(fnCTF,nr_steps,CTF1D);
        Vctfs1D.push_back(CTF1D);

        // Get the number of images contributing to this group
    	int NumberOfImages;
    	ctfdat.getValue(MDL_IMAGE_CLASS_COUNT,NumberOfImages);
        tot_nr_imgs += NumberOfImages;

        // Calculate denominator of the Wiener filter
        if (ii==0)
            sumterm.resize(CTF1D);
        FOR_ALL_ELEMENTS_IN_ARRAY1D(sumterm)
        sumterm(i) += NumberOfImages * CTF1D(i) * CTF1D(i);
        ii++;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D(sumterm)
    {
        // Find min and max values of the sumterm
        if (sumterm(i)>maxsum)
            maxsum=sumterm(i);
        if (sumterm(i)<minsum)
            minsum=sumterm(i);
        // Add (normalized) Wiener filter constant
        sumterm(i) += tot_nr_imgs*wienConst;
    }

    int iimax = ii;
    // Fill the Wiener filter vector
    for (ii=0; ii<iimax; ii++)
    {
    	int NumberOfImages;
    	ctfdat.getValue(MDL_IMAGE_CLASS_COUNT,NumberOfImages);

        FOR_ALL_ELEMENTS_IN_ARRAY1D(CTF1D)
        {
            CTF1D(i) = NumberOfImages * Vctfs1D[ii](i) / sumterm(i);
            if (CTF1D(i)>maxwien)
                maxwien=CTF1D(i);
            if (CTF1D(i)<minwien)
                minwien=CTF1D(i);
        }
        Vwien1D.push_back(CTF1D);

        // Write CTF and Wiener filter curves to disc
        fn_tmp = fnOut + "_wien";
        fn_tmp.compose(fn_tmp, ii+1, "txt");
        fh.open((fn_tmp).c_str(), std::ios::out);
        if (!fh)
            REPORT_ERROR(ERR_IO_NOWRITE, fn_tmp);
        for (int step = 0; step < nr_steps; step++)
        {
            res = (step * sqrt(3.) ) /
                  (OVERSAMPLE * sqrt( (double) (Zdim*Zdim + Ydim*Ydim + Xdim*Xdim) ) );
            if (res<=0.5)
                fh << res << " " << Vwien1D[ii](step) << " " << Vctfs1D[ii](step) << "\n";
        }
        fh.close();
    }

    // Some output to screen
    std::cerr <<" ---------------------------------------------------"<<std::endl;
    std::cerr <<" + Number of defocus groups      = "<<ctfdat.size()<<std::endl;
    std::cerr <<" + Total number of images        = "<<tot_nr_imgs<<std::endl;
    std::cerr <<" + Normalized Wiener constant    = "<<tot_nr_imgs*wienConst<<std::endl;
    std::cerr <<" + Minimum of sum in denominator = "<<minsum<<std::endl;
    std::cerr <<" + Maximum of sum in denominator = "<<maxsum<<std::endl;
    std::cerr <<" + Minimum Wiener filter value   = "<<minwien<<std::endl;
    std::cerr <<" + Maximum Wiener filter value   = "<<maxwien<<std::endl;
    std::cerr <<" ---------------------------------------------------"<<std::endl;
}

/* Make deconvolved volume -------------------------------------------------- */
void CorrectAmplitude3DParams::generateVolumes()
{

    Image<double> V;
    FileName fnVol, fnCTF;
    MultidimArray<std::complex<double> > fft,fft_out;
    Matrix1D<int>    idx(3);
    Matrix1D<double> freq(3);
    int ires;

    int ii = 0;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        ctfdat.getValue(MDL_IMAGE,fnVol);
        ctfdat.getValue(MDL_CTFMODEL,fnCTF);
        V.read(fnVol);
        FourierTransform(V(),fft);
        if (ii == 0)
            fft_out.resize(fft);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(fft)
        {
            XX(idx) = j;
            YY(idx) = i;
            ZZ(idx) = k;
            FFT_idx2digfreq(fft, idx, freq);
            ires= ROUND(OVERSAMPLE*sqrt(XX(freq)*XX(freq)*Xdim*Xdim+YY(freq)*YY(freq)*Ydim*Ydim+ZZ(freq)*ZZ(freq*Zdim*Zdim)));
            fft_out(k,i,j)+=(Vwien1D[ii])(ires)*fft(k,i,j);
        }
        ii++;
    }

    // Inverse Fourier Transform
    InverseFourierTransform(fft_out,V());
    fnVol=fnOut+"_deconvolved.vol";
    V.write(fnVol);

    // Calculate CTF-affected volumes
    ii = 0;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        ctfdat.getValue(MDL_IMAGE,fnVol);
        ctfdat.getValue(MDL_CTFMODEL,fnCTF);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(fft)
        {
            XX(idx) = j;
            YY(idx) = i;
            ZZ(idx) = k;
            FFT_idx2digfreq(fft, idx, freq);
            ires= ROUND(OVERSAMPLE*sqrt(XX(freq)*XX(freq)*Xdim*Xdim+YY(freq)*YY(freq)*Ydim*Ydim+ZZ(freq)*ZZ(freq*Zdim*Zdim)));
            fft(k,i,j)=Vctfs1D[ii](ires)*fft_out(k,i,j);
        }
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
