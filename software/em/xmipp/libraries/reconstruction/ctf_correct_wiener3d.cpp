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

#include "ctf_correct_wiener3d.h"

#include <data/xmipp_fft.h>

#define OVERSAMPLE 8

/* Read parameters from command line. -------------------------------------- */
void ProgCtfCorrectAmplitude3D::readParams()
{
    fnIn  = getParam("-i");
    fnRoot = getParam("--oroot");
    minFreq  = getDoubleParam("--minFreq");
    wienerConstant = getDoubleParam("--wienerConstant");
    isFlipped = checkParam("--phase_flipped");
}

/* Show -------------------------------------------------------------------- */
void ProgCtfCorrectAmplitude3D::show()
{
    if (verbose==0)
        return;
    std::cout
    << "Input file:      " << fnIn << std::endl
    << "Output rootname: " << fnRoot << std::endl
    << "Wiener constant: " << wienerConstant << std::endl
    << "Phase flipped:   " << isFlipped << std::endl
    ;
    if (minFreq>0)
        std::cout << "Apply Wiener filter only beyond:   " << minFreq << std::endl;
}

/* Usage ------------------------------------------------------------------- */
void ProgCtfCorrectAmplitude3D::defineParams()
{
    addUsageLine("Wiener filtering of volumes");
    addUsageLine("+The program combines a set of volumes, each one with its one CTF and produces a deconvolved Wiener volume");

    addParamsLine("  -i <metadataFile>           : Metadata with the volumes, ctfs, number of images in that group");
    addParamsLine("                              :+The metadata labels are _image, _CTFModel, _class_count");
    addParamsLine("  --oroot <file>              : Output rootname ");
    addParamsLine("                              :+oroot+_deconvolved.vol contains the combination of all volumes");
    addParamsLine("                              :+oroot+_ctffiltered_group01.vol contains each volume obtained after filtering the deconvolved one");
    addParamsLine("                              :+With verbose==2: oroot+_wien01.txt contains the Wiener filters in Fourier");
    addParamsLine("  [--minFreq <Ang=-1>]        : Apply Wiener filter only beyond this resolution (in Angstrom)");
    addParamsLine("  [--phase_flipped]           : Use this if the maps were reconstructed from phase corrected images ");
    addParamsLine("  [--wienerConstant <K=0.05>] : Wiener constant (to be multiplied by the total number of images) ");
    addExampleLine("xmipp_ctf_correct_wiener3d -i ctf_correct3d.xmd --oroot volumeCorrected");
    addExampleLine("In the following link you can find an example of input file:",false);
    addExampleLine(" ",false);
    addExampleLine("http://sourceforge.net/p/testxmipp/code/ci/master/tree/input/ctf_correct3d.xmd?format=raw",false);
}

/* Produce Side information ------------------------------------------------ */
void ProgCtfCorrectAmplitude3D::produceSideInfo()
{
    // Read the CTFdat
    ctfdat.read(fnIn);

    // Get dimensions of the volumes
    size_t id = ctfdat.firstObject();
    FileName fnVol, fnCTF;
    ctfdat.getValue(MDL_IMAGE,fnVol, id);
    ctfdat.getValue(MDL_CTF_MODEL,fnCTF, id);
    Image<double> V;
    V.read(fnVol, HEADER);
    size_t Ndim;
    V.getDimensions(Xdim,Ydim,Zdim,Ndim);
}

/* Make 3D CTF ------------------------------------------------------------- */
void ProgCtfCorrectAmplitude3D::generateCTF1D(const FileName &fnCTF, size_t nr_steps,
        MultidimArray<double> &CTF1D)
{
    // Read the CTF
    ctf.FilterBand = CTF;
    ctf.ctf.enable_CTFnoise = false;
    ctf.ctf.read(fnCTF);
    ctf.ctf.produceSideInfo();

    double maxres = ( 0.5 * sqrt(3.) ) / ctf.ctf.Tm;
    double stepsize = maxres / (double)nr_steps;
    CTF1D.resizeNoCopy(nr_steps);
    double freq = 0.;

    freq=0.;
    for (size_t step=0; step < nr_steps; step++)
    {
        if ( (minFreq < 0) || (1./freq < minFreq) )
        {
            ctf.ctf.precomputeValues(freq, 0.0);
            A1D_ELEM(CTF1D,step)=ctf.ctf.getValueAt();
            if (isFlipped)
                A1D_ELEM(CTF1D,step)=fabs(CTF1D(step));
        }
        else
            A1D_ELEM(CTF1D,step)=1.;
        freq+=stepsize;
    }
}

/* Make Wiener filters ------------------------------------------------------------- */
void ProgCtfCorrectAmplitude3D::generateWienerFilters()
{
    MultidimArray<double> CTF1D, sumterm;
    std::ofstream  fh;
    double res;
    double tot_nr_imgs = 0;
    // Oversample the 1D CTF and Wiener filter vectors OVERSAMPLE times
    // Use 0.55*sqrt(3) to make sure all pixels fit in...
    size_t nr_steps= (size_t)ceil(OVERSAMPLE * 0.55 * sqrt((double)(Zdim*Zdim + Ydim*Ydim + Xdim*Xdim)));

    Vctfs1D.clear();
    Vwien1D.clear();
    int ii = 0;
    FileName fnVol, fnCTF;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        // Calculate 1D CTF
        ctfdat.getValue(MDL_IMAGE,fnVol,__iter.objId);
        ctfdat.getValue(MDL_CTF_MODEL,fnCTF,__iter.objId);
        generateCTF1D(fnCTF,nr_steps,CTF1D);
        Vctfs1D.push_back(CTF1D);

        // Get the number of images contributing to this group
        size_t NumberOfImages;
        ctfdat.getValue(MDL_CLASS_COUNT,NumberOfImages,__iter.objId);
        tot_nr_imgs += NumberOfImages;

        // Calculate denominator of the Wiener filter
        if (ii==0)
            sumterm.resizeNoCopy(CTF1D);
        FOR_ALL_ELEMENTS_IN_ARRAY1D(sumterm)
        {
            double CTF1Di=A1D_ELEM(CTF1D,i);
            A1D_ELEM(sumterm,i) += NumberOfImages * CTF1Di * CTF1Di;
        }
        ii++;
    }

    // Add (normalized) Wiener filter constant
    sumterm += tot_nr_imgs*wienerConstant;

    double maxwien=0, minwien=1e38;
    ii=0;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        size_t NumberOfImages;
        ctfdat.getValue(MDL_CLASS_COUNT,NumberOfImages,__iter.objId);
        CTF1D=Vctfs1D[ii];
        CTF1D*=NumberOfImages;
        CTF1D/=sumterm;
        Vwien1D.push_back(CTF1D);

        // Write CTF and Wiener filter curves to disc
        if (verbose>2)
        {
            FileName fn_tmp;
            fn_tmp = fnRoot + "_wien";
            fn_tmp.compose(fn_tmp, ii+1, "txt");
            fh.open((fn_tmp).c_str(), std::ios::out);
            if (!fh)
                REPORT_ERROR(ERR_IO_NOWRITE, fn_tmp);
            for (size_t step = 0; step < nr_steps; step++)
            {
                res = (step * sqrt(3.) ) /
                      (OVERSAMPLE * sqrt( (double) (Zdim*Zdim + Ydim*Ydim + Xdim*Xdim) ) );
                if (res<=0.5)
                    fh << res << " " << Vwien1D[ii](step) << " " << Vctfs1D[ii](step) << "\n";
            }
            fh.close();
        }
        ii++;
    }

    // Some output to screen
    std::cerr <<" ---------------------------------------------------"<<std::endl;
    std::cerr <<" + Number of defocus groups      = "<<ctfdat.size()<<std::endl;
    std::cerr <<" + Total number of images        = "<<tot_nr_imgs<<std::endl;
    std::cerr <<" + Normalized Wiener constant    = "<<tot_nr_imgs*wienerConstant<<std::endl;
    std::cerr <<" + Minimum Wiener filter value   = "<<minwien<<std::endl;
    std::cerr <<" + Maximum Wiener filter value   = "<<maxwien<<std::endl;
    std::cerr <<" ---------------------------------------------------"<<std::endl;
}

/* Make deconvolved volume -------------------------------------------------- */
void ProgCtfCorrectAmplitude3D::generateVolumes()
{
    Image<double> V;
    FileName fnVol, fnCTF;
    MultidimArray<std::complex<double> > fft,fft_out;
    Matrix1D<int>    idx(3);
    Matrix1D<double> freq(3);

    int ii = 0;
    double Xdim2=Xdim*Xdim;
    double Ydim2=Ydim*Ydim;
    double Zdim2=Zdim*Zdim;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        ctfdat.getValue(MDL_IMAGE,fnVol,__iter.objId);
        ctfdat.getValue(MDL_CTF_MODEL,fnCTF,__iter.objId);
        V.read(fnVol);
        FourierTransform(V(),fft);
        if (ii == 0)
            fft_out.initZeros(fft);
        MultidimArray<double>& Vwien1D_ii=Vwien1D[ii];
        FOR_ALL_ELEMENTS_IN_ARRAY3D(fft)
        {
            XX(idx) = j;
            YY(idx) = i;
            ZZ(idx) = k;
            FFT_idx2digfreq(fft, idx, freq);
            int ires= (int)round(OVERSAMPLE*sqrt(XX(freq)*XX(freq)*Xdim2+
                                                 YY(freq)*YY(freq)*Ydim2+
                                                 ZZ(freq)*ZZ(freq)*Zdim2));
            A3D_ELEM(fft_out,k,i,j)+=A1D_ELEM(Vwien1D_ii,ires)*A3D_ELEM(fft,k,i,j);
        }
        ii++;
    }

    // Inverse Fourier Transform
    InverseFourierTransform(fft_out,V());
    fnVol=fnRoot+"_deconvolved.vol";
    V.write(fnVol);

    // Calculate CTF-affected volumes
    ii = 0;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        ctfdat.getValue(MDL_IMAGE,fnVol,__iter.objId);
        ctfdat.getValue(MDL_CTF_MODEL,fnCTF,__iter.objId);
        MultidimArray<double>& Vwien1D_ii=Vwien1D[ii];
        FOR_ALL_ELEMENTS_IN_ARRAY3D(fft)
        {
            XX(idx) = j;
            YY(idx) = i;
            ZZ(idx) = k;
            FFT_idx2digfreq(fft, idx, freq);
            int ires= (int)round(OVERSAMPLE*sqrt(XX(freq)*XX(freq)*Xdim2+
                                                 YY(freq)*YY(freq)*Ydim2+
                                                 ZZ(freq)*ZZ(freq)*Zdim2));
            A3D_ELEM(fft,k,i,j)+=A1D_ELEM(Vwien1D_ii,ires)*A3D_ELEM(fft_out,k,i,j);
        }
        ii++;

        fnVol=fnRoot+"_ctffiltered_group";
        fnVol.compose(fnVol,ii,"vol");
        InverseFourierTransform(fft,V());
        V.write(fnVol);
    }
}

/* Correct a set of images ------------------------------------------------- */
void ProgCtfCorrectAmplitude3D::run()
{
    produceSideInfo();
    generateWienerFilters();
    generateVolumes();
}
