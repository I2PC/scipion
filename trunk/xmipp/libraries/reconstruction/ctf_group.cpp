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
 
#include "ctf_group.h"

/* Read parameters from command line. -------------------------------------- */
void CtfGroupParams::read(int argc, char **argv)
{
    fn_sel        = getParameter(argc, argv, "-i");
    fn_ctfdat     = getParameter(argc, argv, "-ctfdat");
    fn_root       = getParameter(argc, argv, "-o","ctf_group");
    phase_flipped = checkParameter(argc, argv, "-phase_flipped");
    do_discard_anisotropy = checkParameter(argc, argv, "-discard_anisotropy");
    do_auto       = !checkParameter(argc, argv, "-divide_at");
    if (do_auto)
    {
        max_error     = textToFloat(getParameter(argc, argv, "-error","0.3"));
        resol_error   = textToFloat(getParameter(argc, argv, "-resol"));
    }
    else
    {
        fn_divide = getParameter(argc, argv, "-divide_at");
    }

}

/* Show -------------------------------------------------------------------- */
void CtfGroupParams::show()
{

    std::cerr << "  Input sel file          : "<< fn_sel << std::endl;
    std::cerr << "  Input ctfdat file       : "<< fn_ctfdat << std::endl;
    std::cerr << "  Output rootname         : "<< fn_root << std::endl;
    if (do_discard_anisotropy)
    {
        std::cerr << " -> Exclude anisotropic CTFs from the groups"<<std::endl;
    }
    if (do_auto)
    {
	std::cerr << " -> Using automated mode for making groups"<<std::endl;
        std::cerr << "  Maximum allowed error   : "<<max_error
                  <<" at "<<resol_error<<"Ang resolution"<<std::endl; 
    }
    else
    {
        std::cerr << " -> Group based on defocus values in "<<fn_divide<<std::endl;
    }
    if (phase_flipped)
    {
	std::cerr << " -> Assume that data are PHASE FLIPPED"<<std::endl;
    }
    else
    {
	std::cerr << " -> Assume that data are NOT PHASE FLIPPED"<<std::endl;
    }
}

/* Usage ------------------------------------------------------------------- */
void CtfGroupParams::usage()
{
    std::cerr << "   -sel <selfile>           : Input selfile \n"
              << "   -ctfdat <ctfdat file>    : Input CTFdat file for all data\n"
              << "  [-o <oext=\"wien\">]        : Output root name\n"
	      << "  [-phase_flipped]          : Output filters for phase-flipped data\n"
              << "  [-discard_anisotropy]     : Exclude anisotropic CTFs from groups\n"
              << " MODE 1: AUTOMATED: \n"
              << "  [-error <float=0.3> ]     : Maximum allowed error\n"
              << "  [-resol <float> ]         : Resol. (in Ang) for error calculation (default=Nyquist)\n"
              << " MODE 2: MANUAL: \n"
              << "  [-divide_at <docfile> ]   : 1-column docfile with defocus values \n"
    ;
}

/* Produce Side information ------------------------------------------------ */
void CtfGroupParams::produceSideInfo()
{



    FileName fnt_img, fnt_ctf, fnt;
    ImageXmipp img;
    XmippCTF ctf;
    CTFDat ctfdat;
    SelFile SF;
    Matrix2D<double> Mctf;
    Matrix2D<std::complex<double> >  ctfmask;
    std::vector<FileName> dum;
    int dim, ydim, imgno;
    bool is_unique, found, is_first;
    double avgdef;

    SF.read(fn_sel);
    SF.ImgSize(dim,ydim);
    if ( dim != ydim )
	REPORT_ERROR(1,"ctfgroup ERROR%% Only squared images are allowed!");
    Mctf.resize(dim,dim);

    mics_fnctf.clear();
    SF.go_beginning();
    is_first=true;
    imgno = 0;
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
                is_unique=true;
                for (int i = 0; i< mics_fnctf.size(); i++)
                {
                    if (fnt_ctf == mics_fnctf[i])
                    {
                        is_unique=false;
                        mics_fnimgs[i].push_back(fnt_img);
                        mics_count[i]++;
                        break;
                    }
                }
                if (is_unique)
                {
                    // Read CTF in memory
                    ctf.read(fnt_ctf);
                    ctf.enable_CTF = ctf.enable_CTFnoise = true;
                    ctf.Produce_Side_Info();
                    if (is_first)
                        pixel_size = ctf.Tm;
                    else if (pixel_size != ctf.Tm)
                        REPORT_ERROR(1,"ctf_group: Can not mix CTFs with different sampling rates!");

                    if (!do_discard_anisotropy || isIsotropic(ctf))
                    {
                        avgdef = (ctf.DeltafU+ctf.DeltafV)/2.;
                        ctf.DeltafU = avgdef;
                        ctf.DeltafV = avgdef;
                        ctf.Generate_CTF(dim, dim, ctfmask);
                        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mctf)
                        {
                            if (phase_flipped) dMij(Mctf, i, j) = fabs(dMij(ctfmask, i, j).real());
                            else dMij(Mctf, i, j) = dMij(ctfmask, i, j).real();
                        }
                        // Fill vectors
                        mics_fnctf.push_back(fnt_ctf);
                        mics_count.push_back(1);
                        mics_fnimgs.push_back(dum);
                        mics_fnimgs[mics_fnimgs.size()-1].push_back(fnt_img);
                        mics_ctf2d.push_back(Mctf);
                        mics_defocus.push_back(avgdef);
                    }
                    else
                    {
                        std::cerr<<"Discard CTF "<<fnt_ctf<<" because of too large anisotropy"<<std::endl;
                    }
                }
		break;
	    }
	}
	if (!found)
	    REPORT_ERROR(1, "ctf_group ERROR%% Did not find image "+fnt+" in the CTFdat file");
        imgno++;
    }

    // Set resolution limits in dig freq:
    resol_error = pixel_size / resol_error;
    // and in pixels:
    iresol_error = ROUND(resol_error * dim);

}

// Check whether a CTF is anisotropic
bool CtfGroupParams::isIsotropic(XmippCTF &ctf)
{
    double xp, yp, xpp, ypp;
    double cosp, sinp, ctfp, diff;
    Matrix1D<double> freq(2);

    cosp = COSD(ctf.azimuthal_angle);
    sinp = SIND(ctf.azimuthal_angle);

    for (double digres = 0; digres < resol_error; digres+= 0.001)
    {
        XX(freq) = cosp * digres;
        YY(freq) = sinp * digres;
        digfreq2contfreq(freq, freq, pixel_size);
        ctfp = ctf.CTF_at(XX(freq), YY(freq));
        diff = ABS(ctfp - ctf.CTF_at(YY(freq), XX(freq)));
        if (diff > max_error)
        {
            return false;
        }
    }
    return true;

}

// Do the actual work
void CtfGroupParams::autoRun()
{

    double diff, mindiff;
    int iopt_group, nr_groups;
    std::vector<int> dum;

    nr_groups = 0;
    for (int imic=0; imic < mics_ctf2d.size(); imic++)
    {
        mindiff = 99999.;
        iopt_group = -1;
        for (int igroup=0; igroup < nr_groups; igroup++)
        {
            // loop over all mics in this group
            for (int igmic=0; igmic < pointer_group2mic[igroup].size(); igmic++)
            {
                bool is_ok = true;
                for (int iresol=0; iresol<=iresol_error; iresol++)
                {
                    diff = ABS( dMij(mics_ctf2d[imic],iresol,0) - 
                                dMij(mics_ctf2d[pointer_group2mic[igroup][igmic]],iresol,0) );
                    if (diff > max_error)
                    {
                        is_ok = false;
                        break;
                    }
                    else if (diff < mindiff)
                    {
                        mindiff = diff;
                        iopt_group = igroup;
                    }
                }
                if (!is_ok)
                {
                    break;
                }
            }
        }
        if (iopt_group < 0)
        {
            //add to existing group
            pointer_group2mic[iopt_group].push_back(imic);
            pointer_mic2group.push_back(iopt_group);
        }
        else
        {
            // add new group
            pointer_group2mic.push_back(dum);
            pointer_group2mic[nr_groups].push_back(imic);
            pointer_mic2group.push_back(nr_groups);
            nr_groups++;
        }
    }

    // I/O: For each group: 
    //
    // 1. write selfile with images
    // 2. write average Mctf
    // 3. write file with defocus values
    // 4. write 1D profile files
    // 5. write CTF-param file
    FileName fnt;
    SelFile SFo;
    ImageXmipp img;
    Matrix2D<double> Mavg;
    double sumw, avgdef;
    int imic;
    std::ofstream fh, fh2;

    for (int igroup=0; igroup < nr_groups; igroup++)
    {
        SFo.clear();
        Mavg.initZeros(dim,dim);
        sumw = 0.;
        avgdef = 0.;
        fh.open((fnt+".defocus").c_str(), std::ios::out);
        for (int igmic=0; igmic < pointer_group2mic[igroup].size(); igmic++)
        {
            imic = pointer_group2mic[igroup][igmic];
            sumw += (double) mics_count[imic];
            // calculate (weighted) average Mctf
            Mavg += mics_count[imic] * mics_ctf2d[imic];
            // calculate (weighted) average defocus
            avgdef += mics_count[imic] * mics_defocus[imic];
            // Fill SelFile
            for (int iimg=0; iimg < mics_fnimgs[imic].size(); iimg++)
            {
                SFo.insert(mics_fnimgs[imic][iimg]);
            }
            // Fill file with defocus values 
            fh << floatToString(mics_defocus[imic]) <<"\n";
        }
        Mavg /= sumw;
        avgdef /= sumw;
        fnt.compose(fnt+"_group",igroup+1,"");
        // 1. write selfile
        SFo.write(fnt+".sel");
        // 2. write average Mctf
        img() = Mavg;
        img.weight() = sumw;
        img.write(fnt+".fft");
        // 3. Write files with defocus values 
        fh.close();
        // 4. Write file with 1D profiles
        fh2.open((fnt+".profiles").c_str(), std::ios::out);
        for (int i=0; i < dim/2; i++)
        {
            fh2 << floatToString((double)i / (pixel_size*dim));
            fh2.width(10);
            fh2 << floatToString(dMij(Mavg,i,0));
            fh2.width(10);
            for (int igmic=0; igmic < pointer_group2mic[igroup].size(); igmic++)
            {
                imic = pointer_group2mic[igroup][igmic];
                fh2 << floatToString(dMij(mics_ctf2d[imic],i,0));
                fh2.width(10);
            }
        }
        fh2.close();
        // 5. Write CTF parameter file
        // TODO
    }

}




