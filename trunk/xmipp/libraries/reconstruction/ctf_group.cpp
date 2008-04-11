
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
    fn_root       = getParameter(argc, argv, "-o","ctf");
    phase_flipped = checkParameter(argc, argv, "-phase_flipped");
    do_discard_anisotropy = checkParameter(argc, argv, "-discard_anisotropy");
    do_auto       = !checkParameter(argc, argv, "-split");
    if (do_auto)
    {
        max_error     = textToFloat(getParameter(argc, argv, "-error","0.5"));
        resol_error   = textToFloat(getParameter(argc, argv, "-resol"));
    }
    else
    {
        fn_split = getParameter(argc, argv, "-split");
    }
    do_wiener = checkParameter(argc, argv, "-wiener");
    wiener_constant = textToFloat(getParameter(argc, argv, "-wc","-1"));

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
        std::cerr << " -> With a maximum allowed error of "<<max_error
                  <<" at "<<resol_error<<" Ang resolution"<<std::endl; 
    }
    else
    {
        std::cerr << " -> Group based on defocus values in "<<fn_split<<std::endl;
    }
    if (phase_flipped)
    {
	std::cerr << " -> Assume that data are PHASE FLIPPED"<<std::endl;
    }
    else
    {
	std::cerr << " -> Assume that data are NOT PHASE FLIPPED"<<std::endl;
    }
    if (do_wiener)
    {
	std::cerr << " -> Also calculate Wiener filters, with constant= "<<wiener_constant<<std::endl;
    }
    std::cerr << "----------------------------------------------------------"<<std::endl;
}

/* Usage ------------------------------------------------------------------- */
void CtfGroupParams::usage()
{
    std::cerr << "   -i <selfile>             : Input selfile \n"
              << "   -ctfdat <ctfdat file>    : Input CTFdat file for all data\n"
              << "  [-o <oext=\"ctf\">]        : Output root name\n"
	      << "  [-phase_flipped]          : Output filters for phase-flipped data\n"
              << "  [-discard_anisotropy]     : Exclude anisotropic CTFs from groups\n"
              << "  [-wiener]                 : Also calculate Wiener filters\n"
              << "  [-wc <float=-1>]          : Wiener-filter constant (if < 0: use FREALIGN default)\n"
              << " MODE 1: AUTOMATED: \n"
              << "  [-error <float=0.5> ]     : Maximum allowed error\n"
              << "  [-resol <float> ]         : Resol. (in Ang) for error calculation (default=Nyquist)\n"
              << " MODE 2: MANUAL: \n"
              << "  [-split <docfile> ]       : 1-column docfile with defocus values where to split the data \n"
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
    int ydim, imgno;
    bool is_unique, found, is_first;
    double avgdef;
    std::vector<FileName> mics_fnctf;

    SF.read(fn_sel);
    SF.ImgSize(dim,ydim);
    if ( dim != ydim )
	REPORT_ERROR(1,"ctfgroup ERROR%% Only squared images are allowed!");
    Mctf.resize(dim,dim);
    ctfdat.read(fn_ctfdat);

    if (do_wiener)
    {
        Mwien.resize(dim,dim);
        Mwien.initZeros();
    }

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
            ctfdat.nextLine();
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
                        //std::cerr<<" add "<<fnt_img<<" to mic "<<i<<"count= "<<mics_count[i]<<std::endl;
                        break;
                    }
                }
                if (is_unique)
                {
                    // Read CTF in memory
                    ctf.read(fnt_ctf);
                    ctf.enable_CTF = true; 
                    ctf.enable_CTFnoise = false;
                    ctf.Produce_Side_Info();
                    if (is_first)
                    {
                        pixel_size = ctf.Tm;
                        // Set resolution limits in dig freq:
                        resol_error = pixel_size / resol_error;
                        resol_error = XMIPP_MIN(0.5, resol_error);
                        // and in pixels:
                        iresol_error = ROUND(resol_error * dim);
                        std::cerr<<" Resolution for error limit = "<<resol_error<< " (dig. freq.) = "<<iresol_error<<" (pixels)"<<std::endl;
                        is_first=false;
                    }
                    else if (pixel_size != ctf.Tm)
                    {
                        REPORT_ERROR(1,"ctf_group: Can not mix CTFs with different sampling rates!");
                    }
                    if (!do_discard_anisotropy || isIsotropic(ctf))
                    {
                        avgdef = (ctf.DeltafU + ctf.DeltafV)/2.;
                        ctf.DeltafU = avgdef;
                        ctf.DeltafV = avgdef;
                        ctf.Generate_CTF(dim, dim, ctfmask);
                        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(ctfmask)
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
                        //std::cerr<<" uniq "<<fnt_img<<" in mic "<<mics_fnimgs.size()-1<<"def= "<<avgdef<<std::endl;
                    }
                    else
                    {
                        std::cerr<<" Discard CTF "<<fnt_ctf<<" because of too large anisotropy"<<std::endl;
                    }
                }
		break;
	    }
	}
	if (!found)
	    REPORT_ERROR(1, "ctf_group ERROR%% Did not find image "+fnt+" in the CTFdat file");
        imgno++;
    }

    // Precalculate denominator term of the Wiener filter
    if (do_wiener)
    {
        double sumimg = 0.;
        for (int imic=0; imic < mics_count.size(); imic++)
        {
            sumimg += (double)mics_count[imic];
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mwien)
            {
                dMij(Mwien,i,j) += mics_count[imic] * dMij(mics_ctf2d[imic],i,j) * dMij(mics_ctf2d[imic],i,j); 
            }
        }
        // Use Grigorieff's default for Wiener filter constant: 10% of average over all Mwien terms
        // Grigorieff JSB 157(1) (2006), pp 117-125
        if (wiener_constant < 0.) 
        {
            wiener_constant = 0.1 * Mwien.computeAvg();
        }
        // Also divide by sumimg (Wiener filter is for summing images, not averaging!)
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mwien)
        {
            dMij(Mwien,i,j) += wiener_constant;
            dMij(Mwien,i,j) /= sumimg;
        }
        
    }
 
    // Now that we have all information, sort micrographs on defocus value
    // And update all corresponding pointers and data vectors.
    std::vector<Matrix2D<double> > newmics_ctf2d=mics_ctf2d;
    std::vector<int> newmics_count=mics_count;
    std::vector< std::vector <FileName> > newmics_fnimgs=mics_fnimgs;
    std::vector<double> sorted_defocus=mics_defocus;

    std::sort(sorted_defocus.begin(), sorted_defocus.end());
    std::reverse(sorted_defocus.begin(), sorted_defocus.end());
    for (int isort = 0; isort < sorted_defocus.size(); isort++)
        for (int imic=0; imic < mics_defocus.size(); imic++)
            if (ABS(sorted_defocus[isort] -  mics_defocus[imic]) < XMIPP_EQUAL_ACCURACY)
            {
                newmics_ctf2d[isort]  = mics_ctf2d[imic];
                newmics_count[isort]  = mics_count[imic];
                newmics_fnimgs[isort] = mics_fnimgs[imic];
                break;
            }
    mics_ctf2d   = newmics_ctf2d;
    mics_count   = newmics_count;
    mics_fnimgs  = newmics_fnimgs;
    mics_defocus = sorted_defocus;
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
            std::cerr<<" Anisotropy!"<<digres<<" "<<max_error<<" "<<diff<<" "<<ctfp<<" "<<ctf.CTF_at(YY(freq), XX(freq))<<std::endl;
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

    // Make the actual groups
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

                for (int iresol=0; iresol<=iresol_error; iresol++)
                {
                    diff = ABS( dMij(mics_ctf2d[imic],iresol,0) - 
                                dMij(mics_ctf2d[pointer_group2mic[igroup][igmic]],iresol,0) );
                    if (diff > max_error)
                    {
                        break;
                    }
                }
                if (diff <= max_error && diff < mindiff)
                {
                    mindiff = diff;
                    iopt_group = igroup;
                }
            }
        }
        if (iopt_group < 0)
        {
            // add new group
            pointer_group2mic.push_back(dum);
            pointer_group2mic[nr_groups].push_back(imic);
            nr_groups++;
        }
        else
        {
            //add to existing group
            pointer_group2mic[iopt_group].push_back(imic);
        }
    }
    std::cerr<<" Number of CTF groups= "<<nr_groups<<std::endl;
}
   

void CtfGroupParams::manualRun()
{
    DocFile DF;
    double defocus, split, oldsplit=99.e99;
    int count, nr_groups = 0;
    Matrix1D<double> dataline(1);
    std::vector<int> dum;

    pointer_group2mic.clear();
    DF.read(fn_split);    
    // Append -99.e-99 to the end of the docfile
    DF.locate(DF.get_last_key());
    dataline(0)=-99.e99;
    DF.append_data_line(dataline);

    // Loop over all splits
    DF.go_first_data_line();
    while (!DF.eof())
    {
        split=DF(0);
        count = 0;
        for (int imic=0; imic < mics_ctf2d.size(); imic++)
        {
            defocus = mics_defocus[imic];
            if (defocus < oldsplit && defocus >= split)
            {
                if (count == 0)
                {
                    // add new group
                    pointer_group2mic.push_back(dum);
                }
                pointer_group2mic[nr_groups].push_back(imic);
                count++;
            }
        }
        if (count == 0)
        {
            std::cerr<<" Warning: group with defocus values between "<<oldsplit<<" and "<<split<<" is empty!"<<std::endl;
        }
        else
        {
            nr_groups++;
        }
        oldsplit = split;
        DF.next_data_line();
    }

}

void CtfGroupParams::writeOutputToDisc()
{

    FileName fnt;
    SelFile SFo;
    DocFile DFo;
    Matrix1D<double> dataline(1);
    ImageXmipp img;
    Matrix2D<double> Mavg;
    double sumw, avgdef, mindef, maxdef, split, oldmin=99.e99;
    int imic;
    std::ofstream fh, fh2, fh3;

    
    DFo.append_comment("Defocus values to split into "+integerToString(pointer_group2mic.size())+" ctf groups");
    fh.open((fn_root+ "_groups.imgno").c_str(), std::ios::out);
    fh  << "# Number of images in each group \n";
    fh3.open((fn_root+ "_groups.defocus").c_str(), std::ios::out);
    fh3 << "# Defocus values for each group (avg, max & min) \n";

    for (int igroup=0; igroup < pointer_group2mic.size(); igroup++)
    {
        SFo.clear();
        Mavg.initZeros(dim,dim);
        sumw = 0.;
        avgdef = 0.;
        mindef = 99.e99;
        maxdef = -99.e99;
        for (int igmic=0; igmic < pointer_group2mic[igroup].size(); igmic++)
        {
            imic = pointer_group2mic[igroup][igmic];
            sumw += (double) mics_count[imic];
            // calculate (weighted) average Mctf
            Mavg += mics_count[imic] * mics_ctf2d[imic];
            // Calculate avg, min and max defocus values in this group
            avgdef += mics_count[imic] * mics_defocus[imic];
            mindef = XMIPP_MIN(mics_defocus[imic],mindef);
            maxdef = XMIPP_MAX(mics_defocus[imic],maxdef);
            // Fill SelFile
            for (int iimg=0; iimg < mics_fnimgs[imic].size(); iimg++)
            {
                SFo.insert(mics_fnimgs[imic][iimg]);
            }
        }
        Mavg /= sumw;
        avgdef /= sumw;
        fnt.compose(fn_root+"_group",igroup+1,"");
        std::cerr<<" Group "<<fnt <<" contains "<< pointer_group2mic[igroup].size()<<" ctfs and "<<sumw<<" images and has average defocus "<<avgdef<<std::endl;
        
        // 1. write selfile
        SFo.write(fnt+".sel");
        // 2. write average Mctf
        img() = Mavg;
        img.weight() = sumw;
        img.write(fnt+".ctf");
        // 3. Output to file with number of images per group
        fh << integerToString(igroup+1);
        fh.width(10);
        fh << floatToString(sumw)<<std::endl;
        // 4. Output to file with avgdef, mindef and maxdef per group
        fh3 << integerToString(igroup+1);
        fh3.width(10);
        fh3 << floatToString(avgdef);
        fh3.width(10);
        fh3 << floatToString(maxdef);
        fh3.width(10);
        fh3 << floatToString(mindef)<<std::endl;
        // 5. Output to docfile for manual grouping
        if (oldmin < 9.e99)
        {
            split = (oldmin + maxdef) / 2.;
            dataline(0) = split;
            DFo.append_data_line(dataline);
        }
        oldmin = mindef;
        // 6. Write file with 1D profiles
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
            fh2 << std::endl;
        }
        fh2.close();
        // 7. Write Wiener filter
        if (do_wiener)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Mavg)
            {
                dMij(Mavg,i,j) /= dMij(Mwien,i,j);
            }
            img() = Mavg;
            img.weight() = sumw;
            img.write(fnt+".wien");
        }
    }
    
    // 3. Write file with number of images per group
    fh.close();
    // 4. Write file with avgdef, mindef and maxdef per group
    fh3.close();
    // 5. Write docfile with defocus values to split manually
    DFo.write(fn_root+"_groups_split.doc");

}

void CtfGroupParams::run()
{

    if (do_auto)
    {
        autoRun();
    }
    else
    {
        manualRun();
    }

    writeOutputToDisc();

    std::cerr << " Done!" <<std::endl;
    
}

