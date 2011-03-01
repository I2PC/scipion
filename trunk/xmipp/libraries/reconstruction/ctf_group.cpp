
/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
 * Rewritten by Roberto Marabini roberto@cnb.uam.es
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

#include "ctf_group.h"
#include <data/fft.h>
#include <data/metadata_extension.h>
#include <data/metadata_extension.h>

/* Read parameters from command line. -------------------------------------- */
void ProgCtfGroup::readParams()
{
    //    fn_sel        = getParam("-i");
    fn_ctfdat     = getParam("--ctfdat");
    fn_root       = getParam("-o");

    format = fn_root.getFileFormat();
    fn_root = fn_root.removeFileFormat();

    phase_flipped = checkParam("--phase_flipped");
    do_discard_anisotropy = checkParam("--discard_anisotropy");
    do_auto       = !checkParam("--split");
    pad           = XMIPP_MAX(1., getDoubleParam("--pad"));
    if (do_auto)
    {
        max_error     = getDoubleParam("--error");
        resol_error   = getDoubleParam("--resol");
    }
    else
    {
        fn_split = getParam("--split");
    }
    do_wiener = checkParam("--wiener");
    memory = getDoubleParam("--memory");
    do1Dctf   = checkParam("--do1Dctf");
    wiener_constant = getDoubleParam("--wc");
}

/* Show -------------------------------------------------------------------- */
void ProgCtfGroup::show()
{
    //    std::cerr << "  Input sel file          : "<< fn_sel << std::endl;
    std::cerr << "  Input ctfdat file       : "<< fn_ctfdat << std::endl;
    std::cerr << "  Output rootname         : "<< fn_root << std::endl;
    if (pad > 1.)
    {
        std::cerr << "  Padding factor          : "<< pad << std::endl;
    }
    if (do_discard_anisotropy)
    {
        std::cerr << " -> Exclude anisotropic CTFs from the groups"<<std::endl;
    }
    if (do_auto)
    {
        std::cerr << " -> Using automated mode for making groups"<<std::endl;
        std::cerr << " -> With a maximum allowed error of "<<max_error
        <<" at "<<resol_error<<" dig. freq."<<std::endl;
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
    std::cerr << "  Available memory (Gb)        : "<< memory << std::endl;
    if (do1Dctf)
    {
        std::cerr << " -> compute CTF groups using 1D CTFs= "<<std::endl;
    }
    else
    {
        std::cerr << " -> compute CTF groups using 2D CTFs= "<<std::endl;
    }
    std::cerr << "----------------------------------------------------------"<<std::endl;
}

/* Usage ------------------------------------------------------------------- */
void ProgCtfGroup::defineParams()
{
    addUsageLine("Generate CTF (or defocus) groups from a single CTFdat file");
    addUsageLine("Example of use: Sample using automated mode (resolution = 15 Ang.)");
    addUsageLine("   xmipp_ctf_group -i input.sel --ctfdat input.ctfdat --pad 1.5 --wiener --phase_flipped --resol 15");
    addUsageLine("Example of use: Sample using manual mode (after manual editing of ctf_group_split.doc)");
    addUsageLine("   xmipp_ctf_group -i input.sel --ctfdat input.ctfdat --pad 1.5 --wiener --phase_flipped --split ctf_group_split.doc");

    //    addParamsLine("   -i <sel_file>              : Input selfile");
    addParamsLine("   --ctfdat <ctfdat_file>     : Input CTFdat file for all data");
    addParamsLine("   [-o <oext=\"ctf\">]        : Output root name, you may force format ctf:mrc");
    addParamsLine("   [--pad <float=1>]          : Padding factor ");
    addParamsLine("   [--phase_flipped]          : Output filters for phase-flipped data");
    addParamsLine("   [--discard_anisotropy]     : Exclude anisotropic CTFs from groups");
    addParamsLine("   [--wiener]                 : Also calculate Wiener filters");
    addParamsLine("   [--memory <double=1.>]     : Available memory in Gb");
    addParamsLine("   [--do1Dctf]                : Compute Groups using 1D CTF, select this option is you have many \
                  non astismatic CTFs");
    addParamsLine("   [--wc <float=-1>]          : Wiener-filter constant (if < 0: use FREALIGN default)");
    addParamsLine(" == MODE 1: AUTOMATED: == ");
    addParamsLine("   [--error <float=0.5> ]     : Maximum allowed error");
    addParamsLine("   [--resol <float=-1> ]      : Resol. (in Ang) for error calculation. Default (-1) = Nyquist");
    addParamsLine(" == MODE 2: MANUAL: == ");
    addParamsLine("   [--split <docfile> ]       : 1-column docfile with defocus values where to split the data ");

}

/* Produce Side information ------------------------------------------------ */
void ProgCtfGroup::produceSideInfo()
{
    ;
    //#ifdef GGGG

    FileName  fnt_ctf;
    CTFDescription ctf;
    //MetaData ctfdat,
    MetaData SF;
    MultidimArray<double> Mctf;
    MultidimArray<std::complex<double> >  ctfmask;
    int ydim;
    int zdim;
    size_t ndim;
    double avgdef;

    SF.read(fn_ctfdat);
    ImgSize(SF,dim,ydim,zdim,ndim);

    if ( dim != ydim )
        REPORT_ERROR(ERR_MULTIDIM_SIZE,"Only squared images are allowed!");

    paddim=xpaddim = ROUND(pad*dim);
    if(do1Dctf)
    {
        ypaddim=1;
        ctfxpaddim =  sqrt(2.) *  xpaddim + 1;
    }
    else
    {
        //ypaddim = xpaddim;
        //This is ready for the day in which we use anisotropic ctf
        ypaddim=1;
        //ctfxpaddim = xpaddim;
        ctfxpaddim =  sqrt(2.) *  xpaddim + 1;
    }
    Mctf.resize(ypaddim,ctfxpaddim);

    if (do_wiener)
    {
        Mwien.resize(paddim,paddim);
        Mwien.initZeros();
    }

    MetaData ctfMD;
    //number of different CTFs
    ctfMD.aggregate(SF, AGGR_COUNT,MDL_CTFMODEL,MDL_CTFMODEL,MDL_COUNT);
    ctfMD.fillExpand(MDL_CTFMODEL);
    int nCTFs = ctfMD.size();
    //how much memory do I need to store them
    double _sizeGb = (double) ypaddim * xpaddim * sizeof(double) * nCTFs /1073741824.;
    if (_sizeGb > memory)
        mmapOn=true;
    else
        mmapOn=false;
    mics_ctf2d.setMmap(mmapOn);
    mics_ctf2d.resize(nCTFs,1,ypaddim,ctfxpaddim);

    int c = XMIPP_MAX(1, nCTFs / 60);
    init_progress_bar(nCTFs);

    ctf.readFromMetadataRow(ctfMD,ctfMD.firstObject());

    //do not read directly Tm from metadata because it may be not there
    pixel_size = ctf.Tm;

    if (do_auto)
    {
        if (resol_error < 0)
        {
            // Set to Nyquist
            resol_error = 2. * pixel_size;
        }
        // Set resolution limits in dig freq:
        resol_error = pixel_size / resol_error;
        resol_error = XMIPP_MIN(0.5, resol_error);
        // and in pixels:

        iresol_error = ROUND(resol_error * paddim);
    }

    //fill multiarray with ctfs
    int count;
    int counter=0;
    if (verbose!=0)
        std::cout << "\nFill multiarray with ctfs" <<std::endl;
    FOR_ALL_OBJECTS_IN_METADATA(ctfMD)
    {
        ctf.readFromMetadataRow(ctfMD, __iter.objId);
        ctf.enable_CTF = true;
        ctf.enable_CTFnoise = false;
        ctf.Produce_Side_Info();
        if (pixel_size != ctf.Tm)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "Cannot mix CTFs with different sampling rates!");
        ctf.Tm /= sqrt(2.);
        if (!do_discard_anisotropy || isIsotropic(ctf))
        {
            avgdef = (ctf.DeltafU + ctf.DeltafV)/2.;
            ctf.DeltafU = avgdef;
            ctf.DeltafV = avgdef;
            ctf.Generate_CTF(ypaddim, ctfxpaddim, ctfmask);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ctfmask)
            {
                if (phase_flipped)
                    dAij(Mctf, i, j) = fabs(dAij(ctfmask, i, j).real());
                else
                    dAij(Mctf, i, j) =      dAij(ctfmask, i, j).real();
            }

            //#define DEBUG
#ifdef  DEBUG
            {
                MetaData md1;
                size_t id;
                static int counter=0;

                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mctf)
                {
                    id=md1.addObject();
                    md1.setValue(MDL_ORDER,(int)(id-1),id);
                    md1.setValue(MDL_RESOLUTION_FRC,dAij(Mctf, 0, j),id );
                }
                std::stringstream ss;
                ss << counter;//add number to the stream
                md1.write(std::string("NEW/CTF_") +  ss.str());
                counter++;
            }
            ctf.write("new.ctfparam");
#endif
#undef DEBUG
            // Fill vectors
            ctfMD.setValue(MDL_ORDER,counter,__iter.objId);
            ctfMD.setValue(MDL_CTF_DEFOCUSA,avgdef,__iter.objId);
            mics_ctf2d.setSlice(0,Mctf,counter++);
        }
        else
        {
            std::cerr<<" Discard CTF "<<fnt_ctf<<" because of too large anisotropy"<<std::endl;
            ctfMD.removeObject(__iter.objId);
        }
        if (counter % c == 0 && verbose!=0)
            progress_bar(counter);
    }
    int enabled;
    // Precalculate denominator term of the Wiener filter
    int ii,jj;
    if (do_wiener)
    {

        if (verbose!=0)
            std::cout << "\nPrecalculate denominator term of the Wiener filter" <<std::endl;

        double sumimg = 0.;
        int bar=0;
        FOR_ALL_OBJECTS_IN_METADATA(ctfMD)
        {
            ctfMD.getValue(MDL_COUNT,count,__iter.objId);
            ctfMD.getValue(MDL_ORDER,counter,__iter.objId);
            sumimg += (double)count;
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mwien)
            {
                double d;
                if (i <= paddim/2 )
                    ii=i;
                else
                    ii=(paddim-i);

                if (j <= paddim/2)
                    jj=j;
                else
                    jj=(paddim-j);

                int dd = (int)(sqrt(ii*ii+jj*jj)+0.5);
                d = NZYX_ELEM(mics_ctf2d, counter, 0, 0, dd);
                //d = NZYX_ELEM(mics_ctf2d, counter, 0, i,j);
                ////d = NZYX_ELEM(mics_ctf2d, counter, 1, 1, i,j));
                dAij(Mwien,i,j) += count * d * d ;
            }

        }
        // Divide by sumimg (Wiener filter is for summing images, not averaging!)
        Mwien /= sumimg;
        // Add Wiener constant
        if (wiener_constant < 0.)
        {
            // Use Grigorieff's default for Wiener filter constant: 10% of average over all Mwien terms
            // Grigorieff JSB 157(1) (2006), pp 117-125
            wiener_constant = 0.1 * Mwien.computeAvg();
        }
        Mwien += wiener_constant;
        //#define DEBUG
#ifdef DEBUG

        {
            Image<double> save;
            save()=Mwien;
            save.write("vienerA.spi");
        }
#endif
#undef DEBUG

    }
    // Sort by average defocus
    sortedCtfMD.sort(ctfMD,MDL_CTF_DEFOCUSA,false);

}
// Check whether a CTF is anisotropic
bool ProgCtfGroup::isIsotropic(CTFDescription &ctf)
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
        ctf.precomputeValues(XX(freq), YY(freq));
        ctfp = ctf.CTF_at();
        ctf.precomputeValues(YY(freq), XX(freq));
        diff = ABS(ctfp - ctf.CTF_at());
        if (diff > max_error)
        {
            std::cerr<<" Anisotropy!"<<digres<<" "<<max_error<<" "<<diff<<" "<<ctfp
            <<" "<<ctf.CTF_at()<<std::endl;
            return false;
        }
    }
    return true;
}

// Do the actual work
void ProgCtfGroup::autoRun()
{
    double diff,mindiff;
    long int ctfMdSize=sortedCtfMD.size();
    init_progress_bar(ctfMdSize);
    int c = XMIPP_MAX(1, ctfMdSize / 60);
    int counter=0;
    bool newgroup;

    //size_t id;
    int orderOut, orderIn,defocusGroup;
    int groupNumber = 1;
    std::vector<size_t> vectorID;
    sortedCtfMD.findObjects(vectorID);
    sortedCtfMD.setValueCol(MDL_DEFGROUP,-1);
    //iterate
    std::vector<size_t>::iterator itOut;
    std::vector<size_t>::iterator itIn;
    size_t begin=vectorID.at(0);
    sortedCtfMD.setValue(MDL_DEFGROUP,groupNumber,begin);
    init_progress_bar(ctfMdSize);

    if (verbose!=0)
        std::cout << "\nCompute differences between CTFs" <<std::endl;

    for ( itOut=vectorID.begin()+1 ; itOut < vectorID.end(); itOut++ )
    {
        counter++;
        newgroup = true;
        mindiff = 99999.;
        sortedCtfMD.getValue(MDL_DEFGROUP,defocusGroup,*itOut);
        if(defocusGroup!=-1)
            continue;
        sortedCtfMD.getValue(MDL_ORDER,orderOut,*itOut);//index in mics_ctf2d array
        for ( itIn=vectorID.begin() ; itIn < itOut; itIn++ )
        {
            sortedCtfMD.getValue(MDL_ORDER,orderIn,*itIn);//index in mics_ctf2d array
            for (int iresol=0; iresol<=iresol_error; iresol++)
            {
                //NZYX_ELEM(mics_ctf2d, orderIn, 1, 1, iresol);
                diff = fabs( NZYX_ELEM(mics_ctf2d, orderIn,  0, 0, iresol) -
                             NZYX_ELEM(mics_ctf2d, orderOut, 0, 0, iresol) );
                if (diff > max_error)
                {
                    break;
                }
            }
            if (diff < max_error)
            {
                newgroup=false;
                break;
            }
        }
        if(newgroup)
        {
            groupNumber++;

        }
        sortedCtfMD.setValue(MDL_DEFGROUP,groupNumber,*itOut);
        if (counter % c == 0 && verbose!=0)
            progress_bar(counter);
    }

    progress_bar(ctfMdSize);

}
/////////////////////////////////////////////
void ProgCtfGroup::manualRun()
{
    MetaData DF;
    int groupNumber = 0;
    DF.read(fn_split);
    int counter=0;
    DF.setValueCol(MDL_DEFGROUP,-2);
    sortedCtfMD.setValueCol(MDL_DEFGROUP,-1);
    MetaData unionMD;
    DF.unionAll(sortedCtfMD);
    unionMD.sort(DF,MDL_CTF_DEFOCUSA);
    int n = unionMD.size();
    init_progress_bar(n);
    int c = XMIPP_MAX(1, n / 60);
    int defGroup;

    FOR_ALL_OBJECTS_IN_METADATA(unionMD)
    {
        unionMD.getValue(MDL_DEFGROUP,defGroup,__iter.objId);
        if(defGroup==-2)
            groupNumber++;
        else
            unionMD.setValue(MDL_DEFGROUP,groupNumber,__iter.objId);
        if (counter % c == 0)
            progress_bar(counter);
    }
    progress_bar(n);
    sortedCtfMD.importObjects(unionMD, MDValueNE(MDL_DEFGROUP, -2));
}

void ProgCtfGroup::writeOutputToDisc()
{
    //compute no of micrographs, no of images , minimum defocus ,maximum defocus, average defocus per ctf group
    MetaData ctfInfo,ctfImagesGroup;

    const AggregateOperation MyaggregateOperations[] =
        {
            AGGR_COUNT, AGGR_SUM, AGGR_MIN,         AGGR_MAX,          AGGR_AVG
        };
    std::vector<AggregateOperation> aggregateOperations(MyaggregateOperations,MyaggregateOperations+5);

    const MDLabel MyoperateLabels[]       =
        {
            MDL_COUNT,MDL_COUNT, MDL_CTF_DEFOCUSA, MDL_CTF_DEFOCUSA, MDL_CTF_DEFOCUSA
        };
    std::vector<MDLabel> operateLabels(MyoperateLabels,MyoperateLabels+5);

    const MDLabel MyresultLabels[]        =
        {
            MDL_DEFGROUP,MDL_COUNT, MDL_SUM,  MDL_MIN,          MDL_MAX,          MDL_AVG
        };
    std::vector<MDLabel> resultLabels(MyresultLabels,MyresultLabels+6);

    ctfInfo.aggregate(sortedCtfMD,aggregateOperations,operateLabels,resultLabels);
    ctfInfo.setComment("N. of micrographs, N. of particles, min defocus, max defocus and avg defocus");
    ctfInfo.write(fn_root+"Info.xmd");

    // make block-sel per image group
    MetaData ImagesMD, auxMetaData;
    ImagesMD.read(fn_ctfdat);
    ctfImagesGroup.join(ImagesMD,sortedCtfMD,MDL_CTFMODEL,INNER );

    unlink( (fn_root+"s_images.sel").c_str());
    FileName imagesInDefoculGroup;
    auxMetaData.setComment("images (particles) per defocus group, block name is defocusgroup No");
    for(int i=1;i<= ctfInfo.size(); i++)
    {
        auxMetaData.importObjects(ctfImagesGroup,MDValueEQ(MDL_DEFGROUP,i));
        imagesInDefoculGroup.assign( formatString("b%06d@%s_images.sel", i, fn_root.c_str()) );
        auxMetaData.write( imagesInDefoculGroup,APPEND);
    }

    //create average ctf<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check ctfs
    int olddefGroup,defGroup,order,count;
    double sumimg=0.;
    bool changeDefocusGroup;

    MultidimArray<double> ctf2D(Mwien.xdim,Mwien.ydim);
    Image<double> Ictf2D;
    Ictf2D.data.alias(ctf2D);

    olddefGroup=1;
    int jj,ii;
    bool savefile=false;
    FileName outFileName,outFileName2;
    outFileName = fn_root + ".ctf";
    if (format !="")
        outFileName=outFileName + ":" + format;
    std::cerr << "Saving CTFs" <<std::endl;
    FOR_ALL_OBJECTS_IN_METADATA(sortedCtfMD)
    {
        sortedCtfMD.getValue(MDL_DEFGROUP,defGroup,__iter.objId);
        sortedCtfMD.getValue(MDL_ORDER,order,__iter.objId);
        sortedCtfMD.getValue(MDL_COUNT,count,__iter.objId);
        if (defGroup != olddefGroup)
        {
            if (sumimg==0)
                continue;
            ctf2D /= sumimg;
            outFileName2.compose(olddefGroup,outFileName);
            Ictf2D.write(outFileName2);
            ctf2D.initZeros();
            olddefGroup=defGroup;
            sumimg=0.;
        }
        sumimg += (double)count;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ctf2D)
        {
            double d;
            if (i <= paddim/2 )
                ii=i;
            else
                ii=(paddim-i);

            if (j <= paddim/2)
                jj=j;
            else
                jj=(paddim-j);
            int dd = (int)(sqrt(ii*ii+jj*jj)+0.5);
            d = NZYX_ELEM(mics_ctf2d, order, 0, 0, dd);
            dAij(ctf2D,i,j) += count * d ;
        }
    }
    if (defGroup != olddefGroup)
    {
        if (sumimg!=0)
        {
            ctf2D /= sumimg;
            outFileName2.compose(defGroup,outFileName);
            Ictf2D.write(outFileName2);
            ctf2D.initZeros();
            olddefGroup=defGroup;
            sumimg=0.;
        }
    }

    //save ctf per group
    //save winer per group
    // save ctf and wienner profiles
#ifdef NEVER

    FileName fnt;
    MetaData SFo, DFo;
    Image<double> img;
    MultidimArray<double> Mavg;
    double sumw, avgdef, mindef, maxdef, split, oldmin=99.e99;
    int imic;
    std::ofstream fh, fh2, fh3, fh4;

    fh.open((fn_root+ "_groups.imgno").c_str(), std::ios::out);
    fh  << "# Number of images in each group \n";
    fh3.open((fn_root+ "_groups.defocus").c_str(), std::ios::out);
    fh3 << "# Defocus values for each group (avg, max & min) \n";

    for (int igroup=0; igroup < pointer_group2mic.size(); igroup++) //images names for all micrographies
    {
        SFo.clear();
        SFo.setComment("Defocus values to split into "+
                       integerToString(pointer_group2mic.size())+" ctf groups");
        Mavg.initZeros(paddim,paddim);
        sumw = 0.;
        avgdef = 0.;
        mindef = 99.e99;
        maxdef = -99.e99;
        for (int igmic=0; igmic < pointer_group2mic[igroup].size(); igmic++)//number images per micrograph
        {
            imic = pointer_group2mic[igroup][igmic];
            sumw += (double) mics_count[imic];
            // calculate (weighted) average Mctf
            //I need a loop here since mics_ctf2d this very likely will be 1D not 2Dtop

            Mavg += mics_count[imic] * mics_ctf2d[imic];
            // Calculate avg, min and max defocus values in this group
            avgdef += mics_count[imic] * mics_defocus[imic];
            mindef = XMIPP_MIN(mics_defocus[imic],mindef);
            maxdef = XMIPP_MAX(mics_defocus[imic],maxdef);
            // Fill SelFile
            size_t id;
            for (int iimg=0; iimg < mics_fnimgs[imic].size(); iimg++)
            {
                id = SFo.addObject();
                SFo.setValue(MDL_IMAGE,mics_fnimgs[imic][iimg], id);
            }
        }
        Mavg /= sumw;
        avgdef /= sumw;
        fnt.compose(fn_root+"_group",igroup+1,"");
        std::cerr<<" Group "<<fnt <<" contains "<< pointer_group2mic[igroup].size()<<" ctfs and "<<sumw<<" images and has average defocus "<<avgdef<<std::endl;

        // 1. write selfile with images per defocus group -> block now
        SFo.write(fnt+".sel");
        // 2. write average Mctf //skip these since we have profiles
        img() = Mavg;
        img.setWeight(sumw);
        img.write(fnt+".ctf");
        // 3. Output to file with numinsertber of images per group -> _groups.imgno, save the metadafile with this info
        fh << integerToString(igroup+1);
        fh.width(10);
        fh << floatToString(sumw)<<std::endl;
        // 4. Output to file with avgdef, mindef and maxdef per group save with previous file
        fh3 << integerToString(igroup+1);
        fh3.width(10);
        fh3 << floatToString(avgdef);
        fh3.width(10);
        fh3 << floatToString(maxdef);
        fh3.width(10);
        fh3 << floatToString(mindef)<<std::endl;
        // 5. Output to docfile for manual grouping, OK save al metadata
        if (oldmin < 9.e99)
        {
            split = (oldmin + maxdef) / 2.;
            DFo.addObject();
            DFo.setValue(MDL_CTF_DEFOCUSU,split);
        }
        oldmin = mindef;
        // 6. Write file with 1D profiles single metadata file
        fh2.open((fnt+".ctf_profiles").c_str(), std::ios::out);
        for (int i=0; i < paddim/2; i++)
        {
            fh2 << floatToString((double)i / (pixel_size*paddim));
            fh2.width(10);
            fh2 << floatToString(dAij(Mavg,i,0));
            fh2.width(10);
            for (int igmic=0; igmic < pointer_group2mic[igroup].size(); igmic++)
            {
                imic = pointer_group2mic[igroup][igmic];
                fh2 << floatToString(dAij(mics_ctf2d[imic],i,0));
                fh2.width(10);
            }
            fh2 << std::endl;
        }
        fh2.close();
        // 7. Write Wiener filter stack
        if (do_wiener)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mavg)
            {
                dAij(Mavg,i,j) /= dAij(Mwien,i,j);
            }
            img() = Mavg;
            img.setWeight(sumw);
            img.write(fnt+".wien");
            fh4.open((fnt+".wien_profile").c_str(), std::ios::out);
            for (int i=0; i < paddim/2; i++)
            {
                fh4 << floatToString((double)i / (pixel_size*paddim));
                fh4.width(10);
                fh4 << floatToString(dAij(Mavg,i,0));
                fh4.width(10);
                fh4 << std::endl;
            }
            fh4.close();
        }
    }

    // 3. Write file with number of images per group
    fh.close();
    // 4. Write file with avgdef, mindef and maxdef per group
    fh3.close();
    // 5. Write docfile with defocus values to split manually
    DFo.write(fn_root+"_groups_split.doc");
#endif
}

void ProgCtfGroup::run()
{
    produceSideInfo();

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
