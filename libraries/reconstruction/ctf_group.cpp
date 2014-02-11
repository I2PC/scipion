
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
#include <data/xmipp_fft.h>
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
    do_wiener        = checkParam("--wiener");
    replaceSampling  = checkParam("--sampling_rate");
    if(replaceSampling)
        samplingRate       =  getDoubleParam("--sampling_rate");
    memory           = getDoubleParam("--memory");
    do1Dctf          = checkParam("--do1Dctf");
    wiener_constant  = getDoubleParam("--wc");
}

/* Show -------------------------------------------------------------------- */
void ProgCtfGroup::show()
{
    std::cout << "  Input ctfdat file       : "<< fn_ctfdat << std::endl;
    std::cout << "  Output rootname         : "<< fn_root << std::endl;
    if (pad > 1.)
    {
        std::cout << "  Padding factor          : "<< pad << std::endl;
    }
    if (do_discard_anisotropy)
    {
        std::cout << " -> Exclude anisotropic CTFs from the groups"<<std::endl;
    }
    if (do_auto)
    {
        std::cout << " -> Using automated mode for making groups"<<std::endl;
        std::cout << " -> With a maximum allowed error of "<<max_error
        <<" at "<<resol_error<<" dig. freq."<<std::endl;
    }
    else
    {
        std::cout << " -> Group based on defocus values in "<<fn_split<<std::endl;
    }
    if (phase_flipped)
    {
        std::cout << " -> Assume that data are PHASE FLIPPED"<<std::endl;
    }
    else
    {
        std::cout << " -> Assume that data are NOT PHASE FLIPPED"<<std::endl;
    }
    if (do_wiener)
    {
        std::cout << " -> Also calculate Wiener filters, with constant= "<<wiener_constant<<std::endl;
    }
    std::cout << "  Available memory (Gb)        : "<< memory << std::endl;
    if (do1Dctf)
    {
        std::cout << " -> compute CTF groups using 1D CTFs= "<<std::endl;
    }
    else
    {
        std::cout << " -> compute CTF groups using 2D CTFs= "<<std::endl;
    }
    std::cout << "----------------------------------------------------------"<<std::endl;
}

/* Usage ------------------------------------------------------------------- */
void ProgCtfGroup::defineParams()
{
    addUsageLine("Generate CTF (or defocus) groups from a single CTFdat file.");
    addUsageLine("+The automated mode groups all images with absolute difference in");
    addUsageLine("+CTF-values up to a given resolution below a given threshold.");
    addUsageLine("+For example the example bellow groups all images together with");
    addUsageLine("+absolute CTF differences smaller than 0.5 up to 15 Angstroms resolution).");
    addUsageLine("+A complementary manual mode allows to combine different groups or to split)");
    addUsageLine("+groups up even further.");
    addSeeAlsoLine("ctf_create_ctfdat");
    addExampleLine("Example of use: Sample using automated mode (resolution = 15 Ang.)",false);
    addExampleLine("   xmipp_ctf_group --ctfdat all_images_new.ctfdat -o CtfGroupsNew/ctfAuto   --wiener --wc -1 --pad 2 --phase_flipped --error 0.5 --resol 15 ");
    addExampleLine("Example of use: Sample using manual mode (after manual editing of ctf_group_split.doc)",false);
    addExampleLine("   xmipp_ctf_group --ctfdat all_images_new.ctfdat -o CtfGroupsNew/ctfManual --wiener --wc -1 --pad 2 --phase_flipped --split CtfGroupsNew/ctf_group_split.doc");


    //    addParamsLine("   -i <sel_file>              : Input selfile");
    addParamsLine("   --ctfdat <ctfdat_file>     : Input CTFdat file for all data");
    addParamsLine("   [-o <oext=\"ctf:stk\">]        : Output root name, you may force format ctf:mrc");
    addParamsLine("   [--pad <float=1>]          : Padding factor ");
    addParamsLine("   [--phase_flipped]          : Output filters for phase-flipped data");
    addParamsLine("   [--discard_anisotropy]     : Exclude anisotropic CTFs from groups");
    addParamsLine("   [--wiener]                 : Also calculate Wiener filters");
    addParamsLine("   [--sampling_rate <s>]       : This sampling rate overwrites the one in the ctf.param files");
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
    FileName  fnt_ctf,aux;
    CTFDescription ctf;
    //MetaData ctfdat,
    MetaData SF;
    MultidimArray<double> Mctf;
    MultidimArray<std::complex<double> >  ctfmask;
    size_t ydim, zdim, ndim;
    double avgdef;

    SF.read(fn_ctfdat);
    getImageSize(SF,dim,ydim,zdim,ndim);
    //set output format
    SF.getValue(MDL_IMAGE,aux,SF.firstObject());
    if(format=="")
        format=aux.getFileFormat();

    if ( dim != ydim )
        REPORT_ERROR(ERR_MULTIDIM_SIZE,"Only squared images are allowed!");

    paddim=xpaddim = ROUND(pad*dim);
    if(do1Dctf)
    {
        ypaddim=1;
        ctfxpaddim =  (size_t)(sqrt(2.) *  xpaddim + 1);
    }
    else
    {
        //ypaddim = xpaddim;
        //This is ready for the day in which we use anisotropic ctf
        ypaddim=1;
        //ctfxpaddim = xpaddim;
        ctfxpaddim =  (size_t)(sqrt(2.) *  xpaddim + 1);
    }
    Mctf.resize(ypaddim,ctfxpaddim);

    if (do_wiener)
    {
        Mwien.resize(paddim,paddim);
        Mwien.initZeros();
    }

    diff.resize(paddim,paddim);
    dd.resize(paddim,paddim);
    {
        double d;
        int ii,jj;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(diff)
        {
            if (i <= paddim/2 )
                ii=i;
            else
                ii=(paddim-i);

            if (j <= paddim/2)
                jj=j;
            else
                jj=(paddim-j);

            d      = sqrt(ii*ii+jj*jj);
            int idd = (int) d ;
            dAij(diff, i,j)=d-idd;
            dAij(dd, i,j)=idd;
        }
    }


    MetaData ctfMD;
    groupCTFMetaData(SF, ctfMD);

    int nCTFs = ctfMD.size();
    //how much memory do I need to store them
    double _sizeGb = (double) ypaddim * xpaddim * sizeof(double) * nCTFs /1073741824.;
    mmapOn = _sizeGb > memory;
    mics_ctf2d.setMmap(mmapOn);
    mics_ctf2d.resize(nCTFs,1,ypaddim,ctfxpaddim);

    int c = XMIPP_MAX(1, nCTFs / 60);
    init_progress_bar(nCTFs);

    //use this sampling instead of the one in the CTFparam file
    if(replaceSampling)
    {
        ctfMD.setValueCol(MDL_CTF_SAMPLING_RATE, samplingRate);
    }
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
    size_t count;
    size_t counter=0;
    if (verbose!=0)
        std::cout << "\nFill multiarray with ctfs" <<std::endl;
    FOR_ALL_OBJECTS_IN_METADATA(ctfMD)
    {
        ctf.readFromMetadataRow(ctfMD, __iter.objId);
        ctf.enable_CTF = true;
        ctf.enable_CTFnoise = false;
        ctf.produceSideInfo();
        if (pixel_size != ctf.Tm)
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         "Cannot mix CTFs with different sampling rates!");
        ctf.Tm /= sqrt(2.);
        if (!do_discard_anisotropy || isIsotropic(ctf))
        {
            avgdef = (ctf.DeltafU + ctf.DeltafV)/2.;
            ctf.DeltafU = avgdef;
            ctf.DeltafV = avgdef;
            ctf.generateCTF(ypaddim, ctfxpaddim, ctfmask);
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
                    md1.setValue(MDL_ORDER,(size_t)(id-1),id);
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
            std::cout<<" Discard CTF "<<fnt_ctf<<" because of too large anisotropy"<<std::endl;
            ctfMD.removeObject(__iter.objId);
        }
        if (counter % c == 0 && verbose!=0)
            progress_bar(counter);
    }
    // Precalculate denominator term of the Wiener filter
    if (do_wiener)
    {

        if (verbose!=0)
            std::cout << "\nPrecalculate denominator term of the Wiener filter" <<std::endl;

        double sumimg = 0.;
        double result;
        FOR_ALL_OBJECTS_IN_METADATA(ctfMD)
        {
            ctfMD.getValue(MDL_COUNT,count,__iter.objId);
            ctfMD.getValue(MDL_ORDER,counter,__iter.objId);
            double dCount = (double)count;
            sumimg += dCount;
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mwien)
            {
                //change DIRECT_N__X_ELEM by DIRECT_N_YX_ELEM if you want to process a 2D ctf
                result =         dAij(diff, i,j)  * DIRECT_N__X_ELEM(mics_ctf2d, counter, 0, 0, dAij(dd,i,j)+1  ) +
                                 (1.-dAij(diff, i,j)) * DIRECT_N__X_ELEM(mics_ctf2d, counter, 0, 0, dAij(dd,i,j));
                dAij(Mwien,i,j) += dCount * result *result;

            }
            //#define DEBUG
#ifdef DEBUG

            {
                std::cerr << "no_micro_per_Group order " << count
                << " " << counter
                << " " << NZYX_ELEM(mics_ctf2d, counter, 0, 0, 141)
                << " " << dAij(Mwien,100,100)
                << std::endl;
                Image<double> save;
                save()=Mwien;
                save.write("vienertempnew.spi");
                std::cout << "Press any key\n";
                char c;
                std::cin >> c;
            }
#endif
#undef DEBUG

        }
        // Divide by sumimg (Wiener filter is for summing images, not averaging!)
        //#define DEBUG
#ifdef DEBUG

        {
            Image<double> save;
            save()=Mwien;
            save.write("vienerB.spi");
        }
#endif
#undef DEBUG
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
        ctfp = ctf.getValueAt();
        ctf.precomputeValues(YY(freq), XX(freq));
        diff = ABS(ctfp - ctf.getValueAt());
        if (diff > max_error)
        {
            std::cout<<" Anisotropy!"<<digres<<" "<<max_error<<" "<<diff<<" "<<ctfp
            <<" "<<ctf.getValueAt()<<std::endl;
            return false;
        }
    }
    return true;
}

// Do the actual work
void ProgCtfGroup::autoRun()
{
    double diff=0.;
    long int ctfMdSize=sortedCtfMD.size();
    int c = XMIPP_MAX(1, ctfMdSize / 60);
    int counter=0;
    bool newgroup;

    //size_t id;
    size_t orderOut, orderIn;
    int defocusGroup;
    int groupNumber = 1;
    std::vector<size_t> vectorID;
    sortedCtfMD.findObjects(vectorID);
    sortedCtfMD.setValueCol(MDL_DEFGROUP,-1);
    //iterate
    std::vector<size_t>::iterator itOut;
    std::vector<size_t>::iterator itIn;
    size_t begin=vectorID.at(0);
    sortedCtfMD.setValue(MDL_DEFGROUP,groupNumber,begin);

    if (verbose!=0)
    {
        std::cout << "\nCompute differences between CTFs" <<std::endl;
        init_progress_bar(ctfMdSize);
    }

    for ( itOut=vectorID.begin()+1 ; itOut < vectorID.end(); itOut++ )
    {
        counter++;
        newgroup = true;
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
    int groupNumber = 1;
    DF.read(fn_split);
    int counter=0;
    DF.setValueCol(MDL_DEFGROUP,-2);
    sortedCtfMD.setValueCol(MDL_DEFGROUP,-1);
    //#define DEBUG
#ifdef DEBUG

    sortedCtfMD.write("sortedCtfMD1.xmd");
#endif
#undef DEBUG

    MetaData unionMD;
    DF.unionAll(sortedCtfMD);
    unionMD.sort(DF,MDL_CTF_DEFOCUSA,false);
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
    //#define DEBUG
#ifdef DEBUG

    unionMD.write("unionMD.xmd");
#endif
#undef DEBUG

    sortedCtfMD.importObjects(unionMD, MDValueNE(MDL_DEFGROUP, -2));
    //#define DEBUG
#ifdef DEBUG

    sortedCtfMD.write("sortedCtfMD3.xmd");
#endif
#undef DEBUG

}

void ProgCtfGroup::writeOutputToDisc()
{
    //(1) compute no of micrographs, no of images , minimum defocus ,maximum defocus, average defocus per ctf group
    MetaData ctfInfo,ctfImagesGroup,auxMetaData;

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
    ctfInfo.write("groups@"+fn_root+"Info.xmd");
    size_t numberDefGroups=ctfInfo.size();
    MetaData MD;
    size_t idctf = MD.addObject();
    MD.setValue(MDL_COUNT,numberDefGroups,idctf);
    MD.setColumnFormat(false);
    MD.write("numberGroups@"+fn_root+"Info.xmd",MD_APPEND);

    //(2)save auxiliary file for defocus split

    double maxDef,minDef;
    MDIterator it(ctfInfo);
    size_t id1,id2,id;
    auxMetaData.clear();
    auxMetaData.setComment(formatString("Defocus values to split into %lu ctf groups", ctfInfo.size()));
    id1=it.objId;
    while(it.moveNext())
    {
        id2=it.objId;
        ctfInfo.getValue(MDL_MIN,minDef,id1);
        ctfInfo.getValue(MDL_MAX,maxDef,id2);
        id1=id2;

        id=auxMetaData.addObject();
        auxMetaData.setValue(MDL_CTF_DEFOCUSA,(minDef+maxDef)/2.,id);
    }
    auxMetaData.write(fn_root+"_split.doc");


    //(3) make block-sel per image group
    MetaData ImagesMD;
    ImagesMD.read(fn_ctfdat);
    //
    if (ImagesMD.containsLabel(MDL_CTF_MODEL))
    {
        ctfImagesGroup.join1(ImagesMD, sortedCtfMD, MDL_CTF_MODEL, INNER );
    }
    else
    {
    	//ImagesMD.write("/tmp/ImagesMD");
    	//sortedCtfMD.write("/tmp/sortedCtfMD");
        ctfImagesGroup.join1(ImagesMD, sortedCtfMD, MDL_ITEM_ID, INNER);
    }

    //
    unlink( (fn_root+"s_images.sel").c_str());
    FileName imagesInDefoculGroup;
    auxMetaData.clear();
    auxMetaData.setComment("images (particles) per defocus group, block name is defocusgroup No");
    int ctfInfoSize;
    ctfInfoSize = (int) ctfInfo.size();
    for(int i=1;i<= ctfInfoSize; i++)
    {
        auxMetaData.importObjects(ctfImagesGroup,MDValueEQ(MDL_DEFGROUP,i));
        imagesInDefoculGroup.assign( formatString("ctfGroup%06d@%s_images.sel", i, fn_root.c_str()) );
        auxMetaData.write( imagesInDefoculGroup, i > 1 ? MD_APPEND : MD_OVERWRITE);
    }

    //(4)create average ctf
    int olddefGroup,defGroup;
    size_t order, count;
    double sumimg=0.;

    MultidimArray<double> ctf2D(Mwien.xdim,Mwien.ydim);
    Image<double> Ictf2D;
    Ictf2D.data.alias(ctf2D);

    olddefGroup=1;
    //defGroup=-1;
    FileName outFileNameCTF,outFileNameWIEN,outFileName;
    outFileNameCTF = fn_root + "_ctf."+format;
    outFileNameWIEN = fn_root + "_wien."+format;
    if (verbose!=0)
    {
        std::cout << "Saving CTF Images" <<std::endl;
    }
    FOR_ALL_OBJECTS_IN_METADATA(sortedCtfMD)
    {
        sortedCtfMD.getValue(MDL_DEFGROUP,defGroup,__iter.objId);
        sortedCtfMD.getValue(MDL_ORDER,order,__iter.objId);
        sortedCtfMD.getValue(MDL_COUNT,count,__iter.objId);
        double dCount = (double)count;

        if (defGroup != olddefGroup)
        {
            if (sumimg!=0)
            {
                ctf2D /= sumimg;
                outFileName.compose(olddefGroup,outFileNameCTF);
                //save CTF
                Ictf2D.write(outFileName);
                //save winer filter
                if (do_wiener)
                {
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ctf2D)
                    {
                        dAij(ctf2D,i,j) /= dAij(Mwien,i,j);
                    }
                    outFileName.compose(olddefGroup,outFileNameWIEN);
                    Ictf2D.write(outFileName);
                }

                ctf2D.initZeros();
                olddefGroup=defGroup;
                sumimg=0.;
            }
        }

        sumimg += dCount;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ctf2D)
        {
            double result;
            //interpolate ctf point from table
            //change DIRECT_N__X_ELEM by DIRECT_N_YX_ELEM if you want to process a 2D ctf

            result =     dAij(diff, i,j)      * DIRECT_N__X_ELEM(mics_ctf2d, order, 0, 0, dAij(dd,i,j)+1  ) +
                         (1.-dAij(diff, i,j)) * DIRECT_N__X_ELEM(mics_ctf2d, order, 0, 0, dAij(dd,i,j));
            dAij(ctf2D,i,j) += dCount * result ;
        }
    }
    //Save last CTF
    if (defGroup == olddefGroup)
    {
        if (sumimg!=0)
        {
            ctf2D /= sumimg;
            outFileName.compose(defGroup,outFileNameCTF);
            //save CTF
            Ictf2D.write(outFileName);
            //save winer filter
            if (do_wiener)
            {
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ctf2D)
                {
                    dAij(ctf2D,i,j) /= dAij(Mwien,i,j);
                }
                outFileName.compose(defGroup,outFileNameWIEN);
                Ictf2D.write(outFileName);
            }
        }
    }

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

    std::cout << " Done!" <<std::endl;
}
