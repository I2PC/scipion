/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2010)
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
#include <data/filters.h>
#include "xray_import.h"

// Read arguments ==========================================================
void ProgXrayImport::readParams()
{
    fnDirData = getParam("-data");
    fnRoot    = getParam("-oroot");
    fnDirFlat = getParam("-flat");
    cropSize  = getIntParam("-crop");
    thrNum    = getIntParam("-thr");
    fnBPMask  = getParam("--filterBadPixels");
    logFilt   = checkParam("--log");
}

// Show ====================================================================
void ProgXrayImport::show() const
{
    std::cout
    << "Input data directory       : " << fnDirData  << std::endl
    << "Input flatfield directory  : " << fnDirFlat  << std::endl
    << "Output rootname            : " << fnRoot     << std::endl
    << "Crop size                  : " << cropSize   << std::endl
    << "Number of threads          : " << thrNum     << std::endl
    ;
}

// usage ===================================================================
void ProgXrayImport::defineParams()
{
    addUsageLine("Import X-ray micrographs with .spe extension");
    addParamsLine("   -data <inputDataDirectory>       : Directory with SPE images, position files,");
    addParamsLine("                                    : and optionally, darkfields");
    addParamsLine("   -oroot <outputRootname>          : Rootname for output files");
    addParamsLine("  [-flat <inputFlatfieldDirectory=\"\">] : Directory with SPE images, position files,");
    addParamsLine("                                    : and optionally, darkfields");
    addParamsLine("  [-crop <size=0>]                  : Number of pixels to crop from each side");
    addParamsLine("  [-thr  <N=1>]                     : Number of threads");
    addParamsLine("   == Filters                                          ");
    addParamsLine("  [--filterBadPixels  <mask=\"\">]       : Apply a boundaries median filter to bad pixels given in mask.");
    addParamsLine("  alias -f;");
    addParamsLine("  [--log]       : Apply a log10 correction to normalized output images.");
    addParamsLine("  alias -l;");
}

// Really import ==========================================================
void ProgXrayImport::readAndCrop(const FileName &fn, Image<double> &I) const
{
    I.read(fn);
    FileName fnBase=fn.withoutExtension();
    std::ifstream fhPosition;
    fhPosition.open((fnBase+"-positions.txt").c_str());
    if (!fhPosition)
        REPORT_ERROR(ERR_IO_NOTEXIST,fnBase+"-positions.txt");
    int itemsFound=0;
    std::string line;
    std::vector<std::string> tokens;
    double tiltAngle;
    while (!fhPosition.eof())
    {
        getline(fhPosition,line);
        splitString(line," ",tokens);
        if (tokens[0]=="xm:sample:rx")
        {
            tiltAngle=textToFloat(tokens[1]);
            itemsFound++;
        }
        if (itemsFound==1)
            break;
    }
    fhPosition.close();
    if (itemsFound!=1)
        REPORT_ERROR(ERR_VALUE_EMPTY,(std::string)"Cannot find tilt angle in "+
                     fnBase+"-positions.txt");
    if (cropSize>0)
        I().window(cropSize,cropSize,
        (int)(YSIZE(I())-cropSize-1),(int)(XSIZE(I())-cropSize-1));
    STARTINGX(I())=STARTINGY(I())=0;
    I.setTilt(tiltAngle);
}

void ProgXrayImport::correct(Image<double> &I) const
{
    FileName fnBase=I.name().withoutExtension();
    std::ifstream fhPosition;
    fhPosition.open((fnBase+"-positions.txt").c_str());
    if (!fhPosition)
        REPORT_ERROR(ERR_IO_NOTEXIST,fnBase+"-positions.txt");
    int itemsFound=0;
    std::string line;
    std::vector<std::string> tokens;
    double currentBeam, exposureTime, slitWidth;
    while (!fhPosition.eof())
    {
        getline(fhPosition,line);
        splitString(line," ",tokens);
        if (tokens[0]=="sr:current")
        {
            currentBeam=textToFloat(tokens[1]);
            itemsFound++;
        }
        else if (tokens[0]=="xm:ccd:exp_time")
        {
            exposureTime=textToFloat(tokens[1]);
            itemsFound++;
        }
        else if (tokens[0]=="xm:mono:slitwidth")
        {
            slitWidth=textToFloat(tokens[1]);
            itemsFound++;
        }
        if (itemsFound==3)
            break;
    }
    fhPosition.close();
    if (itemsFound!=3)
        REPORT_ERROR(ERR_VALUE_EMPTY,(std::string)"Cannot find all parameters in "+
                     fnBase+"-positions.txt");
    I()/=currentBeam*exposureTime*slitWidth;
}

void ProgXrayImport::getDarkfield(const FileName &fnDir, Image<double> &IavgDark) const
{
    IavgDark.clear();
    std::vector<FileName> listDir;
    getdir(fnDir,listDir);
    bool found=false;
    for (int i=0; i<listDir.size(); i++)
        if (listDir[i]=="darkfields")
        {
            found=true;
            std::vector<FileName> listDirDark;
            getdir(fnDir+"/darkfields",listDirDark);
            int N=0;
            for (int j=0; j<listDirDark.size(); j++)
            {
                FileName extension=((FileName)listDirDark[j]).getExtension();
                if (extension!="spe")
                    continue;
                Image<double> Iaux;
                readAndCrop(fnDir+"/darkfields/"+listDirDark[j],Iaux);
                if (N==0)
                    IavgDark()=Iaux();
                else
                    IavgDark()+=Iaux();
                N++;
            }
            if (N==0)
                REPORT_ERROR(ERR_IO_NOTEXIST,"darkfields directory is empty");
            IavgDark()/=N;
            break;
        }
    if (!found)
        std::cerr << "Cannot find darkfields for "+fnDir << std::endl
        << "Execution continued without darkfield correction" << std::endl;
}

void ProgXrayImport::getFlatfield(const FileName &fnDir,
                                  Image<double> &Iavg) const
{
    // Get Darkfield
    std::cerr << "Getting darkfield from "+fnDir << " ..." << std::endl;
    Image<double> IavgDark;
    getDarkfield(fnDir, IavgDark);
    if (XSIZE(IavgDark())!=0)
        IavgDark.write(fnRoot+"_"+fnDir.removeDirectories()+"_darkfield.xmp");

    // Process the rest of the images
    std::vector<FileName> listDir;
    getdir(fnDir,listDir);
    int N=0;
    Image<double> Iaux;
    for (int i=0; i<listDir.size(); i++)
    {
        FileName extension=((FileName)listDir[i]).getExtension();
        if (extension!="spe")
            continue;
        FileName fnImg=fnDir+"/"+listDir[i];
        readAndCrop(fnImg,Iaux);
        if (XSIZE(IavgDark())!=0)
        {
            Iaux()-=IavgDark();
            forcePositive(Iaux());
        }
        correct(Iaux);

        if (N==0)
            Iavg()=Iaux();
        else
            Iavg()+=Iaux();
        N++;
    }
    if (N==0)
        REPORT_ERROR(ERR_IO_NOTEXIST,fnDir+" is empty");
    Iavg()/=N;

    /* Create a mask with zero valued pixels to apply boundaries median filter
     * to avoid dividing by zero when normalizing */
    MultidimArray<char> mask;
    Iavg().equal(0,mask);
    if (XSIZE(bpMask()) != 0)
        mask += bpMask();
    boundMedianFilter(Iavg(),mask);
}

void runThread(ThreadArgument &thArg)
{
    int thread_id = thArg.thread_id;
    ProgXrayImport * ptrProg= (ProgXrayImport *)thArg.workClass;

    MetaData localMD;
    Image<double> Iaux;
    FileName fnImg;
    long long int first = -1, last = -1;
    while (ptrProg->td->getTasks(first, last))
    {
        for (int i=first; i<=last; i++)
        {
            ptrProg->readAndCrop(ptrProg->filenames[i],Iaux);
            if (XSIZE(ptrProg->IavgDark())!=0)
            {
                Iaux()-=ptrProg->IavgDark();
                forcePositive(Iaux());
            }
            ptrProg->correct(Iaux);

            if (XSIZE(ptrProg->IavgFlat())!=0)
                Iaux()/=ptrProg->IavgFlat();

            // Assign median filter to zero valued pixels to avoid -inf when applying log10
            MultidimArray<char> mask;
            Iaux().equal(0,mask);
            if (XSIZE(ptrProg->bpMask()) != 0)
                mask += ptrProg->bpMask();
            boundMedianFilter(Iaux(), mask);

            if (ptrProg->logFilt)
                Iaux().selfLog10();

            fnImg.compose(i, ptrProg->fnRoot);
            fnImg = fnImg.addExtension("mrcs");

            localMD.addObject();
            localMD.setValue(MDL_IMAGE,fnImg);
            localMD.setValue(MDL_ANGLETILT,Iaux.tilt());
            Iaux.write(fnImg,i,true,WRITE_APPEND);
            if (thread_id==0)
                progress_bar(i);
        }
    }
    //Lock for update the total counter
    ptrProg->mutex.lock();
    ptrProg->MD.unionAll(localMD);
    ptrProg->mutex.unlock();
}

void ProgXrayImport::run()
{
    // Delete output stack if it exists
    if (exists((fnRoot+".mrcs").c_str()))
        unlink((fnRoot+".mrcs").c_str());

    // Reading bad pixels mask
    if (fnBPMask != "")
    {
        std::cerr << "Reading bad pixels mask from "+fnBPMask << "." << std::endl;
        bpMask.read(fnBPMask);
        if (cropSize>0)
            bpMask().window(cropSize,cropSize,
            (int)(YSIZE(bpMask())-cropSize-1),(int)(XSIZE(bpMask())-cropSize-1));
        STARTINGX(bpMask())=STARTINGY(bpMask())=0;
    }

    // Get the flatfield
    if (fnDirFlat!="")
    {
        std::cerr << "Getting flatfield from "+fnDirFlat << " ..." << std::endl;
        getFlatfield(fnDirFlat,IavgFlat);
        if (XSIZE(IavgFlat())!=0)
            IavgFlat.write(fnRoot+"_"+fnDirFlat.removeDirectories()+"_avg.xmp");
    }

    // Get the darkfield
    std::cerr << "Getting darkfield from "+fnDirData << " ..." << std::endl;
    getDarkfield(fnDirData, IavgDark);
    if (XSIZE(IavgDark())!=0)
        IavgDark.write(fnRoot+"_darkfield.xmp");

    // Count the number of images
    std::vector<FileName> listDir;
    getdir(fnDirData,listDir);
    int N=0;
    for (int i=0; i<listDir.size(); i++)
    {
        FileName extension=((FileName)listDir[i]).getExtension();
        if (extension!="spe")
            continue;
        FileName fnImg=fnDirData+"/"+listDir[i];
        filenames.push_back(fnImg);
        N++;
    }

    // Process images
    td=new ThreadTaskDistributor(filenames.size(),XMIPP_MAX(1,filenames.size()/30));
    tm=new ThreadManager(thrNum,this);
    std::cerr << "Getting data from " << fnDirData << " ...\n";
    init_progress_bar(filenames.size());
    tm->run(runThread);
    progress_bar(filenames.size());

    // Write Selfile and angles
    MetaData MDSorted;
    MDSorted.sort(MD,MDL_ANGLETILT);
    MDSorted.write(fnRoot+".sel");
    std::ofstream fhTlt;
    fhTlt.open((fnRoot+".tlt").c_str());
    if (!fhTlt)
        REPORT_ERROR(ERR_IO_NOWRITE,fnRoot+".tlt");
    FOR_ALL_OBJECTS_IN_METADATA(MDSorted)
    {
        double tilt;
        MDSorted.getValue(MDL_ANGLETILT,tilt);
        fhTlt << tilt << std::endl;
    }
    fhTlt.close();
}
