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
    fnDirData = getParam("--data");
    fnRoot    = getParam("--oroot");
    fnDirFlat = getParam("--flat");
    cropSize  = getIntParam("--crop");
    thrNum    = getIntParam("--thr");
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
    addUsageLine("Imports X-ray micrographs with .spe extension");
    addUsageLine(" ");
    addUsageLine("+This program converts a single-axis tilt series coming from a X-ray");
    addUsageLine("+microscope in SPE format into a MRC stack. Several corrections are");
    addUsageLine("+applied during this conversion including correction by the beam");
    addUsageLine("+current, the exposure time, and the slitwidth. The flatfield and a");
    addUsageLine("+possible darkfield are also corrected. The exact formula applied is");
    addUsageLine("+ ");
    addUsageLine("+                       (I-darkfield)/(expTime*beamCurrent*slitWidth)",true);
    addUsageLine("+Icorrected=log10 -----------------------------------------------------------",true);
    addUsageLine("+                 Avg{(Iflatfield-darkfield)/(expTime*beamCurrent*slitWidth)}",true);
    addUsageLine("+ ");
    addUsageLine("+If the darkfield is not available, it is not used for the correction.");
    addUsageLine("+The flatfields and the the tilt series may be in any directory. If the");
    addUsageLine("+darkfield are available, they must be within the flatfield and the tilt");
    addUsageLine("+series under a directory named darkfields (in small letters). For each");
    addUsageLine("+SPE file there must be a positions file. For instance, the file myfile.spe");
    addUsageLine("+must have in the same directory the file myfile-positions.txt. The");
    addUsageLine("+exposure time, the beam current and the slit width are taken from this file.");
    addUsageLine("+ ");
    addUsageLine("+The exposure time is the parameter xm:ccd:exp_time, the beam current");
    addUsageLine("+is sr:current, and the slit width is xm:mono:slitwidth. The tilt angle");
    addUsageLine("+xm:sample:rx.");
    addParamsLine("   --data <inputDataDirectory>       : Directory with SPE images, position files,");
    addParamsLine("                                     : and optionally, darkfields (in a directory");
    addParamsLine("                                     : called darkfields)");
    addParamsLine("   --oroot <outputRootname>          : Rootname for output files");
    addParamsLine("                                     :+ The following files are created:");
    addParamsLine("                                     :+ rootname_darkfield.xmp: darkfield for the tilt series");
    addParamsLine("                                     :+ rootname_flatfields_darkfield.xmp: darkfield of the flatfields");
    addParamsLine("                                     :+ rootname_flatfields_avg.xmp: average flatfield corrected by its darkfield");
    addParamsLine("                                     :+ rootname.mrcs: MRC stack with tilt series fully corrected as described above");
    addParamsLine("                                     :+ rootname.sel: selfile with the images in the stack");
    addParamsLine("                                     :+ rootname.tlt: list of angles");
    addParamsLine("  [--flat <inputFlatfieldDirectory=\"\">] : Directory with SPE images, position files,");
    addParamsLine("                                     : and optionally, darkfields (in a directory");
    addParamsLine("                                     : called darkfields)");
    addParamsLine("  [--crop <size=0>]                  : Number of pixels to crop from each side");
    addParamsLine("                                     : This is used to avoid some black pixels of some cameras");
    addParamsLine("  [--thr  <N=1>]                     : Number of threads");
    addParamsLine("   == Filters                                          ");
    addParamsLine("  [--filterBadPixels  <mask=\"\">]        : Apply a boundaries median filter to bad pixels given in mask.");
    addParamsLine("  alias -f;");
    addParamsLine("  [--log]                            : Apply a log10 correction to normalized output images.");
    addParamsLine("  alias -l;");
    addExampleLine("The most standard call is",false);
    addExampleLine("xmipp_xray_import --data 10s --flat flatfields --oroot ProcessedData/img --crop 7");
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
        I().selfWindow(cropSize,cropSize,
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
    I()*=1.0/(currentBeam*exposureTime*slitWidth);
}

void ProgXrayImport::getDarkfield(const FileName &fnDir, Image<double> &IavgDark) const
{
    IavgDark.clear();
    std::vector<FileName> listDir;
    fnDir.getFiles(listDir);
    bool found = false;

    for (size_t i=0; i<listDir.size(); i++)
        if (listDir[i]=="darkfields")
        {
            found=true;
            std::vector<FileName> listDirDark;
            FileName(fnDir+"/darkfields").getFiles(listDirDark);
            int N = 0;

            for (size_t j=0; j<listDirDark.size(); j++)
            {
                if (!listDirDark[j].hasImageExtension())
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
            IavgDark()*=1.0/N;
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
    fnDir.getFiles(listDir);
    int N=0;
    Image<double> Iaux;
    for (size_t i=0; i<listDir.size(); i++)
    {
        if (!listDir[i].hasImageExtension())
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
    Iavg()*=1.0/N;

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
    size_t first = 0, last = 0;
    MultidimArray<char> mask;

    while (ptrProg->td->getTasks(first, last))
    {
        for (size_t i=first; i<=last; i++)
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
            Iaux().equal(0,mask);
            if (XSIZE(ptrProg->bpMask()) != 0)
                mask += ptrProg->bpMask();
            boundMedianFilter(Iaux(), mask);

            if (ptrProg->logFilt)
                Iaux().selfLog10();

            fnImg.compose(i+1, ptrProg->fnOut);

            size_t objId = localMD.addObject();
            localMD.setValue(MDL_IMAGE,fnImg,objId);
            localMD.setValue(MDL_ANGLE_TILT,Iaux.tilt(),objId);
            Iaux.write(fnImg);
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
    fnOut = fnRoot + ".mrcs";
    fnOut.deleteFile();

    // Reading bad pixels mask
    if (fnBPMask != "")
    {
        std::cerr << "Reading bad pixels mask from "+fnBPMask << "." << std::endl;
        bpMask.read(fnBPMask);
        if (cropSize>0)
            bpMask().selfWindow(cropSize,cropSize,
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
    fnDirData.getFiles(listDir);
    int N=0;
    for (size_t i=0; i<listDir.size(); i++)
    {
        if (!listDir[i].hasImageExtension())
            continue;
        FileName fnImg=fnDirData+"/"+listDir[i];
        filenames.push_back(fnImg);
        N++;
    }

    // Create empty output stack file
    size_t Xdim, Ydim, Zdim, Ndim;
    getImageSizeFromFilename(filenames[0], Xdim, Ydim, Zdim, Ndim);
    createEmptyFile(fnOut, Xdim-2*cropSize, Ydim-2*cropSize, 1, filenames.size());

    // Process images
    td = new ThreadTaskDistributor(filenames.size(),XMIPP_MAX(1,filenames.size()/30));
    tm = new ThreadManager(thrNum,this);
    std::cerr << "Getting data from " << fnDirData << " ...\n";
    init_progress_bar(filenames.size());
    tm->run(runThread);
    progress_bar(filenames.size());

    // Write Selfile and angles
    MetaData MDSorted;
    MDSorted.sort(MD,MDL_ANGLE_TILT);
    MDSorted.write(fnRoot+".sel");
    std::ofstream fhTlt;
    fhTlt.open((fnRoot+".tlt").c_str());
    if (!fhTlt)
        REPORT_ERROR(ERR_IO_NOWRITE,fnRoot+".tlt");
    FOR_ALL_OBJECTS_IN_METADATA(MDSorted)
    {
        double tilt;
        MDSorted.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        fhTlt << tilt << std::endl;
    }
    fhTlt.close();
    delete td;
    delete tm;
}
