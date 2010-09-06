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
#include "xray_import.h"
#include <sort.h>

// Read arguments ==========================================================
void Prog_xray_import_prm::read(int argc, char **argv)
{
    fnDirData = getParameter(argc, argv, "-data");
    fnRoot = getParameter(argc, argv, "-oroot");
    fnDirFlat = getParameter(argc, argv, "-flat","");
    cropSize = textToInteger(getParameter(argc, argv, "-crop", "0"));
}

// Show ====================================================================
std::ostream & operator << (std::ostream &out, const Prog_xray_import_prm &prm)
{
    out << "Input data directory       : " << prm.fnDirData  << std::endl
    << "Input flatfield directory  : " << prm.fnDirFlat  << std::endl
    << "Output rootname            : " << prm.fnRoot     << std::endl
    << "Crop size                  : " << prm.cropSize   << std::endl
    ;
    return out;
}

// usage ===================================================================
void Prog_xray_import_prm::usage() const
{
    std::cerr
    << "   -data <input data directory>      : Directory with SPE images, position files,\n"
    << "                                       and optionally, darkfields\n"
    << "  [-flat <input flatfield directory>]: Directory with SPE images, position files,\n"
    << "                                       and optionally, darkfields\n"
    << "   -oroot <output rootname>          : Rootname for output files\n"
    << "  [-crop <size=0>]                   : Number of pixels to crop from each side\n"
    ;
}

// Really import ==========================================================
void Prog_xray_import_prm::readAndCorrect(const FileName &fn, Image<double> &I)
{
    I.read(fn);
    FileName fnBase=fn.without_extension();
    std::ifstream fhPosition;
    fhPosition.open((fnBase+"-positions.txt").c_str());
    if (!fhPosition)
        REPORT_ERROR(ERR_IO_NOTEXIST,fnBase+"-positions.txt");
    int itemsFound=0;
    std::string line;
    std::vector<std::string> tokens;
    double currentBeam, exposureTime, slitWidth, tiltAngle;
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
        else if (tokens[0]=="xm:sample:rx")
        {
            tiltAngle=textToFloat(tokens[1]);
            itemsFound++;
        }
        if (itemsFound==4)
            break;
    }
    fhPosition.close();
    if (itemsFound!=4)
        REPORT_ERROR(ERR_VALUE_EMPTY,(std::string)"Cannot find all parameters in "+
                     fnBase+"-positions.txt");
    I().window(cropSize,cropSize,YSIZE(I())-cropSize-1,XSIZE(I())-cropSize-1);
    STARTINGX(I())=STARTINGY(I())=0;
    I()/=currentBeam*exposureTime*slitWidth;
    I.setTilt(tiltAngle);
}

void Prog_xray_import_prm::getDarkfield(const FileName &fnDir, Image<double> &IavgDark)
{
    IavgDark.clear();
    std::vector<std::string> listDir;
    getdir(fnDir,listDir);
    bool found=false;
    for (int i=0; i<listDir.size(); i++)
        if (listDir[i]=="darkfields")
        {
            found=true;
            std::vector<std::string> listDirDark;
            getdir(fnDir+"/darkfields",listDirDark);
            int N=0;
            for (int j=0; j<listDirDark.size(); j++)
            {
                FileName extension=((FileName)listDirDark[j]).get_extension();
                if (extension!="spe")
                    continue;
                Image<double> Iaux;
                readAndCorrect(fnDir+"/darkfields/"+listDirDark[j],Iaux);
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

void Prog_xray_import_prm::getCorrectedAvgOfDirectory(const FileName &fnDir,
        Image<double> &Iavg)
{
    // Get Darkfield
    std::cerr << "Getting darkfield from "+fnDir << " ..." << std::endl;
    Image<double> IavgDark;
    getDarkfield(fnDir, IavgDark);
    if (XSIZE(IavgDark())!=0)
        IavgDark.write(fnRoot+"_"+fnDir+"_darkfield.xmp");

    // Process the rest of the images
    std::vector<std::string> listDir;
    getdir(fnDir,listDir);
    int N=0;
    Image<double> Iaux;
    for (int i=0; i<listDir.size(); i++)
    {
        FileName extension=((FileName)listDir[i]).get_extension();
        if (extension!="spe")
            continue;
        FileName fnImg=fnDir+"/"+listDir[i];
        std::cout << "Processing " << fnImg << " ..." << std::endl;
        readAndCorrect(fnImg,Iaux);
        if (XSIZE(IavgDark())!=0) {
            Iaux()-=IavgDark();
            forcePositive(Iaux());
        }

        if (N==0)
            Iavg()=Iaux();
        else
            Iavg()+=Iaux();
        N++;
    }
    if (N==0)
        REPORT_ERROR(ERR_IO_NOTEXIST,fnDir+" is empty");
    Iavg()/=N;
}

void Prog_xray_import_prm::run()
{
    // Get the flatfield
    Image<double> IavgFlat;
    if (fnDirFlat!="")
    {
        std::cerr << "Getting flatfield from "+fnDirFlat << " ..." << std::endl;
        getCorrectedAvgOfDirectory(fnDirFlat,IavgFlat);
        if (XSIZE(IavgFlat())!=0)
            IavgFlat.write(fnRoot+"_"+fnDirFlat+"_avg.xmp");
    }

    // Get the darkfield
    std::cerr << "Getting darkfield from "+fnDirData << " ..." << std::endl;
    Image<double> IavgDark;
    getDarkfield(fnDirData, IavgDark);
    if (XSIZE(IavgDark())!=0)
        IavgDark.write(fnRoot+"_"+fnDirData+"_darkfield.xmp");

    // Count the number of images
    std::vector<std::string> listDir;
    std::vector<FileName> filenames;
    getdir(fnDirData,listDir);
    int N=0;
    Image<double> Iaux;
    for (int i=0; i<listDir.size(); i++)
    {
        FileName extension=((FileName)listDir[i]).get_extension();
        if (extension!="spe")
            continue;
        FileName fnImg=fnDirData+"/"+listDir[i];
        filenames.push_back(fnImg);
        N++;
    }
    std::sort(filenames.begin(),filenames.end());

    // Process images
    for (int i=0; i<filenames.size(); i++)
    {
        std::cout << "Processing " << filenames[i] << " ..." << std::endl;
        readAndCorrect(filenames[i],Iaux);
        if (XSIZE(IavgDark())!=0)
        {
            Iaux()-=IavgDark();
            forcePositive(Iaux());
        }
        if (XSIZE(IavgFlat())!=0)
            Iaux()/=IavgFlat();
    }
}
