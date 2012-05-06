/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

#include <data/metadata_extension.h>
#include <data/xmipp_program.h>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>

typedef enum
{
    EMX_COORDINATES,      //particle picking
    EMX_ALIGNMENT,        //aligment info
    EMX_CLASS,            //classification info
    EMX_CTFMICROGRAPH,    //CTF per micrograph
    EMX_CTFPARTICLE       //CTF per particle
} EMXConversion;


class ProgMetadataConvert22Emx: public XmippProgram
{
private:
    FileName fn_in, fn_out;
    MetaData mdParticle;
    MetaData mdProcessedParticle;
    MetaData mdMicrograph;
    std::map<String, EMXConversion> string2EMXconversion;
    std::map<EMXConversion, String> EMXconversion2string;
    EMXConversion conversionType;
    bool toemx;
    bool toxmipp;
    //tmpname is required just in case /dev/stdout is selected
    char *tmpname;
    String comment;
public:

    /* EMX standard has label with dots, xmipp does not support it
     * so first a temporary file is created changing "." by "____"
     */
    void restoreDots(void)
    {
        std::ifstream ifile(tmpname,std::ios::binary);
        ifile.seekg(0,std::ios_base::end);
        long s=ifile.tellg();
        char *buffer=new char[s];
        ifile.seekg(0);
        ifile.read(buffer,s);
        ifile.close();
        std::string txt(buffer,s);
        delete[] buffer;
        size_t off=0;
        while ((off=txt.find("____",off))!=std::string::npos)
            txt.replace(off,sizeof("____")-1,".");
        std::ofstream ofile(fn_out.c_str());
        ofile.write(txt.c_str(),txt.size());
        ofile.close();
    }

    void removeDots(void)
    {

        std::ifstream ifile(fn_in.c_str(),std::ios::binary);
        ifile.seekg(0,std::ios_base::end);
        long s=ifile.tellg();
        char *buffer=new char[s];
        ifile.seekg(0);
        ifile.read(buffer,s);
        ifile.close();
        std::string txt(buffer,s);
        delete[] buffer;
        size_t off=0;
        while ((off=txt.find("_emx_particle.",off))!=std::string::npos)
            txt.replace(off,sizeof("_emx_particle.")-1,"_emx_particle____");
        off=0;
        while ((off=txt.find("_emx_micrograph.",off))!=std::string::npos)
            txt.replace(off,sizeof("_emx_micrograph.")-1,"_emx_micrograph____");
        off=0;
        while ((off=txt.find("_emx_class.",off))!=std::string::npos)
            txt.replace(off,sizeof("_emx_class.")-1,"_emx_class____");
        std::ofstream ofile(tmpname);
        ofile.write(txt.c_str(),txt.size());
        ofile.close();
    }

    ProgMetadataConvert22Emx()
    {
        toemx=false;
        toxmipp=false;

        string2EMXconversion["coordinates"]   = EMX_COORDINATES;
        string2EMXconversion["alignment"]     = EMX_ALIGNMENT;
        string2EMXconversion["class"]         = EMX_CLASS;
        string2EMXconversion["ctfMicrograph"] = EMX_CTFMICROGRAPH;
        string2EMXconversion["ctfParticle"]   = EMX_CTFPARTICLE;

        EMXconversion2string[EMX_COORDINATES]   = (String)"coordinates";
        EMXconversion2string[EMX_ALIGNMENT]     = (String)"alignment";
        EMXconversion2string[EMX_CLASS]         = (String)"class";
        EMXconversion2string[EMX_CTFMICROGRAPH] = (String)"ctfMicrograph";
        EMXconversion2string[EMX_CTFPARTICLE]   = (String)"ctfParticle";
        //temporary file name, since emx standard uses dots
        //and xmipp metadata does not support them I change them to ____
        tmpname = strdup("/tmp/EMXtmpfileXXXXXX");
        mkstemp(tmpname);
    }

    ~ProgMetadataConvert22Emx()
    {
        unlink(tmpname);
        free(tmpname);
    }
protected:
    void defineParams()
    {
        addUsageLine("Convert to and from EMX metadata files");
        addParamsLine("  -i <text_file>                   :Input metadata file either xmipp or emx      ");
        addParamsLine("     alias --input;");
        addParamsLine("  [-o <output_metadata=\"/dev/stdout\">]  :Output metadata file, default screen");
        addParamsLine("  --conversion_type <mode>            :merge the imported metadata to an existing one");
        addParamsLine("         where <mode>");
        addParamsLine("             coordinates   : particle picking");
        addParamsLine("             alignment     : aligment info");
        addParamsLine("             class         : classification info");
        addParamsLine("             ctfMicrograph : CTF per micrograph");
        addParamsLine("             ctfParticle   : CTF per particle");
        addParamsLine("     alias -t;");
        addKeywords("convert emx");
        addExampleLine("Convert coordinates file from EMX to Xmipp", false);
        addExampleLine("xmipp_emx_convert -i particlePicking.emx -o particlePicking.xmd -t coordinates");

    }

    void readParams()
    {
        fn_in  = getParam("-i");
        fn_out = getParam("-o");
        conversionType = string2EMXconversion[getParam("-t")];

        setMetadataVersion("XMIPP_STAR_1");
        if (fn_in.isStar1(false))
        {
            toemx=true;
            comment = "##########################################################################\n"
                      "#               EMX Exchange file \n"
                      "# \n"
                      "#  Information on this file format is available at \n"
                      "#  http://i2pc.cnb.csic.es/emx\n"
                      "##########################################################################\n";
        }
        else
        {
            setMetadataVersion("EMX1.0");
            if (fn_in.isStar1(false))
            {
                toxmipp=true;
                comment = "##########################################################################\n"
                          "#               EMX Exchange file \n"
                          "# \n"
                          "# this metadata file has been imported from a EMX metadata file"
                          "#  Information on emx file format is available at \n"
                          "#  http://i2pc.cnb.csic.es/emx\n"
                          "##########################################################################\n";
            }
        }
        if(toemx==toxmipp)//both false
            REPORT_ERROR(ERR_IO_NOTFILE,(std::string)"File ");
    }

    void show()
    {
        std::cout << "\n\nxmipp_convert_metadata_22_emx " <<std::endl;
        std::cerr << "In file         : "  << fn_in << std::endl;
        std::cerr << "Out file        : " << fn_out << std::endl;
        std::cerr << "Conversion Type : " << EMXconversion2string[conversionType] << std::endl;
        if(toxmipp)
            std::cerr << "Converting from emx to xmipp:" << std::endl;
        else if (toemx)
            std::cerr << "Converting from xmipp to emx:" << std::endl;
    }
public:

    void convertXmipp2EmxCoordinates(void)
    {
        // Loop through all blocks
        StringVector blockList;
        MetaData mdCoordinateXmipp;
        MetaData mdCoordinateEMX;
        FileName fnCoordinateXmipp;
        FileName micrographName;
        MDRow rowIn, rowOut;
        int x,y;
        String particleName;

        getBlocksInMetaDataFile(fn_in,blockList);
        for (StringVector::iterator it= blockList.begin();
             it!=blockList.end(); it++)
        {
            micrographName=*it;
            fnCoordinateXmipp.compose(*it, fn_in);
            mdCoordinateXmipp.read(fnCoordinateXmipp);
            if (!mdCoordinateXmipp.containsLabel(MDL_XINT)  )
                REPORT_ERROR(ERR_MD_BADLABEL,(String)"Label: " + MDL::label2Str(MDL_XINT) + "missing.");
            if (!mdCoordinateXmipp.containsLabel(MDL_YINT)  )
                REPORT_ERROR(ERR_MD_BADLABEL,(String)"Label: " + MDL::label2Str(MDL_YINT) + "missing.");
            FOR_ALL_OBJECTS_IN_METADATA(mdCoordinateXmipp)
            {
                mdCoordinateXmipp.getRow(rowIn, __iter.objId);
                rowIn.getValue(MDL_XINT,x);
                rowIn.getValue(MDL_YINT,y);
                particleName=formatString("%s_%04d",micrographName.c_str(),__iter.objId);
                rowOut.setValue(MDL_EMX_PARTICLE_URL,particleName);
                rowOut.setValue(MDL_EMX_MICROGRAPH_URL,micrographName);
                rowOut.setValue(MDL_EMX_PARTICLE_COORDINATE_X,(double)x);
                rowOut.setValue(MDL_EMX_PARTICLE_COORDINATE_Y,(double)y);
                mdCoordinateEMX.addRow(rowOut);
            }
        }
        FileName tmpFn;
        tmpFn.compose("particle",tmpname);
        setMetadataVersion("EMX1.0");
        mdCoordinateEMX.setComment(comment);
        mdCoordinateEMX.write(tmpFn);
    }
    void convertXmipp2EmxCTFMicrograph(void)
    {
        std::cerr << "convertXmipp2EmxCTFMicrograph" <<std::endl;
        MetaData mdMicrographXmipp;
        MetaData mdMicrographEMX;
        String micrographXmipp;
        MDRow  rowIn, rowOut;
        double samplingRate;
        double defocusU;
        double defocusV;
        double defocusAngle;
        double voltage;
        double sphericalAberration;
        double Q0;
        mdMicrographXmipp.read(fn_in);
        FOR_ALL_OBJECTS_IN_METADATA(mdMicrographXmipp)
        {
            mdMicrographXmipp.getRow(rowIn, __iter.objId);
            rowIn.getValue(MDL_CTF_SAMPLING_RATE,samplingRate);
            rowIn.getValue(MDL_CTF_DEFOCUSU,defocusU);
            rowIn.getValue(MDL_CTF_DEFOCUSV,defocusV);
            rowIn.getValue(MDL_CTF_DEFOCUS_ANGLE,defocusAngle);
            rowIn.getValue(MDL_CTF_VOLTAGE,voltage);
            rowIn.getValue(MDL_CTF_CS,sphericalAberration);
            rowIn.getValue(MDL_CTF_Q0,Q0);
            rowIn.getValue(MDL_IMAGE,micrographXmipp);

            rowOut.setValue(MDL_EMX_MICROGRAPH_SAMPLING,samplingRate);
            rowOut.setValue(MDL_EMX_MICROGRAPH_DEFOCUSU, defocusU);
            rowOut.setValue(MDL_EMX_MICROGRAPH_DEFOCUSV, defocusV);
            rowOut.setValue(MDL_EMX_MICROGRAPH_ASTIGMATISM_ANGLE,defocusAngle);
            rowOut.setValue(MDL_EMX_MICROGRAPH_VOLTAGE, voltage);
            rowOut.setValue(MDL_EMX_MICROGRAPH_CS, sphericalAberration);
            rowOut.setValue(MDL_EMX_MICROGRAPH_AMPLITUDE_CONTRAST,Q0);
            rowOut.setValue(MDL_EMX_MICROGRAPH_URL,micrographXmipp);

            mdMicrographEMX.addRow(rowOut);
        }
        FileName tmpFn;
        tmpFn.compose("micrograph",tmpname);
        setMetadataVersion("EMX1.0");
        mdMicrographEMX.setComment(comment);
        mdMicrographEMX.write(tmpFn);
    }
    void convertXmipp2EmxClass(void)
    {
        MetaData mdClassXmipp;
        MetaData mdClassEMX;
        //        FileName fnCoordinateXmipp;
        //        FileName micrographName;
        MDRow rowIn, rowOut;
        int ref;
        String particleName;
        //
        mdClassXmipp.read(fn_in);
        std::stringstream out;
        FOR_ALL_OBJECTS_IN_METADATA(mdClassXmipp)
        {
            mdClassXmipp.getRow(rowIn, __iter.objId);
            rowIn.getValue(MDL_IMAGE,particleName);
            rowIn.getValue(MDL_REF,ref);
            rowOut.setValue(MDL_EMX_PARTICLE_URL,particleName);
            rowOut.setValue(MDL_EMX_PARTICLE_CLASS_ID, formatString("%04d",ref));
            mdClassEMX.addRow(rowOut);
        }
        FileName tmpFn;
        tmpFn.compose("processedParticle",tmpname);
        setMetadataVersion("EMX1.0");
        mdClassEMX.setComment(comment);
        mdClassEMX.write(tmpFn);
    }

    void convertEmx2XmippCoordinates(void)
    {
        MetaData mdCoordinateEMX;
        MetaData mdMicAggregate;
        MetaData mdSingleMicrograph;
        MetaData mdCoordinateXmipp;
        FileName micrographName;
        FileName particleName;
        FileName blockOut;
        MDRow rowIn, rowOut;
        double x,y;

        mdCoordinateEMX.read(tmpname);
        const MDLabel MyLabels[]       =
            {
                MDL_EMX_PARTICLE_COORDINATE_X,
                MDL_EMX_PARTICLE_COORDINATE_Y
            };
        std::vector<MDLabel> myVector(MyLabels,MyLabels+8);
        for (std::vector<MDLabel>::iterator it = myVector.begin(); it != myVector.end(); ++it)
        {
            if (!mdCoordinateEMX.containsLabel(*it) )
                REPORT_ERROR(ERR_MD_BADLABEL,
                             (String)"Label: " + MDL::label2Str(*it) + " missing.");
        }

        mdMicAggregate.aggregate(mdCoordinateEMX,AGGR_COUNT,MDL_EMX_MICROGRAPH_URL,
                                 MDL_UNDEFINED,MDL_COUNT);
        //must delete, append multiple blocks
        unlink(fn_out.c_str());
        FOR_ALL_OBJECTS_IN_METADATA(mdMicAggregate)
        {
            mdMicAggregate.getValue(MDL_EMX_MICROGRAPH_URL,
                                    micrographName,__iter.objId);
            mdSingleMicrograph.clear();
            MDValueEQ eq(MDL_EMX_MICROGRAPH_URL,micrographName);
            mdSingleMicrograph.importObjects(mdCoordinateEMX, eq);
            mdCoordinateXmipp.clear();
            FOR_ALL_OBJECTS_IN_METADATA(mdSingleMicrograph)
            {
                mdSingleMicrograph.getRow(rowIn, __iter.objId);
                rowIn.getValue(MDL_EMX_PARTICLE_COORDINATE_X,x);
                rowIn.getValue(MDL_EMX_PARTICLE_COORDINATE_Y,y);
                rowIn.getValue(MDL_EMX_PARTICLE_URL,particleName);
                rowOut.setValue(MDL_XINT,ROUND(x));
                rowOut.setValue(MDL_YINT,ROUND(y));
                rowOut.setValue(MDL_IMAGE,particleName);
                mdCoordinateXmipp.addRow(rowOut);
            }
            blockOut.compose(micrographName,fn_out);
            setMetadataVersion("XMIPP_STAR_1");
            mdCoordinateXmipp.setComment(comment);
            mdCoordinateXmipp.write(blockOut,MD_APPEND);
        }
    }

    void convertEmx2XmippCTFMicrograph(void)
    {
        MetaData mdCTFMicrograph;
        MetaData mdCTFMicrographXmipp;
        String micrographXmipp;
        MDRow  rowIn, rowOut;
        double samplingRate;
        double defocusU;
        double defocusV;
        double defocusAngle;
        double voltage;
        double sphericalAberration;
        double Q0;

        const MDLabel MyLabels[]       =
            {
                MDL_EMX_MICROGRAPH_SAMPLING,
                MDL_EMX_MICROGRAPH_DEFOCUSU,
                MDL_EMX_MICROGRAPH_DEFOCUSV,
                MDL_EMX_MICROGRAPH_ASTIGMATISM_ANGLE,
                MDL_EMX_MICROGRAPH_VOLTAGE,
                MDL_EMX_MICROGRAPH_CS,
                MDL_EMX_MICROGRAPH_AMPLITUDE_CONTRAST,
                MDL_EMX_MICROGRAPH_URL
            };
        std::vector<MDLabel> myVector(MyLabels,MyLabels+8);
        mdCTFMicrograph.read(tmpname);
        for (std::vector<MDLabel>::iterator it = myVector.begin(); it != myVector.end(); ++it)
        {
            if (!mdCTFMicrograph.containsLabel(*it) )
                REPORT_ERROR(ERR_MD_BADLABEL,
                             (String)"Label: " + MDL::label2Str(*it) + " missing.");
        }
        FOR_ALL_OBJECTS_IN_METADATA(mdCTFMicrograph)
        {
            mdCTFMicrograph.getRow(rowIn, __iter.objId);
            rowIn.getValue(MDL_EMX_MICROGRAPH_URL,micrographXmipp);
            rowIn.getValue(MDL_EMX_MICROGRAPH_SAMPLING,samplingRate);
            rowIn.getValue(MDL_EMX_MICROGRAPH_DEFOCUSU,defocusU);
            rowIn.getValue(MDL_EMX_MICROGRAPH_DEFOCUSV,defocusV);
            rowIn.getValue(MDL_EMX_MICROGRAPH_ASTIGMATISM_ANGLE,defocusAngle);
            rowIn.getValue(MDL_EMX_MICROGRAPH_VOLTAGE,voltage);
            rowIn.getValue(MDL_EMX_MICROGRAPH_CS,sphericalAberration);
            rowIn.getValue(MDL_EMX_MICROGRAPH_AMPLITUDE_CONTRAST,Q0);

            rowOut.setValue(MDL_IMAGE,micrographXmipp);
            rowOut.setValue(MDL_CTF_SAMPLING_RATE,samplingRate);
            rowOut.setValue(MDL_CTF_DEFOCUSU,defocusU);
            rowOut.setValue(MDL_CTF_DEFOCUSV,defocusV);
            rowOut.setValue(MDL_CTF_DEFOCUS_ANGLE,defocusAngle);
            rowOut.setValue(MDL_CTF_VOLTAGE,voltage);
            rowOut.setValue(MDL_CTF_CS,sphericalAberration);
            rowOut.setValue(MDL_CTF_Q0,Q0);

            mdCTFMicrographXmipp.addRow(rowOut);
        }
        setMetadataVersion("XMIPP_STAR_1");
        mdCTFMicrographXmipp.setComment(comment);
        mdCTFMicrographXmipp.write(fn_out);
    }
    void convertEmx2XmippClass(void)
    {
        MetaData mdClassEMX;
        MetaData mdClassXmipp;
        MetaData mdMicAggregate;
        String particleName;
        String _class;
        MDRow  rowIn, rowOut;

        const MDLabel MyLabels[]       =
            {
                MDL_EMX_PARTICLE_URL,
                MDL_EMX_PARTICLE_CLASS_ID,
            };
        std::vector<MDLabel> myVector(MyLabels,MyLabels+2);
        mdClassEMX.read(tmpname);
        for (std::vector<MDLabel>::iterator it = myVector.begin(); it != myVector.end(); ++it)
        {
            if (!mdClassEMX.containsLabel(*it) )
                REPORT_ERROR(ERR_MD_BADLABEL,
                             (String)"Label: " + MDL::label2Str(*it) + " missing.");
        }
        mdMicAggregate.aggregate(mdClassEMX,AGGR_COUNT,MDL_EMX_PARTICLE_CLASS_ID,
                                 MDL_UNDEFINED,MDL_COUNT);

        std::map< String, int> mapReferences;
        FOR_ALL_OBJECTS_IN_METADATA(mdMicAggregate)
        {
        	mdMicAggregate.getValue(MDL_EMX_PARTICLE_CLASS_ID,_class,__iter.objId);
            mapReferences[_class] = (int)__iter.objId;
        }
        FOR_ALL_OBJECTS_IN_METADATA(mdClassEMX)
        {
        	mdClassEMX.getRow(rowIn, __iter.objId);
            rowIn.getValue(MDL_EMX_PARTICLE_URL,particleName);
            rowIn.getValue(MDL_EMX_PARTICLE_CLASS_ID,_class);

            rowOut.setValue(MDL_IMAGE,particleName);
            rowOut.setValue(MDL_REF,mapReferences[_class]);
            mdClassXmipp.addRow(rowOut);
        }
        setMetadataVersion("XMIPP_STAR_1");
        mdClassXmipp.setComment(comment);
        mdClassXmipp.write(fn_out);
    }

    void run()
    {
        show();

        if(toxmipp)
            switch (conversionType)
            {
            case EMX_CLASS:
                removeDots();
                convertEmx2XmippClass();
                break;
            case EMX_COORDINATES:
                removeDots();
                convertEmx2XmippCoordinates();
                break;
            case EMX_CTFMICROGRAPH:
                removeDots();
                convertEmx2XmippCTFMicrograph();
                break;
            default:
                REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,"Congratulations you have found a bug in convertXmipp2EmxCoordinates");
                break;
            }
        else if (toemx)
            switch (conversionType)
            {
            case EMX_CLASS:
                convertXmipp2EmxClass();
                restoreDots();
                break;
            case EMX_COORDINATES:
                convertXmipp2EmxCoordinates();
                restoreDots();
                break;
            case EMX_CTFMICROGRAPH:
                convertXmipp2EmxCTFMicrograph();
                restoreDots();
                break;
            default:
                REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,"Congratulations you have found a bug in convertXmipp2EmxCoordinates");
                break;
            }
        else
            REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,"Congratulations you have found a bug in convertXmipp2EmxCoordinates");

    }
}
;//end of class ProgMetadataConvert22Emx

RUN_XMIPP_PROGRAM(ProgMetadataConvert22Emx);
