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

public:
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
            toemx=true;
        else
        {
            setMetadataVersion("EMX1.0");
            if (fn_in.isStar1(false))
                toxmipp=true;
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
    /*xmipp coordenate file example
    # XMIPP_STAR_1 *
    #
    data_kk1
    loop_
    _Xcoor
    _Ycoor
    549        136
    549        204
    data_kk2
    loop_
    _Xcoor
    _Ycoor
    266        325
    277        255
    */
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
        mdCoordinateEMX.write(fn_out);
    }

    void run()
    {
        show();

        if(toxmipp)
            int a=0;
        else if (toemx)
            switch (conversionType)
            {
            case EMX_COORDINATES:
                convertXmipp2EmxCoordinates()
                ;
                break;
            default:
                break;
            }

    }
}
;//end of class ProgMetadataConvert22Emx

RUN_XMIPP_PROGRAM(ProgMetadataConvert22Emx);
