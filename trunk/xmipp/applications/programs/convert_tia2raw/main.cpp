/***************************************************************************
 *
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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
#include <data/selfile.h>
#include <data/image.h>


struct DataHeader
{
    short int endianess, SeriesID, SeriesVersion;
    int DATA_TYPE_ID, TagTypeID, TotalNumberElements, NUMBER_IMAGES,
        OFFSET_ARRAY_OFFSET, numberdimensions, * pDATA_OFFSET;
};

struct ImDataHeader
{
    short int DATA_TYPE,                //DataType
          DATA_TYPE_SIZE;
    int        CalibrationElementX,    //CalibrationElementX
      CalibrationElementY,    //CalibrationElementY
      IMAGE_WIDTH,            //ArraySizeX
      IMAGE_HEIGHT;           //ArraySizeY
    double    CalibrationOffsetX,   //CalibrationOffsetX
         PIXEL_WIDTH,            //CalibrationDeltaX
         CalibrationOffsetY,     //CalibrationOffsetY
         PIXEL_HEIGHT;           //CalibrationDeltaY
    std::string     DATA_TYPE_SIZE_STRING;
    bool      isSigned;

};


void Usage(char **argv);
void tia2raw(const FileName &, const FileName &, float);
void openImage(const FileName &, FileName &, int, float);
void printHeader(DataHeader);
void setDataType(ImDataHeader &);



int main(int argc, char *argv[])
{
    FileName fn_in; // input file
    FileName fn_out; // output file
    FileName fn_sel; // input selfile
    FileName fn_oext; // output extension
    float    dStddev; // Standard Desv


    if (IsBigEndian())
        EXIT_ERROR(1, "tia2raw: This program only works for little endian boxes");

    /* Parameters ============================================================== */
    try
    {
        if (argc == 1)
            Usage(argv);
        if (checkParameter(argc, argv, "-i"))
        {
            fn_in = getParameter(argc, argv, "-i");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "tia2raw: -i option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "--ImageIn"))
        {
            fn_in = getParameter(argc, argv, "--ImageIn");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "tia2raw: --ImageIn option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "-o"))
        {
            fn_out = getParameter(argc, argv, "-o");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "tia2raw: -o option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "--OutRootName"))
        {
            fn_out = getParameter(argc, argv, "--OutRootName", "OUT");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "tia2raw: --OutRootName option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "-sel"))
        {
            fn_sel = getParameter(argc, argv, "-sel");
            fn_oext = getParameter(argc, argv, "-oext", "raw");
            if (checkParameter(argc, argv, "-i") || checkParameter(argc, argv, "-o"))
                EXIT_ERROR(1, "tia2raw: -sel option is not compatible with -i or -o");
        }
        if (checkParameter(argc, argv, "-s"))
            dStddev = textToFloat(getParameter(argc, argv, "-s", "5"));
        else
            dStddev = textToFloat(getParameter(argc, argv, "--stddev", "5"));
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage(argv);
    }

    try
    {
        /* Perform conversion ====================================================== */

        /* input is a sel file*/
        if (fn_sel != "")
        {
            SelFile SF(fn_sel), SF_out;
            std::cerr << "Converting from TIA to RAW ...\n";
            init_progress_bar(SF.ImgNo());
            int i = 0;
            while (!SF.eof())
            {
                FileName in_name = SF.NextImg();
                FileName out_name = in_name.without_extension() + "." + fn_oext;
                SF_out.insert(out_name);
                tia2raw(in_name, out_name, dStddev);
                progress_bar(i++);
            }
            progress_bar(SF.ImgNo());
            SF_out.write(fn_sel.without_extension() + "_raw.sel");
        }

        /* input/output are single files */
        else if (fn_in != "")
        {
            if (fn_out == "")
                tia2raw(fn_in, fn_in.without_extension() + ".raw", dStddev);
            else
                tia2raw(fn_in, fn_out, dStddev);
        }

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    exit(0);
}

/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf("Usage: %s [Purpose and Parameters]"
           "\nPurpose: Convert from Tecnai imaging and analysis (TIA) images to raw format"
           "\nParameter Values: (note space before value)"
           "\nESPECIFIC PARAMETERS FOR SINGLE-FILE CONVERSION"
           "\n    -i,--ImageIn file_in        input TIA file"
           "\n    -o,--OutRootName file_out   output RAW file"
           "\nESPECIFIC PARAMETERS FOR SEL-FILE  CONVERSION"
           "\n    -sel  input_file            input sel file"
           "\n    -oext input_file            extension for the output files"
           "\nGENERAL PARAMETERS"
           "\n   [-s, --stddev dStddev]       Cut values above (below) s standard deviation (default=5)"
           "\n", argv[0]);
}

void tia2raw(const FileName &fn_in, const FileName &fn_out, float dStddev)
{

    FILE *fh_in;
    FileName fn_outF;

    /* Read Main Header ============================================================== */

    DataHeader Header;

    if ((fh_in = fopen(fn_in.c_str(), "r")) == NULL)
        REPORT_ERROR(6001, "Cannot open file (tia2raw).");

    FREAD(&Header.endianess, sizeof(short int), 1, fh_in, false);

    if (Header.endianess != 18761)
        REPORT_ERROR(6001, "Wrong endianess, I only work with little endian");

    FREAD(&Header.SeriesID, sizeof(short int), 1, fh_in, false);
    FREAD(&Header.SeriesVersion, sizeof(short int), 1, fh_in, false);
    FREAD(&Header.DATA_TYPE_ID, sizeof(int), 1, fh_in, false);
    FREAD(&Header.TagTypeID, sizeof(int), 1, fh_in, false);
    FREAD(&Header.TotalNumberElements, sizeof(int), 1, fh_in, false);
    FREAD(&Header.NUMBER_IMAGES, sizeof(int), 1, fh_in, false);
    FREAD(&Header.OFFSET_ARRAY_OFFSET, sizeof(int), 1, fh_in, false);
    FREAD(&Header.numberdimensions, sizeof(int), 1, fh_in, false);

    if (Header.DATA_TYPE_ID != 16674)
        REPORT_ERROR(6001, "ERROR: This script only process images in real space");

    fseek(fh_in, Header.OFFSET_ARRAY_OFFSET, SEEK_SET);

    Header.pDATA_OFFSET = (int *) malloc(Header.NUMBER_IMAGES * sizeof(int));

    FREAD(Header.pDATA_OFFSET, sizeof(int), Header.NUMBER_IMAGES, fh_in, false);

    fclose(fh_in);


    /* Print Header ===================================================== */

    // printHeader(Header);


    /* Loop through the images  ===================================================== */


    for (int n=0; n<Header.NUMBER_IMAGES; n++)
    {
        fprintf(stderr,"\nOPEN IMAGE No: %d of %d",n+1, Header.NUMBER_IMAGES);

        /* Composing filename for multiple images file */

        if (Header.NUMBER_IMAGES == 1)
            fn_outF = fn_out;
        else
            fn_outF.compose(fn_out.without_extension(), n + 1, fn_out.get_extension());

        openImage(fn_in, fn_outF,Header.pDATA_OFFSET[n],dStddev);
    }
}

void openImage(const FileName &fn_in, FileName &fn_out, int Position, float dStddev)
{
    FILE *fh_in, *fh_out;

    if ((fh_in = fopen(fn_in.c_str(), "r")) == NULL)
        REPORT_ERROR(6001, "Cannot open file (tia2raw).");
    fseek(fh_in, Position, SEEK_SET);

    /* Reading calibration values ================================================= */

    ImDataHeader DataH;

    FREAD(&DataH.CalibrationOffsetX, sizeof(double), 1, fh_in, false);
    FREAD(&DataH.PIXEL_WIDTH, sizeof(double), 1, fh_in, false);
    FREAD(&DataH.CalibrationElementX, sizeof(int), 1, fh_in, false);
    FREAD(&DataH.CalibrationOffsetY, sizeof(double), 1, fh_in, false);
    FREAD(&DataH.PIXEL_HEIGHT, sizeof(double), 1, fh_in, false);
    FREAD(&DataH.CalibrationElementY, sizeof(int), 1, fh_in, false);
    FREAD(&DataH.DATA_TYPE, sizeof(short int), 1, fh_in, false);
    FREAD(&DataH.IMAGE_WIDTH, sizeof(int), 1, fh_in, false);
    FREAD(&DataH.IMAGE_HEIGHT, sizeof(int), 1, fh_in, false);

    setDataType(DataH);

    /* Save Raw Image ===============================================================*/

#define BUFFSIZE 1024*1024*2

    int size = DataH.IMAGE_HEIGHT * DataH.IMAGE_WIDTH, low, high;
    int valuesLeft = size;
    double avg=0, stddev=0, temp=0;

    short int * imBuffer;

    imBuffer = (short int *) malloc(BUFFSIZE * sizeof(short int));


    while (valuesLeft > BUFFSIZE)
    {
        FREAD(imBuffer, sizeof(short int), BUFFSIZE,fh_in, false);

        for (int n=0; n<BUFFSIZE; n++)
        {
            temp = (double) imBuffer[n];
            avg += temp;
            stddev += temp * temp;
        }
        valuesLeft -= BUFFSIZE;
    }
    if (valuesLeft > 0)
    {
        FREAD(imBuffer, sizeof(short int), valuesLeft,fh_in, false);

        for (int n=0; n<valuesLeft; n++)
        {
            temp = (double) imBuffer[n];
            avg += temp;
            stddev += temp * temp;
        }
    }

    avg /= size;
    stddev = stddev/size - avg * avg;
    stddev *= size/(size -1);
    stddev = sqrt(stddev);

    fprintf(stderr, "\nstd = %5.4f,   avg= %5.4f\n", stddev, avg);

    low  = int(avg - dStddev * stddev);
    high = int(avg + dStddev * stddev);

    fseek(fh_in, Position+50, SEEK_SET);


    if ((fh_out = fopen(fn_out.c_str(), "wb")) == NULL)
        REPORT_ERROR(6001, "Cannot create output file (tia2raw).");

    valuesLeft = size;

    while (valuesLeft > BUFFSIZE)
    {
        FREAD(imBuffer, sizeof(short int), BUFFSIZE,fh_in, false);

        for (int n=0; n<BUFFSIZE; n++)
        {
            if(imBuffer[n]<low)
                imBuffer[n]=low;
            else if (imBuffer[n]>high)
                imBuffer[n]=high;
        }

        FWRITE(imBuffer, sizeof(short int),BUFFSIZE, fh_out, false);

        valuesLeft -= BUFFSIZE;
    }
    if (valuesLeft > 0)
    {
        FREAD(imBuffer, sizeof(short int), valuesLeft,fh_in, false);

        for (int n=0; n<valuesLeft; n++)
        {
            if(imBuffer[n]<low)
                imBuffer[n]=low;
            else if (imBuffer[n]>high)
                imBuffer[n]=high;
        }

        FWRITE(imBuffer, sizeof(short int),valuesLeft, fh_out, false);
    }

    if (fclose(fh_out) != 0)
        REPORT_ERROR(6001, "Error creating output file (tia2raw).");


    /* Write INF file ============================================================== */

    fn_out = fn_out.add_extension("inf");
    if ((fh_out = fopen(fn_out.c_str(), "w")) == NULL)
        REPORT_ERROR(6001, "Cannot create output info file (tia2raw).");

    fprintf(fh_out, "# Bits per sample\n");
    fprintf(fh_out, "bitspersample= %d\n", DataH.DATA_TYPE_SIZE * 8);
    fprintf(fh_out, "# Samples per pixel\n");
    fprintf(fh_out, "samplesperpixel= 1\n");
    fprintf(fh_out, "# Image width\n");
    fprintf(fh_out, "Xdim= %d\n", DataH.IMAGE_WIDTH);
    fprintf(fh_out, "# Image length\n");
    fprintf(fh_out, "Ydim= %d\n", DataH.IMAGE_HEIGHT);
    fprintf(fh_out, "# offset in bytes (Optional, zero by default)\n");
    fprintf(fh_out, "offset= 0\n");
    fprintf(fh_out, "# Is a signed or Unsigned int (Optional, by default true)\n");
    if (DataH.isSigned)
        fprintf(fh_out, "is_signed = true\n");
    else
        fprintf(fh_out, "is_signed = false\n");
    fprintf(fh_out, "CalibrationOffsetX %5.3e\n", DataH.CalibrationOffsetX);
    fprintf(fh_out, "PIXEL_WIDTH %5.3e\n",DataH.PIXEL_WIDTH);
    fprintf(fh_out, "CalibrationElementX %d\n", DataH.CalibrationElementX);
    fprintf(fh_out, "CalibrationOffsetY %5.3e\n", DataH.CalibrationOffsetY);
    fprintf(fh_out, "PIXEL_HEIGHT %5.3e\n",DataH.PIXEL_HEIGHT);
    fprintf(fh_out, "CalibrationElementY %d\n", DataH.CalibrationElementX);
    fprintf(fh_out, "DATA_TYPE %s\n",DataH.DATA_TYPE_SIZE_STRING.c_str());
    fprintf(fh_out, "IMAGE_WIDTH %d\n",DataH.IMAGE_WIDTH);
    fprintf(fh_out, "IMAGE_HEIGHT %d\n",DataH.IMAGE_HEIGHT);
    fprintf(fh_out, "IMAGE_DATA_TYPE_SIZE %d\n",DataH.DATA_TYPE_SIZE);

    if (fclose(fh_out) != 0)
        REPORT_ERROR(6001, "Error creating output info file (tia2raw).");
}


void printHeader(DataHeader Header)
{
    fprintf(stderr,"\n SeriesID: %d", Header.SeriesID);
    fprintf(stderr,"\n SeriesVersion: %d", Header.SeriesVersion);
    fprintf(stderr,"\n DATA_TYPE_ID: %d", Header.DATA_TYPE_ID);
    fprintf(stderr,"\n TagTypeID: %d", Header.TagTypeID);
    fprintf(stderr,"\n TotalNumberElements: %d", Header.TotalNumberElements);
    fprintf(stderr,"\n NUMBER_IMAGES: %d", Header.NUMBER_IMAGES);
    fprintf(stderr,"\n OFFSET_ARRAY_OFFSET: %d ", Header.OFFSET_ARRAY_OFFSET);
    fprintf(stderr,"\n numberdimensions: %d ", Header.numberdimensions);

    for (int n=0; n<Header.NUMBER_IMAGES; n++)
        fprintf(stderr,"\n DATA_OFFSET [%d]= %d",n, Header.pDATA_OFFSET[n]);


}


void setDataType(ImDataHeader & DataH)
{
    std::string auxString;
    short int dataType=0;

    DataH.isSigned = 0;

    if(DataH.DATA_TYPE==1)
    {
        auxString = "GRAY8";
        dataType =8;
    }
    else if (DataH.DATA_TYPE==2)
    {
        auxString = "GRAY16_UNSIGNED";
        dataType =16;
    }
    else if(DataH.DATA_TYPE==3)
    {
        auxString = "GRAY32_UNSIGNED";
        dataType =32;
    }
    else if(DataH.DATA_TYPE==4)
    {
        auxString = "GRAY8";
        dataType =8;
    }
    else if(DataH.DATA_TYPE==5)
    {
        auxString = "GRAY16_SIGNED";
        dataType =16;
        DataH.isSigned=1;
    }
    else if(DataH.DATA_TYPE==6)
    {
        auxString = "GRAY32_INT";
        dataType =32;
    }
    else if(DataH.DATA_TYPE==7)
    {
        auxString = "GRAY32_FLOAT";
        dataType =32;
    }
    else if(DataH.DATA_TYPE==8)
    {
        auxString = "GRAY64_DOUBLE";
        dataType =64;
    }
    DataH.DATA_TYPE_SIZE = dataType/8;
    DataH.DATA_TYPE_SIZE_STRING = auxString;
}

