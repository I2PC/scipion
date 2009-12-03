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

void Usage(char **argv);
void dm2spi(const FileName &, const FileName &, bool);


struct image_data{
	 int Xdim, Ydim, header, dataType;
 };


 int readTagDM3(FILE *fh_in, bool bigEndian, int depLevel, int index[], image_data* &imData, int &imCount);
 void FREADTagValue(void *fieldValue, int numberType, int n, FILE* fh_in, bool bigEndian);
 char* sprintfTagValue(void *value, int numberType);


int main(int argc, char *argv[])
{
    FileName       fn_in;    // input file
    FileName       fn_out;   // output file
    FileName       fn_sel;   // input selfile
    FileName       fn_oext;  // output extension
    
    bool           reverse_endian;

    if( IsBigEndian())
    	EXIT_ERROR(1, "md2spi: This program only works for little endian boxes");

    /* Parameters ============================================================== */
    try
    {
        if (argc == 1) Usage(argv);
        if (checkParameter(argc, argv, "-i"))
        {
            fn_in = getParameter(argc, argv, "-i");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "md2spi: -i option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "-o"))
        {
            fn_out = getParameter(argc, argv, "-o");
            if (checkParameter(argc, argv, "-sel") || checkParameter(argc, argv, "-oext"))
                EXIT_ERROR(1, "md2spi: -o option is not compatible with -sel or -oext");
        }
        if (checkParameter(argc, argv, "-sel"))
        {
            fn_sel = getParameter(argc, argv, "-sel");
            fn_oext  = getParameter(argc, argv, "-oext", "xmp");
            if (checkParameter(argc, argv, "-i") || checkParameter(argc, argv, "-o"))
                EXIT_ERROR(1, "md2spi: -sel option is not compatible with -i or -o");
        }

        reverse_endian = checkParameter(argc, argv, "-reverse_endian");
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
        if (fn_sel!="")
        {
            SelFile SF(fn_sel), SF_out;
            std::cerr << "Converting from MD to SPI ...\n";
            init_progress_bar(SF.ImgNo());
            int i=0;
            while (!SF.eof())
            {
                FileName in_name = SF.NextImg();
                FileName out_name = in_name.without_extension()+"."+fn_oext;
                SF_out.insert(out_name);
                dm2spi(in_name, out_name, reverse_endian);
                progress_bar(i++);
            }
            progress_bar(SF.ImgNo());
            SF_out.write(fn_sel.without_extension()+"_spider.sel");
        }

        /* input/output are single files */
        else if (fn_in!="" && fn_out!="")
            dm2spi(fn_in, fn_out, reverse_endian);
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
    printf(
        "Usage: %s [Purpose and Parameters]"
        "\nPurpose: Convert from MD to Spider format "
        "\nParameter Values: (note space before value)"
        "\nESPECIFIC PARAMETERS FOR SINGLE-FILE CONVERSION"
        "\n    -i    file_in        input md file"
        "\n    -o    file_out       output Spider file"
        "\nESPECIFIC PARAMETERS FOR SEL-FILE CONVERSION"
        "\n    -sel  input_file     input sel file"
        "\n    -oext input_file     extension for the output files"
        "\nGENERAL PARAMETERS"
        "\n   [-reverse_endian]     by default, output has the same endiannes as input"
        "\n"
        , argv[0]);
}

void dm2spi(const FileName &fn_in,
		const FileName &fn_out,
		bool reverse_endian)
{

	unsigned int  fileVersion, dummy, byteOrder, nrTags;
	FILE *fh_in, *fh_out;


	if ( ( fh_in = fopen(fn_in.c_str(), "r") ) == NULL )
		REPORT_ERROR(6001, "Cannot open file (md2spi).");

	/* See
	/* Header ============================================================== */

	FREAD(&fileVersion,sizeof(int),1,fh_in,1);
	FREAD(&dummy,sizeof(int) ,1,fh_in,1); //file length - 16 = size of root tag directory
	FREAD(&byteOrder,sizeof(int)  ,1,fh_in,1); //byte order, 0 = big endian (Mac) order,1 = little endian (PC) order
	bool bigEndian;
	if(byteOrder)
		bigEndian = false;
	else
		bigEndian = true;

	/* Root tag directory ===================================================== */

	FREAD(&dummy,1,2,fh_in,0); // skip sorted and open
	FREAD(&nrTags,sizeof(int),1,fh_in,1); //  number of tags in root directory (12h = 18)

	printf( "fileVersion: %d\n", fileVersion);
	printf( "littleEndian: %d\n", byteOrder);
	printf( "Number of tags: %d\n\n", nrTags);


	/* Read of the root tags ===================================================== */

	int depLevel=0, imCount=0, imCountF=0;
	image_data* imData;

	for (int n=1;n<=nrTags;n++)
		readTagDM3(fh_in,bigEndian, depLevel, &n, imData, imCount);


	/* Select the located images =================================================== */

//	if (imCount>2)
//		REPORT_ERROR(6001, "More than one image in file, this option is not implemented.");

	int imIndex[imCount];

	for (int n=0;n<imCount;n++) {			// Select all images except thumbnails

//		printf("%d		%d		%d		%d \n",(imData+n)->dataType,(imData+n)->Xdim,(imData+n)->Ydim,(imData+n)->header );
		if ((imData+n)->dataType==7){ // (thumbnail=23 / image=7)
			imIndex[imCountF] = n;
			imCountF++;
		}
	}

	/* Save the located images =================================================== */

	FileName fn_outF;

#define BUFFSIZE 1024*1024

	for (int n=0;n<imCountF;n++)
	{
		if (imCountF==1)
			fn_outF = fn_out;
		else
			fn_outF.compose(fn_out.without_extension(),n+1,fn_out.get_extension());

		/* Write image to file ==================================*/

		if ( ( fh_out = fopen(fn_outF.c_str(), "wb") ) == NULL )
					REPORT_ERROR(6001, "Cannot create out file (md2spi).");

		int * imBuffer;
		float * imFloatBuffer;
		imBuffer = (int *) malloc (BUFFSIZE * sizeof(int));
		imFloatBuffer = (float *) malloc (BUFFSIZE * sizeof(float));

		fseek(fh_in,(imData+imIndex[n])->header,SEEK_SET);


		unsigned int bytesLeft;
		unsigned int end_header=(sizeof(int)*(imData+imIndex[n])->Xdim*(imData+imIndex[n])->Ydim)+(imData+imIndex[n])->header;
		bytesLeft = end_header - ftell(fh_in);

//		init_progress_bar(bytesLeft);

		while (bytesLeft > BUFFSIZE*4)
		{
			FREAD(imBuffer, sizeof(int),BUFFSIZE, fh_in, bigEndian);
			for( int l=0 ; l < BUFFSIZE; l++)
				imFloatBuffer[l]=(float) imBuffer[l];
			FWRITE((void *) imFloatBuffer, sizeof(int), BUFFSIZE, fh_out, bigEndian);
			bytesLeft = end_header - ftell(fh_in);

//			progress_bar(ftell(fh_in)-(imData+imIndex[n])->header);
		}
		if (bytesLeft > 0)
		{
			FREAD(imBuffer, sizeof(int),bytesLeft/sizeof(int), fh_in, reverse_endian);
			for( int l=0 ; l < bytesLeft/sizeof(int); l++)
				imFloatBuffer[l]=(float) imBuffer[l];
			FWRITE((void *) imFloatBuffer, sizeof(int), bytesLeft/sizeof(int), fh_out, false);

//			progress_bar(ftell(fh_in)-(imData+imIndex[n])->header);

			if (fclose(fh_out)!=0)
				REPORT_ERROR(6001, "Error creating OUT file (md2spi).");


			/* Write INF file ==================================*/

			fn_outF=fn_outF.add_extension("inf");
			if ( ( fh_out = fopen(fn_outF.c_str(), "w") ) == NULL )
								REPORT_ERROR(6001, "Cannot create INF file (md2spi).");

			fprintf(fh_out,"# Bits per sample\n");
			fprintf(fh_out,"bitspersample= %d\n",sizeof(int)*8);
			fprintf(fh_out,"# Samples per pixel\n");
			fprintf(fh_out,"samplesperpixel= 1\n");
			fprintf(fh_out,"# Image width\n");
			fprintf(fh_out,"Xdim= %d\n",(imData+imIndex[n])->Xdim);
			fprintf(fh_out,"# Image length\n");
			fprintf(fh_out,"Ydim= %d\n",(imData+imIndex[n])->Ydim);
			fprintf(fh_out,"# offset in bytes (zero by default)\n");
			fprintf(fh_out,"offset= 0\n");
			fprintf(fh_out,"# Is a signed or Unsigned int (by default true)\n");
			fprintf(fh_out,"is_signed = true\n");

			if (fclose(fh_out)!=0)
							REPORT_ERROR(6001, "Error creating INF file (md2spi).");

		}
	}
}



int readTagDM3(FILE *fh_in,
		         bool bigEndian,
		         int depLevel,
		         int index[],
		         image_data * & imData, ////////////////////////////////////////////
		         int &imCount)
{
	depLevel++;

	/* Header Tag ============================================================== */
	
	unsigned char cdTag;
	unsigned int  idTag;
	unsigned short int ltName;
	FREAD(&cdTag,sizeof (unsigned char),1,fh_in,false); // Identification tag: 20 = tag dir,  21 = tag
	FREAD(&ltName,sizeof(unsigned short int), 1,fh_in,true); // Length of the tag name
	idTag = int(cdTag);

	char * tagName;
	tagName =  new char[ltName+1];
	FREAD(tagName,ltName,1,fh_in,false); // Tag name
	tagName[ltName] = '\0';

	for (int n=1;n<=depLevel;n++)
		printf("%d.",index[n-1]);


	/* Reading tags ===================================================================*/
	if (idTag == 20)		// Tag directory
	{
		printf("- Dir: %s\n",tagName);
		unsigned char dummy;
		unsigned int nTags;
		FREAD(&dummy,sizeof(unsigned char),1,fh_in,false); // 1 = sorted (normally = 1)
		FREAD(&dummy,sizeof(unsigned char),1,fh_in,false); //  0 = closed, 1 = open (normally = 0)
		FREAD(&nTags,sizeof(int),1,fh_in,true);             //  number of tags in tag directory
	
		int * newIndex;
		newIndex = new int[depLevel+1];
		for (int k=0;k<depLevel;k++)
			newIndex[k]=index[k];

		if (strcmp(tagName,"ImageList")==0)    // Number of images
		{
			imData = new image_data[nTags];
		}
		else if (strcmp(tagName,"Dimensions")==0)
		{
			newIndex[depLevel] = 1;
			(imData+imCount)->Ydim = readTagDM3(fh_in, bigEndian, depLevel,newIndex,imData,imCount);

			newIndex[depLevel] = 2;
			(imData+imCount)->Xdim = readTagDM3(fh_in, bigEndian, depLevel,newIndex,imData,imCount);

			imCount++;

			return 0;
		 }

		
		for (int n=1;n<=nTags;n++)
		{
			newIndex[depLevel] = n;
			readTagDM3(fh_in, bigEndian, depLevel,newIndex,imData,imCount);
		}

		return 0;
	}
	else if (idTag == 21)    // Tag
	{
		printf("- Tag: %s ",tagName);


		unsigned int nnum;
		char	buf[4]; // to read %%%% symbols

		FREAD(&buf,1,4,fh_in,false); // To read %%%% symbols
		FREAD(&nnum,sizeof(unsigned int),1,fh_in,true); // Size of info array

		unsigned int * info;
		info = new unsigned int[nnum];
		FREAD(info,sizeof(unsigned int),nnum,fh_in,true); // Reading of Info


		/* Tag classification  =======================================*/

		if (nnum == 1)   // Single entry tag
		{
			int tagValue=0;

			FREADTagValue(&tagValue,info[0],1,fh_in,bigEndian);
			printf(" = %s\n", sprintfTagValue(&tagValue,info[0]));

			if (strcmp(tagName,"DataType")==0)
			{
				(imData+imCount)->dataType = tagValue;
			}

			return tagValue;
		}
		else if(nnum == 3 && info[0]==20)			// Tag array
		{ 										  /*nnum = 3
													info(0) = 20
													info(1) = number type for all values
													info(2) = info(nnum) = size of array*/

			if (strcmp(tagName,"Data")==0)    // Locating the image data
			{
				(imData+imCount)->header = ftell(fh_in);
			}

				// Jump the array values
				int k;
				if(info[1] == 2 || info[1] == 4) k = 2;
				else if(info[1] == 3 || info[1] == 5) k = 4;
				else if(info[1] == 10 ) k = 1;

	//			printf("Position before array = %d \n ", ftell(fh_in));
//				printf("(Array size = %d) \n", info[nnum-1]*k);

				fseek( fh_in, ftell(fh_in)+(info[nnum-1])*k , SEEK_SET );

	//			printf("Position = %d \n ", ftell(fh_in));
	//		}else{
	//			int arrayValue[info[2]];
	//			memset(arrayValue,0,sizeof(arrayValue));
	//			FREADTagValue(&arrayValue,info[1],info[2],fh_in,bigEndian);
	//			for (int n=0;info[2];n++){
	//				printf("%d  ", arrayValue[n]);
	//			} std::cout << "\n\n";

				return 0;
		
		}
		else if (info[0]==20 && info[1] == 15)    // Tag Group array
		{
													/*nnum = size of array
													info(0) = 20 (array)
													info(1) = 15 (group)
													info(2) = 0 (always 0)
													info(3) = number of values in group
													info(2*i+3) = number type for value i
													info(nnum) = size of info array*/

//			printf(" = (");

			int nBytes=0, k, fieldValue;


			for (int n=1;n<=info[3];n++)
			{
				fieldValue=0;

				FREADTagValue(&fieldValue,info[3+2*n],1,fh_in,bigEndian);

				if(info[3+2*n] == 2 || info[3+2*n] == 4) k = 2;
				else if(info[3+2*n] == 3 || info[3+2*n] == 5) k = 4;
				else if(info[3+2*n] == 10 ) k = 1;
				nBytes+=k;
//				printf( "%s", sprintfTagValue(&fieldValue,info[3+2*n]));
//				if(n<info[3]) printf(","); else printf(")\n");
			}

			// Jump the array values

//			printf("Reading of array values\n");
	//		printf("Position = %d \n ", ftell(fh_in));

			fseek( fh_in, ftell(fh_in)+(info[nnum-1]-1)*nBytes , SEEK_SET );
			return 0;
		}
		else if (info[0] == 15)    // Tag Group  (struct)
		{
									/* info[1] = length group name (always =0)
									   info[2] = number of entries in group */
			/*nnum = size of info array
			info(1) = 0fh
			info(3) = number of values in group
			info(2*i+3) = number type for value i
			Other info entries are always zero*/

			printf(" = [");

			for (int n=1;n<=info[2];n++)
			{
				int fieldValue=0;
				FREADTagValue(&fieldValue,info[2+2*n],1,fh_in,bigEndian);

				printf( "%s", sprintfTagValue(&fieldValue,info[2+2*n]));
				if(n<info[2]) printf(","); else printf("]\n");
			}
//					printf("Position = %d \n ", ftell(fh_in));

			return 0;
		}
	}
}


void FREADTagValue(void *fieldValue,
					  int numberType,
					  int n, FILE* fh_in,
					  bool bigEndian)
{

	switch(numberType)
	{
		case 2:						// (02h =  2  i2* signed    (short)
			short* shortValue;
			shortValue = (short*)fieldValue;
			FREAD(shortValue,sizeof(short),n,fh_in,bigEndian);
			break;
		case 3:								  // 03h =  3  i4* signed    (long)
	//		int* intValue[n];
	//		*intValue = (int*) malloc(sizeof(intValue));
	//		FREAD(intValue,sizeof(int),n,fh_in,bigEndian);
	//		fieldValue = (int*)intValue;
			FREAD(fieldValue,sizeof(int),n,fh_in,bigEndian);
			break;
		case 4:							//  04h =  4  i2* unsigned  (ushort) or unicode string
			unsigned short* uShortValue;
			uShortValue = (unsigned short*)fieldValue;
			FREAD(shortValue,sizeof(unsigned short),n,fh_in,bigEndian);
			break;
		case 5:								//  05h =  5  i4* unsigned  (ulong)
			unsigned int* uIntValue;
			uIntValue = (unsigned int*)fieldValue;
			FREAD(uIntValue,sizeof(unsigned int),n,fh_in,bigEndian);
			break;
		case 6:								//  06h =  6  f4*           (float)
			float* floatValue;
			floatValue = (float *) fieldValue;
			FREAD(floatValue,sizeof(float),n,fh_in,bigEndian);
			break;
		case 7:								//  07h =  7  f8*           (double)
			double* doubValue;
			doubValue = (double*) fieldValue;
			FREAD(doubValue,sizeof(double),n,fh_in,bigEndian);
			break;
		case 8:								//  08h =  8  i1            (boolean)
			bool* boolValue;
			boolValue = (bool*) fieldValue;
			FREAD(boolValue,sizeof(bool),n,fh_in,bigEndian);
			break;
		case 10:								//  0ah = 10  i1
			char* cValue;
			cValue = (char*) fieldValue;
			FREAD(cValue,sizeof(char),n,fh_in,bigEndian);
			break;
	}

}




char* sprintfTagValue(void *fieldValue, int numberType)
{

	char str[10];
	int len;

	if (numberType==2){
		short* tempValue; tempValue = (short*)fieldValue;
		len = sprintf(str,"%d",*tempValue);

	} else if (numberType==3){
		int* tempValue; tempValue = (int*)fieldValue;
		len = sprintf(str,"%d",*tempValue);
		
	} else if (numberType==4){
		short* tempValue; tempValue = (short*) fieldValue;
		len = sprintf(str,"%d",*tempValue);
	
	} else if (numberType==5){
		unsigned int* tempValue; tempValue = (unsigned int*) fieldValue;
		len = sprintf(str,"%u",*tempValue);

	} else if (numberType==6){
		float* tempValue; tempValue = (float *) fieldValue;
		len = sprintf(str,"%4.3e",*tempValue);

	} else if (numberType==7){
		double* tempValue;tempValue = (double*) fieldValue;
		len = sprintf(str,"%4.3e",*tempValue);

	} else if (numberType==8){
		bool* tempValue;tempValue = (bool*) fieldValue;
		len = sprintf(str,"%d",*tempValue);

	} else if (numberType==10){
		char* tempValue; tempValue = (char*) fieldValue; ;
		len = sprintf(str,"%s",*tempValue);
	}	


char out[len+1];

for (int k=0;k<len;k++) out[k]=str[k];
out[len]='\0';

return out;

}



 /********************************
  *
  * Digital Micrograph file format

Digital Micrograph is an image processing program produced commercially by Gatan.

Gatan does not publish the file format for Digital Micrograph. This information has been obtained by examining the structure of files, thus it may be inaccurate or wrong and is definitely incomplete.

See also the Greg Jefferis' Digital Micrograph 3 file format page on which some of the information here is based.

DM3 info updated March 2006
Digital Micrograph 2 file format

The files examined were written by Digital Micrograph 2.1.5. For DM 2.5 the main difference is the version tag

A ? means a guess or I've not bothered to check or I'm not sure.

Mac file type: GSHN
Mac creator: GCCD
Resource fork:
Seems to contain ~286 bytes, but no resources as seen by ResEdit.
Data fork:
File is arranged in a series of fields, each field containing tag (2 bytes), data length (4 bytes), data (data length bytes).

The 2 byte tag identifies the type of data contained in the field. These are, approximately in the order they appear in the file (is the order significant?)

tag value
(hex)

  3d  DM version multiplied by 100 and stored as i4 (ie length 4). For DM 2.1.5
      the version is 200, for DM 2.5 it is 250.

ffff  The image itself. First 8 bytes specify size and type as 4i2.
      This is the same as the "small header format", except that only data
      types 1 to 7 are listed in the manual for this format.

      1  width
      2  height
      3  bytes/pixel (eg float=4, complex=8)
      4  data type.  1  2 byte integer signed ("short")
                     2  4 byte real (IEEE 754)
                     3  8 byte complex (real, imaginary)
                     4  ?
                     5  4 byte packed complex (see below)
                     6  1 byte integer unsigned ("byte")
                     7  4 byte integer signed ("long")
                     8  rgb view, 4 bytes/pixel, unused, red, green, blue?
                     9  1 byte integer signed
                    10  2 byte integer unsigned
                    11  4 byte integer unsigned
                    12  8 byte real
                    13 16 byte complex
                    14  1 byte binary (ie 0 or 1)

      The first 3 multiplied together should give the total number of bytes in
      the picture.
      The rest (ie all but first 8 bytes) is the image.

      Packed complex (data type 5)
      This is used for the Fourier transform of real images, which have
      symmetric real parts and antisymmetric imaginary parts and thus can
      be stored in half the number of bytes that the equivalent complex
      picture would take. The format is somewhat strange.
      I have confused things further by using semper's coordinate system.
      If the equivalent full complex picture of size n by n would look like
      x1 = -n/2,          x2 = int((n-1)/2)
      y1 = -int((n-1)/2), y2 = n/2

      real part

            x1         ...         -1       0       1      ...      x2
      y1   rx1,y1                          r0,y1                   rx2,y1
      ...
      -1                           r-1,-1  r0,-1   r1,-1
       0   rx1,0                   r-1,0   r0,0    r1,-1           rx2,0
       1                           r-1,1   r0,0    r1,-1
      ...
      y2   rx1,y2                          r0,y2                   rx2,y2

      imaginary part likewise but i-1,-1 etc

      packed complex

            x1     x1+1  ... -2      -1       0       1    ...  x2-1    x2
      y1   rx1,0  *rx1,y1    r1,y1   i1,y1   r2,y1   i2,y1     rx2,y1  ix2,y1
      ...
       1   rx1,y2  ix1,y2    r1,1    i1,1    r2,1    i2,1      rx2,1   ix2,1
       0   r0,0   *r0,y1     r1,0    i1,0    r2,0    i2,0      rx2,0   ix2,0
      -1   r0,-1   i0,-1     r1,-1   i1,-1   r2,-1   i2,-1     rx2,-1  ix2,-1
      ...
      y2   r0,y2   i0,y2     r1,y2   i1,y2   r2,y2   i2,y2     rx2,y2  ix2,y2

      The top of the x1 and x1+1 columns contain what would be in the bottom
      of the x1 column, with two imaginary parts containing real parts (marked
      with *)

  3b  Contains the local info saved with the picture eg mictroscope cs.
      First 4 bytes - number of tags (i4)
      Each tag has the format
      4i2, string, 8i2, string, 10i2
      The integer before each string is the string length
      The integer before this is the string length + 2
      Integer 8 seems to be the type of the tag
        2  string
        3  number
        4  keyword
        5  unknown
      All the rest of the integers were the same in all tags examined.

  3c  Contents of notes box. First 4 bytes are number of characters. Rest is
      text of notes box. There is no trailing null.

  2d  Display type = raster image if present? Length=0
      Also has 16 and 3e set

  2e  Display type = surface plot if present? Length=0
      Also has 2f, 30, 31, 32, 33, 34 set

 1f4  Display type = line plot if present? Length=0
      Also has 1f5, 1f6, 1f7, 1f8, 1f9 and others set

  16  Display magnification (screen pixels/pixel) (real)

  3e  Position of top left of picture with respect to top left of window
      (2i2)

  1b  Picture maximum value (real)

  1c  Picture minimum value (real)

  35  Units for pixel size (null terminated (eg 1/um for fft) plus other
      stuff to total of 16 bytes (or is everything after the null junk?)

  1f  Pixel size in um? (real)

  20  Pixel size in um? (real)

  23  0 = normal contrast, 1 = inverted contrast (i1) set in display info

   d  Colour mode (i2) set in display info
         1  Gray-scale
         2  Linear
         3  Rainbow?
         4  Temperature?
         5  Custom?

   c  Contrast mode (i2) set in display info
         1  Linear
         2  Equalized
         3  Pseudo-contour?
         4  Custom?

  27  0 = survey off, 1 = survey on (i1) set in display info

  28  0 = survey cross-wires, 1 = survey entire image (i2) set in display info

  11  Value to display as black (contrast limits) (real)

  12  Value to display as white (contrast limits) (real)

  26  Minimum contrast (real) set in display info

  25  Annotation, eg text or lines on screen. First 4 bytes is probably number
      of annotations (i4)

  19  position & size of window on screen, top left = 0,0. (4i2)
      top, left, bot, right

   0  End of file (length 0)

Digital Micrograph 3 file format

Mac file type: GTGI
Mac creator: GCCD

Files examined were written from DM 3.3.1 on a PC and a Mac and later versions.
Notation used

The notation is loosely based on Fortran.

  i1   char    1 byte integer
  i2   short   2 byte integer
  i4   long    4 byte integer

  f4   float   4 byte floating point
  f8   double  8 byte floating point

   a   char      string

Byte order

  i4be   big endian, Motorola, Mac, eg 00 00 01 02 for 258
  i4le   little endian, Intel, PC,  eg 02 01 00 00 for 258
  i4*    order depends on byte order flag (3rd i4 integer in file)

Hex values are written eg 14h, ie 14h = 20
Overall file structure

File consists of a header, a tag directory and a group of nulls marking the end of the file. The tag directory contains both tags and more tag directories in a hierarchical structure.

The image itself is in a tag directory called "ImageList". More than one image can be stored in Imagelist.

All numbers relating to the header and tag structure are in big endian byte order. Tag values are in the native order of the machine writing the file.

Example, Mac DM3 file

  00 00 00 03  00 22 59 b9  00 00 00 00
  01 00 00 00  00 12 15 00  11 41 70 70  6c 69 63 61  74 69 6f 6e
  ......
  00 00 00 00 00 00 00 00

Header

  00 00 00 03    i4be    DM version = 3
  00 22 59 b9    i4be    file length - 16 = size of root tag directory
  00 00 00 00    i4be    byte order, 0 = big endian (Mac) order,
                                     1 = little endian (PC) order

Root tag directory

  01             i1      1 = sorted (normally = 1)
  00             i1      0 = closed, 1 = open (normally = 0)
  00 00 00 12    i4be    number of tags in root directory (12h = 18)
  ......

The root tag directory contains both tags and more tag directories (see below).
End of file

  00 00 00 00 00 00 00 00 End of file is marked with 8 nulls

Tag structure in root tag directory

Tag directories and tags are identified by their first byte

  14h = 20      tag directory
  15h = 21      tag
  00            end of file

Tag directory

Example

  14   00 12   44 6f 63 75 6d 65 6e 74 4f 62 6a 65 63 74 4c 69 73 74
  00   00   00 00 00 01
  ......

  14             i1      identifies tag directory (14h = 20)
  00 12          i2be    bytes in tag name (ltname), may be 0
  44 6f 63 75  6d 65 6e 74  4f 62 6a 65  63 74 4c 69  73 74
                 a       tag name, length ltname "DocumentObjectList"

  00             i1      1 = sorted? (can be 0 or 1)
  00             i1      0 = closed?, 1 = open (normally = 0)
  00 00 00 01    i4be    number of tags in tag directory (01h = 1). Can be 0

Tags
Overall structure

         i1      	 identifies tag (15h = 21)
         i2be    	 ltname, bytes in tag name, may be 0
         a       	 tag name, length ltname

      	 a4      	 string "%%%%"
      	 i4be    	 nnum, size of info array following (=1)
      	 i4be * nnum   	 info(nnum), array of nnum integers
      	         	 contains number type(s) for tag values
      	 xx* * nnum    	 tag values (byte order set by byte order flag)
                         byte lengths specified in info(nnum)

Single entry tag

Example

  15   00 0e   41 6e 6e 6f 74 61 74 69 6f 6e 54 79 70 65
  25 25 25 25   00 00 00 01   00 00 00 03   00 00 00 14


  15             i1      identifies tag (15h = 21)
  00 0e          i2be    bytes in tag name (ltname), may be 0
  41 6e 6e 6f  74 61 74 69  6f 6e 54 79  70 65
                 a       tag name, length ltname "AnnotationType"

  25 25 25 25    a4      "%%%%"
  00 00 00 01    i4be    nnum, size of info array following (=1)
  00 00 00 03    i4be    info(nnum), array of nnum i4 integers, in this case just 1
                         contains number type (3 = signed i4*)
  00 00 00 14    i4*     tag value, 14h = 20

For single entry tags:

  nnum = 1
  info(1) = number type

Tag containing a group of data (struct)

Example

  15 0006 4f6666736574  25252525 00000007 0000000f 00000000 00000002
		        00000000 00000006 00000000 00000006
		        00000000 00000000


  15             i1      identifies tag (15h = 21)
  00 06          i2be    bytes in tag name (ltname), may be 0
  4f 66 66 73 65 74
                 a       tag name, length ltname "Offset"

  25 25 25 25    a4      "%%%%"
  00 00 00 07    i4be    nnum, size of info array following (=7)
                         info(nnum)
  00 00 00 0f    i4be    info(1) number type (0fh = group of data)
  00 00 00 00    i4be    info(2) length of groupname? (always = 0)
  00 00 00 02    i4be    info(3) number of entries in group (=2)
  00 00 00 00    i4be    info(4) length of fieldname? (always = 0)
  00 00 00 06    i4be    info(5) number type for value 1 (06h = f4)
  00 00 00 00    i4be    info(6) length of fieldname? (always = 0)
  00 00 00 06    i4be    info(7) number type for value 2 (06h = f4)
  00 00 00 00    i4*     value(1)
  00 00 00 00    i4*     value(2)

For group tags

nnum = size of info array
info(1) = 0fh
info(3) = number of values in group
info(2*i+3) = number type for value i
Other info entries are always zero

Tag containing an array

Example, an image tag

15 0004 44617461 25252525 00000003 00000014 00000002 00000024
		 fdff feff ffff 0000 0100 0200 0300 0400 0500
		 fdff feff ffff 0000 0100 0200 0300 0400 0500
		 fdff feff ffff 0000 0100 0200 0300 0400 0500
		 fdff feff ffff 0000 0100 0200 0300 0400 0500


  15             i1      identifies tag (15h = 21)
  00 04          i2be    bytes in tag name (ltname)
  44 61 74 61    a       tag name, length ltname "Data"

  25 25 25 25    a4      "%%%%"
  00 00 00 03    i4be    nnum, size of info array following (=3)
                         info(nnum)
  00 00 00 14    i4be    info(1), number type (14h = array)
  00 00 00 02    i4be    info(2), number type (02h = i2 signed)
  00 00 00 24    i4be    info(3) = info(nnum), size of array (=36)
  fd ff          i2*     value(1)
  fe ff          i2*     value(2)
  ....                   etc to value(36)

For array tags

nnum = 3
info(1) = 14h
info(2) = number type for all values
info(3) = info(nnum) = size of array

Tag containing an array of groups

Example

15 0004 434c5554 25252525 0000000b 00000014 0000000f 00000000 00000003
                 00000000 00000002 00000000 00000002 00000000 00000002
                 00000100
                 0000 0000 0000
                 0101 0101 0101
		 0202 0202 0202
		 0303 0303 0303
                 .....

  15             i1      identifies tag (15h = 21)
  00 04          i2be    bytes in tag name (ltname)
  43 4c 55 54    a       tag name, length ltname "CLUT"

  25 25 25 25    a4      "%%%%"
  00 00 00 0b    i4be    nnum, size of info array following (=11)
                         info(nnum)
  00 00 00 14    i4be    info(1), number type (14h = array)
  00 00 00 0f    i4be    info(2), number type (0fh = group)
  00 00 00 00    i4be    info(3), length of groupname? (always = 0)
  00 00 00 03    i4be    info(4), number of entries in group (=3)
  00 00 00 00    i4be    info(5), length of fieldname? (always = 0)
  00 00 00 02    i4be    info(6), number type for value 1 (02h = i2)
  00 00 00 00    i4be    info(7), length of fieldname? (always = 0)
  00 00 00 02    i4be    info(8), number type for value 2 (02h = i2)
  00 00 00 00    i4be    info(9), length of fieldname? (always = 0)
  00 00 00 02    i4be    info(10), number type for value 3 (02h = i2)
  00 00 01 00    i4be    info(11) = info(nnum), size of array (=256)
  0000 0000 0000 3i2*    3 values for first element of array
  0101 0101 0101 3i2*    3 values for second element of array
  ....

For arrays of groups

nnum = size of array
info(1) = 14h
info(2) = 0fh
info(4) = number of values in group
info(2*i+4) = number type for value i
info(nnum) = size of info array

Number types

  02h =  2  i2* signed    (short)
  03h =  3  i4* signed    (long)
  04h =  4  i2* unsigned  (ushort) or unicode string
  05h =  5  i4* unsigned  (ulong)
  06h =  6  f4*           (float)
  07h =  7  f8*           (double)
  08h =  8  i1            (boolean)
  09h =  9  a1            (char)
  0ah = 10  i1
  0fh = 15  group of data (struct)
	    info(2) = 0
	    info(3) = number in group
	    info(2*n+4) = 0
	    info(2*n+5) data type for each value in group
  12h = 18  a             (string)
  14h = 20  array of numbers or groups
            info(nnum) = number = ngroup
            info(2) is then treated as info(1) above

General

There is no simple way of finding the length of a type 15 tag without completely decoding it and working out the number of bytes in each data type.

The image itself is in a type 15 tag with name "Data" about half way through the tags. It is thus difficult to find the image as the length and number of the preceeding tags can change between images. One possible lazy way is to search for the string "15h 0004h Data%%%%", the image will start 16 bytes beyond this.

There may be more than one image in the file. Each image will have its own Data tag. Images are numbered from 0.

There may be a "thumbnail" image which can be either before or after the main images in the file. The image number of the thumbnail Data tag is given in the tag with name Thumbnails::ImageIndex.

Useful tags in order of appearance:

Description        info in the notes box (not always present)
Data               the image itself
DataType           as in DM2. Note this is different from the number type above.
                   These values are only for the image data and must be
                   consistent with the number type for the Data tag.
                   There are a number of other DataTypes defined that I've
                   never seen in images

   0            null
   1     i2     2 byte integer signed ("short")
   2     f4     4 byte real (IEEE 754)
   3     c8     8 byte complex (real, imaginary)
   4            obsolete
   5     c4     4 byte packed complex (see DM2)
   6    ui1     1 byte integer unsigned ("byte")
   7     i4     4 byte integer signed ("long")
   8  4*ui1     rgb, 4 bytes/pixel, unused, red, green, blue
   9     i1     1 byte integer signed
  10    ui2     2 byte integer unsigned
  11    ui4     4 byte integer unsigned
  12     f8     8 byte real
  13     c16   16 byte complex
  14     i1     1 byte binary (ie 0 or 1)
  23  4*ui1     rgba, 4 bytes/pixel, 0, red, green, blue. Used for thumbnail images


Dimensions         a type 14 tag containing 2 type 15s with no names
                   (irritatingly) which are image width and height
PixelDepth         bytes/pixel

For CCD images (these follow the tags above)

Acquisition Date   image acquisition date and time, unfortunatley both
Acquisition Time   as strings. Worse still, the date string can be
                   in either UK/international or US order depending on the
                   date settings on the mac or PC. It is thus impossible
                   to unambiguously determine the date from the date string.
ImageIndex         Image number of thumbnail image

Unfortunately the tags describing the image are after the image itself, this is particularly annoying for the image dimensions.

*/
