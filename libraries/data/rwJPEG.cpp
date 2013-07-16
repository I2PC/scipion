/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include "xmipp_image_base.h"
#include "../../external/jpeg-8c/jpeglib.h"


//#include <jpeglib.h>

int ImageBase::readJPEG(size_t select_img)
{
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);

    jpeg_stdio_src(&cinfo, fimg);
    jpeg_read_header(&cinfo, TRUE);

    MDMainHeader.setValue(MDL_MIN,0.);
    MDMainHeader.setValue(MDL_MAX,0.);
    MDMainHeader.setValue(MDL_AVG,0.);
    MDMainHeader.setValue(MDL_STDDEV,0.);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,1.);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,1.);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Z,1.);
    MDMainHeader.setValue(MDL_DATATYPE,(int) DT_UChar);

    ArrayDim aDim;
    aDim.xdim = cinfo.image_width;
    aDim.ydim = cinfo.image_height;
    aDim.zdim = 1;
    aDim.ndim = 1;

    setDimensions(aDim);

    replaceNsize = aDim.ndim;

    //Read header only
    if (dataMode == HEADER || (dataMode == _HEADER_ALL )) // Stop reading if not necessary
    {
        jpeg_destroy_decompress(&cinfo);
        return 0;
    }

    /* As we cannot mmap a TIFF File, when this option is passed we are going to mmap
     * the multidimarray of Image
     */

    if (mmapOnRead)
    {
        mmapOnRead = false;
        if (aDim.nzyxdim*gettypesize(DT_UChar) > tiff_map_min_size)
            mdaBase->setMmap(true);
    }

    // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
    //if memory already allocated use it (no resize allowed)
    mdaBase->coreAllocateReuse();

    MD.clear();
    MD.resize(aDim.ndim, MDL::emptyHeader);

    /* Start decompression jpeg here */
    jpeg_start_decompress( &cinfo );

    /* allocate memory to hold the uncompressed image */
    char * buffer = new char [aDim.xdim];
    /* now actually read the jpeg into the raw buffer */
    JSAMPROW row_pointer[1];
    row_pointer[0] = new unsigned char [aDim.xdim*cinfo.num_components];
    /* read one scan line at a time */
    while( cinfo.output_scanline < cinfo.image_height )
    {
        jpeg_read_scanlines( &cinfo, row_pointer, 1 );
        for(size_t i=0; i<cinfo.image_width;i++)
            buffer[i] = row_pointer[0][i*cinfo.num_components];
        setPage2T((cinfo.output_scanline - 1)*cinfo.image_width, buffer, DT_UChar, aDim.xdim);
    }
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress( &cinfo );
    jpeg_destroy_decompress( &cinfo );
    free( row_pointer[0] );

    return 0;
}

/* setup the buffer but we did that in the main function */
void init_buffer(jpeg_compress_struct* cinfo)
{}

/* what to do when the buffer is full; this should almost never
 * happen since we allocated our buffer to be big to start with
 */
boolean empty_buffer(jpeg_compress_struct* cinfo)
{
    return TRUE;
}

/* finalize the buffer and do any cleanup stuff */
void term_buffer(jpeg_compress_struct* cinfo)
{}

int ImageBase::writeJPEG(size_t select_img, bool isStack, int mode, String bitDepth, CastWriteMode castMode)
{
    if (isComplexT())
    {
        REPORT_ERROR(ERR_TYPE_INCORRECT,"rwJPEG: Complex images are not supported by JPEG format.");
        return 0;
    }

    ArrayDim aDim;
    mdaBase->getDimensions(aDim);

    // Volumes are not supported
    if (aDim.zdim > 1)
        REPORT_ERROR(ERR_MULTIDIM_DIM, "rwJPEG: volumes are not supported.");

    //Selection of output datatype
    DataType myTypeID = myT();

    castMode = CW_CONVERT;

    if (mmapOnWrite)
    {
        /* As we cannot mmap a JPEG File, when this option is passed we are going to mmap
         * the multidimarray of Image
         */
        mmapOnWrite = false;
        dataMode = DATA;
        MDMainHeader.setValue(MDL_DATATYPE,(int) myTypeID);

        if (aDim.nzyxdim*gettypesize(myTypeID) > tiff_map_min_size)
            mdaBase->setMmap(true);

        // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
        //if memory already allocated use it (no resize allowed)

        mdaBase->coreAllocateReuse();

        return 0;
    }

    jpeg_compress_struct cinfo;
    jpeg_error_mgr jerr;

    //    struct jpeg_destination_mgr dmgr;
    //
    //    /* create our in-memory output buffer to hold the jpeg */
    //    JOCTET * out_buffer   = new JOCTET[aDim.xdim * aDim.ydim];
    //
    //    /* here is the magic */
    //    dmgr.init_destination    = init_buffer;
    //    dmgr.empty_output_buffer = empty_buffer;
    //    dmgr.term_destination    = term_buffer;
    //    dmgr.next_output_byte    = out_buffer;
    //    dmgr.free_in_buffer      = aDim.xdim * aDim.ydim;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    jpeg_stdio_dest(&cinfo, fimg);

    /* Set JPEG image properties */
    cinfo.image_width = aDim.xdim;      /* image width and height, in pixels */
    cinfo.image_height = aDim.ydim;
    cinfo.input_components = 1;     /* # of color components per pixel */
    cinfo.in_color_space = JCS_GRAYSCALE; /* colorspace of input image */

    jpeg_set_defaults(&cinfo);

    jpeg_start_compress(&cinfo, TRUE);

    JSAMPROW row_pointer[1];        /* pointer to a single row */
    row_pointer[0] = new unsigned char [aDim.xdim];
    char * buffer = (char*) row_pointer[0];
    double min0, max0;
    mdaBase->computeDoubleMinMaxRange(min0, max0, 0, aDim.xdim*aDim.ydim);


    while (cinfo.next_scanline < cinfo.image_height)
    {
        getCastConvertPageFromT((cinfo.next_scanline)*cinfo.image_width, buffer,
                                DT_UChar, cinfo.image_width, min0, max0, castMode);

        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

    //    fwrite(out_buffer, cinfo.dest->next_output_byte - out_buffer,1 ,fimg);

    return(0);
}
