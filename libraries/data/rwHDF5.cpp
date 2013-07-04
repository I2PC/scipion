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


DataType ImageBase::datatypeHDF5(hid_t h5datatype)
{
    H5T_sign_t h5sign = H5Tget_sign(h5datatype);

    //    if (h5sign == H5T_SGN_ERROR)
    //        REPORT_ERROR(ERR_IO, "datatypeHDF5: Integer sign error in dataset.");
    bool sign = (h5sign > H5T_SGN_NONE);
    size_t size = H5Tget_size(h5datatype);

    DataType dt;
    switch(H5Tget_class(h5datatype))
    {
    case H5T_FLOAT:
        {
            switch(size)
            {
            case 4:
                dt = DT_Float;
                break;
            case 8:
                dt = DT_Double;
                break;
            default:
                REPORT_ERROR(ERR_IO_SIZE, "datatypeHDF5: bad datatype size");
            }
        }
        break;
    case H5T_INTEGER:
        {
            switch(size)
            {
            case 1:
                dt = (sign)? DT_SChar : DT_UChar;
                break;
            case 2:
                dt = (sign)? DT_Short : DT_UShort;
                break;
            case 4:
                dt = (sign)? DT_Int : DT_UInt;
                break;
            case 8:
                dt = (sign)? DT_Long : DT_ULong;
                break;
            default:
                REPORT_ERROR(ERR_IO_SIZE, "datatypeHDF5: bad datatype size");
            }
        }
        break;
    case H5T_NO_CLASS:
    default:
        dt = DT_Unknown;
        break;
    }
    return dt;
}

std::string ImageBase::getDefaultDataset(hid_t fhdf5)
{
    size_t maxSize = 1024;
    char groupName[maxSize];
    char memName[maxSize];

    hid_t gid;
    ssize_t len;

    gid = H5Gopen(fhdf5,"/", H5P_DEFAULT);

    len = H5Iget_name(gid, groupName, maxSize);

    if (len == 0)
        REPORT_ERROR(ERR_VALUE_EMPTY, "rwHDF5: Empty structure in file.");

    len = H5Gget_objname_by_idx(gid, 0, memName, maxSize);

    if ( strcmp(memName,"NXtomo") == 0 )
        return (std::string) "/NXtomo/instrument/sample/data";
    else if ( strcmp(memName,"MDF") == 0)
        return (std::string) "/MDF/images/0/image";
    else
        REPORT_ERROR(ERR_IO, "rwHDF5: Unknown file provider. Default dataset unknown.");

}

int ImageBase::readHDF5(size_t select_img)
{
    int errCode = 0;

    hid_t dataset;    /* Dataset and datatype identifiers */
    hid_t filespace;
    hsize_t dims[4]; // We are not going to support more than 4 dimensions, at this moment.
    herr_t status, status_n;
    hid_t       cparms;
    int rank;

    String dsname = filename.getBlockName();

    if (dsname.empty())
        dsname = getDefaultDataset(fhdf5); // Dataset name


    dataset = H5Dopen2(fhdf5, dsname.c_str(), H5P_DEFAULT);

    if( dataset < 0)
        REPORT_ERROR(ERR_IO_NOTEXIST, formatString("readHDF5: Dataset '%s' not found",dsname.c_str()));

    cparms = H5Dget_create_plist(dataset); /* Get properties handle first. */

    // Get dataset rank and dimension.
    filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
    rank      = H5Sget_simple_extent_ndims(filespace);
    status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);

    // Offset only set when it is possible to access to data directly
    offset = (H5D_CONTIGUOUS == H5Pget_layout(cparms))? H5Dget_offset(dataset) : 0;


    //    status = H5Dread(dataset, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, bm_out);

    hid_t h5datatype = H5Dget_type(dataset);

    // Reading byte order
    switch(H5Tget_order(h5datatype))
    {
    case H5T_ORDER_ERROR:
        REPORT_ERROR(ERR_IO, "readHDF5: error reading endianess.");
        break;
    case H5T_ORDER_LE:
        swap = IsBigEndian();
        break;
    case H5T_ORDER_BE:
        swap = IsLittleEndian();
        break;
    default:
        REPORT_ERROR(ERR_IO, "readHDF5: unkonwn endianess type, maybe mixed types.");
        break;
    }

    DataType datatype = datatypeHDF5(h5datatype);
    MDMainHeader.setValue(MDL_DATATYPE,(int) datatype);


    bool isStack = ( rank > 2 );


    ArrayDim aDim;
    aDim.xdim = dims[rank-1];
    aDim.ydim = (rank>1)?dims[rank-2]:1;
    aDim.zdim = (rank == 4)?dims[1]:1;
    size_t nDimFile = (rank>2)?dims[0]:1 ;

    if (select_img > nDimFile)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readHDF5 (%s): Image number %lu exceeds stack size %lu", filename.c_str(), select_img, nDimFile));

    aDim.ndim = replaceNsize = (select_img == ALL_IMAGES)? nDimFile :1 ;
    setDimensions(aDim);


    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : aDim.ndim;


    //Read header only
    if(dataMode == HEADER || (dataMode == _HEADER_ALL && aDim.ndim > 1))
        return errCode;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);


    hid_t       memspace;
    hsize_t     chunk_dims[4];
    hsize_t     count[4];
    hsize_t     offset[4];

    int         rank_chunk;


    if (dataMode < DATA)   // Don't read  data if not necessary but read the header
        return errCode;

    if (H5D_CHUNKED == H5Pget_layout(cparms))
    {
        // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
        //if memory already allocated use it (no resize allowed)
        mdaBase->coreAllocateReuse();

        /*
         * Get chunking information: rank and dimensions
         */
        rank_chunk = H5Pget_chunk(cparms, rank, chunk_dims);

        /*
         * Define the memory space to read a chunk.
         */
        memspace = H5Screate_simple(rank_chunk,chunk_dims,NULL);

        /*
         * Define chunk in the file (hyperslab) to read.
         */
        offset[0] = imgStart;
        offset[1] = 0;
        offset[2] = 0;
        count[0]  = 1;
        count[1]  = chunk_dims[1];
        count[2]  = chunk_dims[2];
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                     count, NULL);
        /*
         * Read chunk back and display.
         */

        status = H5Dread(dataset, h5datatype, memspace, filespace,
                         H5P_DEFAULT, this->mdaBase->getArrayPointer());

        /*
         * Close/release resources.
         */
        H5Sclose(memspace);
    }
    else
        readData(fimg, select_img, datatype, 0);


    H5Pclose(cparms);
    H5Sclose(filespace);
    H5Dclose(dataset);


    return errCode;

    // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
    //if memory already allocated use it (no resize allowed)
    //    mdaBase->coreAllocateReuse();

}

int ImageBase::writeHDF5(size_t select_img, bool isStack, int mode, String bitDepth, CastWriteMode castMode)
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "writeHDF5: Not implemented yet.");
}
