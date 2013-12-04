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
#include "xmipp_hdf5.h"



DataType ImageBase::datatypeH5(hid_t h5datatype)
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

hid_t ImageBase::H5Datatype(DataType datatype)
{
    switch (datatype)
    {
    case DT_Float:
        return H5T_NATIVE_FLOAT;
    case DT_ULong:
        return H5T_NATIVE_ULONG;
    case DT_Long:
        return H5T_NATIVE_LONG;
    case DT_UInt:
        return H5T_NATIVE_UINT;
    case DT_Int:
        return H5T_NATIVE_INT;
    case DT_UShort:
        return H5T_NATIVE_USHORT;
    case DT_Short:
        return H5T_NATIVE_SHORT;
    case DT_UChar:
        return H5T_NATIVE_UCHAR;
    case DT_SChar:
        return H5T_NATIVE_CHAR;
    case DT_Double:
        return H5T_NATIVE_DOUBLE;
    default:
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, formatString("rwHDF5: %s datatype not implemented for " \
                     " HDF5 datatype", datatype2Str(datatype).c_str()));
        break;
    }

}



int ImageBase::readHDF5(size_t select_img)
{
    bool isStack = false;

    H5infoProvider provider = getProvider(fhdf5); // Dataset name

    int errCode = 0;

    hid_t dataset;    /* Dataset and datatype identifiers */
    hid_t filespace;
    hsize_t dims[4]; // We are not going to support more than 4 dimensions, at this moment.
    hsize_t nobjEman;
    hid_t cparms;
    int rank;

    String dsname = filename.getBlockName();

    // Setting default dataset name
    if (dsname.empty())
    {
        dsname = provider.second;

        switch (provider.first)
        {
        case EMAN: // Images in stack are stored in separated groups
            hid_t grpid;
            grpid = H5Gopen(fhdf5,"/MDF/images/", H5P_DEFAULT);
            /*herr_t err = */ H5Gget_num_objs(grpid, &nobjEman);

            dsname = formatString(dsname.c_str(), IMG_INDEX(select_img));
            break;
        default:
        	break;
        }
    }

    dataset = H5Dopen2(fhdf5, dsname.c_str(), H5P_DEFAULT);

    if( dataset < 0)
        REPORT_ERROR(ERR_IO_NOTEXIST, formatString("readHDF5: Dataset '%s' not found",dsname.c_str()));

    cparms = H5Dget_create_plist(dataset); /* Get properties handle first. */

    // Get dataset rank and dimension.
    filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
    //    rank      = H5Sget_simple_extent_ndims(filespace);
    rank  = H5Sget_simple_extent_dims(filespace, dims, NULL);

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

    DataType datatype = datatypeH5(h5datatype);
    MDMainHeader.setValue(MDL_DATATYPE,(int) datatype);

    // Setting isStack depending on provider
    switch (provider.first)
    {
    case MISTRAL: // rank 3 arrays are stacks
        isStack = true;
        break;
        //    case EMAN: // Images in stack are stored in separated groups
    default:
    	break;
    }


    ArrayDim aDim;
    size_t nDimFile;
    aDim.xdim = dims[rank-1];
    aDim.ydim = (rank>1)?dims[rank-2]:1;
    aDim.zdim = (rank>3 || (rank==3 && !isStack))?dims[rank-3]:1;
    if ( provider.first == EMAN )
        nDimFile = nobjEman;
    else
        nDimFile = ( rank<3 || !isStack )?1:dims[0] ;

    if (select_img > nDimFile)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readHDF5 (%s): Image number %lu exceeds stack size %lu", filename.c_str(), select_img, nDimFile));

    aDim.ndim = replaceNsize = (select_img == ALL_IMAGES)? nDimFile :1 ;
    setDimensions(aDim);

    //Read header only
    if(dataMode == HEADER || (dataMode == _HEADER_ALL && aDim.ndim > 1))
        return errCode;


    // EMAN stores each image in a separate dataset
    if ( provider.first == EMAN )
        select_img = 1;

    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : aDim.ndim;



    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);

    if (dataMode < DATA)   // Don't read  data if not necessary but read the header
        return errCode;

    if ( H5Pget_layout(cparms) == H5D_CONTIGUOUS ) //We can read it directly
        readData(fimg, select_img, datatype, 0);
    else // We read it by hyperslabs
    {
        // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
        //if memory already allocated use it (no resize allowed)
        mdaBase->coreAllocateReuse();

        hid_t       memspace;

        hsize_t offset[4]; // Hyperslab offset in the file
        hsize_t  count[4]; // Size of the hyperslab in the file

        // Define the offset and count of the hyperslab to be read.

        switch (rank)
        {
        case 4:
            count[0] = 1;
        case 3:
            //            if (stack)
            count[rank-3] = aDim.zdim;
            offset[rank-2]  = 0;
        case 2:
            count[rank-2]  = aDim.ydim;
            offset[rank-2]  = 0;
            break;
        }
        count[rank-1]  = aDim.xdim;
        offset[rank-1]  = 0;

        aDim.xdim = dims[rank-1];
        aDim.ydim = (rank>1)?dims[rank-2]:1;
        aDim.zdim = (rank == 4)?dims[1]:1;
        // size_t nDimFile = (rank>2)?dims[0]:1 ;

        // Define the memory space to read a hyperslab.
        memspace = H5Screate_simple(rank,count,NULL);

        size_t data = (size_t) this->mdaBase->getArrayPointer();
        size_t pad = aDim.zyxdim*gettypesize(myT());


        for (size_t idx = imgStart, imN = 0; idx < imgEnd; ++idx, ++imN)
        {

            // Set the offset of the hyperslab to be read
            offset[0] = idx;

            if ( H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                                     count, NULL) < 0 )
                REPORT_ERROR(ERR_IO_NOREAD, formatString("readHDF5: Error selecting hyperslab %d from filename %s",
                             imgStart, filename.c_str()));

            //            movePointerTo(ALL_SLICES,imN);
            // Read
            if ( H5Dread(dataset, H5Datatype(myT()), memspace, filespace,
                         H5P_DEFAULT, (void*)(data + pad*imN)) < 0 )
                REPORT_ERROR(ERR_IO_NOREAD,formatString("readHDF5: Error reading hyperslab %d from filename %s",
                                                        imgStart, filename.c_str()));
        }
        H5Sclose(memspace);
    }

    H5Pclose(cparms);
    H5Sclose(filespace);
    H5Dclose(dataset);

    return errCode;
}

int ImageBase::writeHDF5(size_t select_img, bool isStack, int mode, String bitDepth, CastWriteMode castMode)
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "writeHDF5: Not implemented yet.");
}
