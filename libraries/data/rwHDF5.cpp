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

int ImageBase::readHDF5(size_t select_img)
{
    int errCode = 0;

    hid_t dataset, tid, sid, grp;    /* Dataset and datatype identifiers */
    hid_t filespace;
    hsize_t* dims;
    herr_t status, status_n;
    int rank;

    /*
     * Open the file and the dataset.
     */
    //    fhdf5 = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    //    H5Fget_info;

    //    grp  = H5Gopen2(fhdf5, "/exchange", H5P_DEFAULT);

    String dsname = "exchange/data"; // Dataset name

    dataset = H5Dopen2(fhdf5, dsname.c_str(), H5P_DEFAULT);
    if( dataset < 0)
        REPORT_ERROR(ERR_IO_NOTEXIST, formatString("readHDF5: Dataset '%s' not found",dsname.c_str()));

    // Get dataset rank and dimension.
    filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
    rank      = H5Sget_simple_extent_ndims(filespace);
    dims = new hsize_t[3];
    status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);
    offset = H5Dget_offset(dataset);

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
    aDim.ydim = (rank>2)?dims[rank-2]:1;
    aDim.zdim = (rank == 4)?dims[1]:1;
    aDim.ndim = (select_img == ALL_IMAGES)? ( (rank>2)?dims[0]:1 ) :1;
    replaceNsize = aDim.ndim;
    setDimensions(aDim);


    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : aDim.ndim;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);


    //Read header only
    if(dataMode == HEADER || (dataMode == _HEADER_ALL && aDim.ndim > 1))
        return errCode;


    readData(fimg, select_img, datatype, 0);


    return errCode;

    // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
    //if memory already allocated use it (no resize allowed)
    //    mdaBase->coreAllocateReuse();

}

int ImageBase::writeHDF5(size_t select_img, bool isStack, int mode, String bitDepth, CastWriteMode castMode)
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "writeHDF5: Not implemented yet.");
}
