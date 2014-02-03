/***************************************************************************
 * Authors:     joaquin Oton (joton@cnb.csic.es)
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

#include "xmipp_hdf5.h"
#include "xmipp_strings.h"
#include "xmipp_error.h"

struct H5TreeInfo
{
    std::string rootname;
    std::ostream * out;
};

herr_t showObjectInfo(hid_t objId, const char *name, void *op_data)
{
    H5TreeInfo &h5Info = *((H5TreeInfo*)op_data);
    std::ostream &out = *(h5Info.out);

    hsize_t nobj;
    herr_t err=0;
    hid_t grpid, dsid;

    H5G_stat_t statbuf;

    // Print the object name
    out <<  formatString("%s%s - ", h5Info.rootname.c_str(), name);

    H5Gget_objinfo(objId, name, 0, &statbuf);
    /*
     * process each object according to its type
     */
    switch(statbuf.type)
    {
    case H5G_LINK:
        out <<  " SYM_LINK:\n";
        //        do_link(gid,memb_name);
        break;
    case H5G_GROUP:
        {
            String rootname = h5Info.rootname;
            h5Info.rootname += (String)name + "/";

            grpid = H5Gopen(objId,name, H5P_DEFAULT);
            err = H5Gget_num_objs(grpid, &nobj);
            out <<  formatString("Group {%d elements}\n", nobj);

            H5Giterate(objId, name, NULL, showObjectInfo, &h5Info);
            h5Info.rootname = rootname;
            break;
        }
    case H5G_DATASET:
        out <<  "Dataset {";
        dsid = H5Dopen(objId, name, H5P_DEFAULT);
        hsize_t dims[4];
        hid_t filespace;
        int rank;
        filespace = H5Dget_space(dsid);    /* Get filespace handle first. */
        rank  = H5Sget_simple_extent_dims(filespace, dims, NULL);

        for (int k = 0; k < rank-1; ++k)
            out << dims[k] << ", ";

        out << dims[rank-1] << "}\n";
        H5Dclose(dsid);
        break;
    case H5G_TYPE:
            out <<  " Data type:\n";
        //        idType = H5Topen(objId,memb_name, H5P_DEFAULT);
        //        H5Tclose(idType);
        break;
    default:
            out << " unknown?\n";
        break;
    }
    return err;
}



std::map<String, H5infoProvider > createProviderMap()
{
    std::map<String, H5infoProvider > m;
    m["NXtomo"] = std::make_pair(MISTRAL, "/NXtomo/instrument/sample/data");
    m["MDF"]  = std::make_pair(EMAN,    "/MDF/images/%i/image");
    return m;
}


void XmippH5File::openFile(const H5std_string& name, unsigned int flags,
                           const H5::FileAccPropList& access_plist)
{
    if ( isHdf5(name.c_str()) )
        H5::H5File::openFile(name.c_str(), flags, access_plist);
    else
        REPORT_ERROR(ERR_IMG_UNKNOWN, formatString("XmippH5File: Format of %s is not HDF5.",name.c_str()));
}

bool XmippH5File::checkDataset(const char* dsname) const
{
    H5::DataSet dataset;
    // Open the dataset
    try
    {
        dataset = openDataSet(dsname);
    }
    catch (H5::Exception &h5e)
    {
        return false;
    }
    dataset.close();
    return true;
}

int XmippH5File::getDataset(const char* dsname, Matrix1D<double> &data, bool reportError) const
{
    H5::DataSet dataset;
    // Open the dataset
    try
    {
        dataset = openDataSet(dsname);
    }
    catch (H5::Exception &h5e)
    {
        if ( reportError )
            REPORT_ERROR(ERR_ARG_MISSING,formatString("getDataset: %s dataset " \
                         "does not exist in file %s.", dsname, this->getFileName().c_str()));
        //        std::cerr << "getDataset Error: " << h5e.getCDetailMsg() << std::endl;
        return -1;
    }

    //Get dataspace of the dataset.
    H5::DataSpace filespace = dataset.getSpace();

    // Check the number of dimensions in the dataspace is one.
    if (filespace.getSimpleExtentNdims()!= 1 )
        REPORT_ERROR(ERR_ARG_INCORRECT, formatString("getDataset: Dataset %s has "\
                     "more than 1 dimension", dsname));

    // Get the size of the dimension in the dataspace
    hsize_t dim[1];
    filespace.getSimpleExtentDims(dim);


    hsize_t offset[1]; // Hyperslab offset in the file
    //    hsize_t  count[1]; // Size of the hyperslab in the file

    // Define the offset and count of the hyperslab to be read.
    offset[0] = 0;
    //    count[0]  = dim[0];

    filespace.selectHyperslab( H5S_SELECT_SET, dim, offset );

    // Allocate space for data and define the memspace
    data.resizeNoCopy((int)*dim);

    H5::DataSpace memspace(1, dim);

    // Read data from hyperslab in the file into the hyperslab in memory
    dataset.read(MATRIX1D_ARRAY(data), H5::PredType::NATIVE_DOUBLE, memspace, filespace);

    filespace.close();
    memspace.close();
    dataset.close();

    return 0;
}

void XmippH5File::showTree(std::ostream &out)
{
    H5TreeInfo h5Info;

    h5Info.out = &out;
    h5Info.rootname = "";

    this->iterateElems("/", NULL, showObjectInfo, &h5Info);
}


H5infoProvider getProvider(hid_t fhdf5)
{
    H5infoProvider provider;

    size_t maxSize = 1024;
    char groupName[1024];
    char memName[1024];

    hid_t gid;
    ssize_t len;

    gid = H5Gopen(fhdf5,"/", H5P_DEFAULT);

    len = H5Iget_name(gid, groupName, maxSize);

    if (len == 0)
        REPORT_ERROR(ERR_VALUE_EMPTY, "rwHDF5: Empty structure in file.");

    len = H5Gget_objname_by_idx(gid, 0, memName, maxSize);

    typedef std::map<String, H5infoProvider >::const_iterator it_type;

    for ( it_type it = H5ProviderMap.begin(); it != H5ProviderMap.end(); it++)
    {
        if ( strcmp(memName,it->first.c_str() ) == 0 )
            return it->second;
    }

    return std::make_pair(NONE , "");

//    REPORT_ERROR(ERR_IO, "rwHDF5: Unknown file provider. Default dataset unknown.");

}
