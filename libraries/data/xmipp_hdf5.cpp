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
    herr_t err;
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

void XmippH5File::showTree(std::ostream &out)
{
    H5TreeInfo h5Info;

    h5Info.out = &out;
    h5Info.rootname = "";

    this->iterateElems("/", NULL, showObjectInfo, &h5Info);
}

