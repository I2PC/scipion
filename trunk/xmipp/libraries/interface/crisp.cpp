/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "crisp.h"

void CrispVolume::read(const FileName &fn)
{
    FILE *fp;

    // Clear Image first
    V().clear();
    name = fn;

    // Open file
    if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "CrispVolume::read: File " + fn + " not found");

    // Read header
    char signature[2];
    fread(signature, sizeof(char), 2, fp);
    if (strcmp(signature, "SH") != 0)
        REPORT_ERROR(1502, "CrispVolume::read: File " + fn + " is not a valid Crisp file");

    FREAD(&flags, sizeof(short),  1, fp, true);
    FREAD(&Adata, sizeof(short),  1, fp, true);
    FREAD(&Bdata, sizeof(short),  1, fp, true);
    FREAD(&Cdata, sizeof(short),  1, fp, true);
    FREAD(&Afile, sizeof(short),  1, fp, true);
    FREAD(&Bfile, sizeof(short),  1, fp, true);
    FREAD(&Cfile, sizeof(short),  1, fp, true);
    FREAD(&Asize, sizeof(short),  1, fp, true);
    FREAD(&Bsize, sizeof(short),  1, fp, true);
    FREAD(&Csize, sizeof(short),  1, fp, true);
    FREAD(&Gamma, sizeof(short),  1, fp, true);
    FREAD(dummy0, sizeof(short), 20, fp, true);
    FREAD(Name,  sizeof(char) , 32, fp);
    FREAD(dummy1, sizeof(short), 208, fp, true);

    // Set correct image size readed from header
    V().resize(Cfile, Bfile, Afile);

    // Read Volume voxels
    for (unsigned z = 0; z < Cfile; z++)
        for (unsigned y = 0; y < Bfile; y++)
            for (unsigned x = 0; x < Afile; x++)
            {
                short tmp;
                FREAD(&tmp, sizeof(short), 1, fp, true);
                DIRECT_VOLVOXEL(V, z, y, x) = tmp;
            }

    // Close file
    fclose(fp);
}

std::ostream & operator << (std::ostream &out, const CrispVolume &cv)
{
    out << "Crisp volume =============================\n";
    out << "flags= " << cv.flags << std::endl
        << "Adata= " << cv.Adata << std::endl
        << "Bdata= " << cv.Bdata << std::endl
        << "Cdata= " << cv.Cdata << std::endl
        << "Afile= " << cv.Afile << std::endl
        << "Bfile= " << cv.Bfile << std::endl
        << "Cfile= " << cv.Cfile << std::endl
        << "Asize= " << cv.Asize << std::endl
        << "Bsize= " << cv.Bsize << std::endl
        << "Csize= " << cv.Csize << std::endl
        << "Gamma= " << cv.Gamma << std::endl
        << "Name = " << cv.Name  << std::endl;
    return out;
}

void CrispVolume::write_as_spider(const FileName &fn)
{
    VolumeXmipp aux;
    aux = V;
    aux.write(fn);
}

