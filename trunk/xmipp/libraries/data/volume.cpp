/***************************************************************************
 *
 * Authors: Pedro Antonio de Alarc�n (pedro@cnb.uam.es)
 *          Carlos Oscar S. Sorzano
 *          Alberto Pascual Montano
 *          Roberto Marabini
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "volume.h"

/* Write as CCP4 ----------------------------------------------------------- */
void write_as_CCP4(VolumeT<double> *V, const FileName &fn, double Ts)
{
    FILE *fp;
    if (fn == "")
        REPORT_ERROR(1, "write_as_CCP4:: Must give a filename");

    if ((fp = fopen(fn.c_str(), "wb")) == NULL)
        REPORT_ERROR(1503, "write_as_CCP4: File " + fn + " cannot be saved");

    float faux;
    int iaux;
#define WRITE_FAUX fwrite(&faux,sizeof(float),1,fp)
#define WRITE_IAUX fwrite(&iaux,sizeof(int),1,fp)

    iaux =       XSIZE((*V)());
    WRITE_IAUX; //  1: NC
    iaux =       YSIZE((*V)());
    WRITE_IAUX; //  2: NR
    iaux =       ZSIZE((*V)());
    WRITE_IAUX; //  3: NS
    iaux =       2;
    WRITE_IAUX; //  4: Mode=Real Image
    iaux =       STARTINGX((*V)());
    WRITE_IAUX; //  5: NCSTART
    iaux =       STARTINGY((*V)());
    WRITE_IAUX; //  6: NRSTART
    iaux =       STARTINGZ((*V)());
    WRITE_IAUX; //  7: NSSTART
    iaux =       XSIZE((*V)());
    WRITE_IAUX; //  8: NX
    iaux =       YSIZE((*V)());
    WRITE_IAUX; //  9: NY
    iaux =       ZSIZE((*V)());
    WRITE_IAUX; // 10: NZ
    faux = (float)XSIZE((*V)()) * Ts;
    WRITE_FAUX; // 11: X length
    faux = (float)YSIZE((*V)()) * Ts;
    WRITE_FAUX; // 12: Y length
    faux = (float)ZSIZE((*V)()) * Ts;
    WRITE_FAUX; // 13: Z length
    faux = (float)90.0;
    WRITE_FAUX; // 14: Alpha, Cell angle
    WRITE_FAUX; // 15: Beta,  Cell angle
    WRITE_FAUX; // 16: Gamma, Cell angle
    faux = (float)1.0;
    WRITE_FAUX; // 17: MAPC
    faux = (float)2.0;
    WRITE_FAUX; // 18: MAPR
    faux = (float)3.0;
    WRITE_FAUX; // 19: MAPS

    double minval, maxval, avgval, stddev;
    (*V)().compute_stats(minval, maxval, avgval, stddev);
    faux = (float)minval;
    WRITE_FAUX; // 20: Minimum density value
    faux = (float)maxval;
    WRITE_FAUX; // 21: Maximum density value
    faux = (float)avgval;
    WRITE_FAUX; // 22: Mean density
    iaux =       1;
    WRITE_IAUX; // 23: Space group number
    iaux =       0;
    WRITE_IAUX; // 24: Bytes to store symmetry
    iaux =       0;
    WRITE_IAUX; // 25: Skew flag (no skew)
    for (int i = 26; i <= 34; i++)
        WRITE_IAUX; // 26-34: Skew matrix
    for (int i = 35; i <= 37; i++)
        WRITE_IAUX; // 35-37: Skew translation
    for (int i = 38; i <= 52; i++)
        WRITE_IAUX; // 38-52: Future use
    char map[5] = "MAP ";
    fwrite(map, 1, 4, fp); // 53: "MAP " to identify the file
    faux = (float)1.0;
    WRITE_FAUX; // 54: Machine stamp
    faux = (float)0.0;
    WRITE_FAUX; // 55: RMS deviation
    iaux =       0;
    WRITE_IAUX; // 56: Number of labels
    for (int i = 57; i <= 256; i++)
        WRITE_IAUX; // 57-256: Labels

    V->write(fp, false, VFLOAT);
    fclose(fp);
}

/* Is Xmipp image? --------------------------------------------------------- */
int Is_VolumeXmipp(const FileName &fn, bool skip_type_check, bool force_reversed)
{
    FILE *fp;
    int result = 0;
    headerXmipp header(headerXmipp::VOL_XMIPP);
    headerXmipp header1(headerXmipp::VOL_INT);

    // Open file
    if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "Is_VolumeXmipp: File " + fn + " not found");

    // Read header
    result = header.read(fp, skip_type_check, force_reversed);
    if (result == 1)
    {
        fclose(fp);
        return (headerXmipp::VOL_XMIPP);
    }


    result = header1.read(fp, skip_type_check, force_reversed);
    fclose(fp);
    return (headerXmipp::VOL_INT * result);
}

/* Is Xmipp image? --------------------------------------------------------- */
int Is_FourierVolumeXmipp(const FileName &fn, bool skip_type_check, bool force_reversed)
{
    FILE *fp;
    int result = 0;
    headerXmipp header(headerXmipp::VOL_FOURIER);

    // Open file
    if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "Is_FourierVolumeXmipp: File " + fn + " not found");

    // Read header
    result = header.read(fp, skip_type_check, force_reversed);
    fclose(fp);
    return result;
}

/* Get Volume size ---------------------------------------------------------- */
void GetXmippVolumeSize(const FileName &fn, int &Zdim, int &Ydim, int &Xdim)
{
    FILE *fp;
    int result;
    headerXmipp header(headerXmipp::VOL_XMIPP);

    // Open file
    if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "Is_ImageXmipp: File " + fn + " not found");

    // Read header
    result = header.read(fp, false, false);

    fclose(fp);
    if (result)
    {
        Zdim = header.iZdim();
        Ydim = header.iYdim();
        Xdim = header.iXdim();
    }
    else
    {
        Zdim = Ydim = Xdim = -1;
    }
}
