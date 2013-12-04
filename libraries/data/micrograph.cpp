/***************************************************************************
 *
 * Authors: Carlos Oscar (coss@cnb.csic.es)
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

#include "micrograph.h"
#include "args.h"
#include "metadata.h"
#include "mask.h"
#include "geometry.h"

#include <fstream>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
//#ifdef LINUX
//#include <unistd.h>
//#endif
Micrograph::Micrograph()
{
    auxI = new (Image<char> );
    IUChar = NULL;
    IShort = NULL;
    IUShort = NULL;
    IInt = NULL;
    IUInt = NULL;
    IFloat = NULL;
    stdevFilter = -1;
}
Micrograph::~Micrograph()
{
    delete (auxI);
    delete (IUChar);
    delete (IShort);
    delete (IUShort);
    delete (IInt);
    delete (IUInt);
    delete (IFloat);
}
/* Clear ------------------------------------------------------------------- */
void Micrograph::clear()
{
    single_particle.clear();
    coords.clear();
    fn_coords = fn_micrograph = "";
    X_window_size = Y_window_size = -1;
    fh_micrograph = -1;
    Xdim = Ydim = -1;
    datatype = -1;
    compute_transmitance = false;
    compute_inverse = false;
    delete (IUChar);
    delete (IShort);
    delete (IUShort);
    delete (IInt);
    delete (IUInt);
    delete (IFloat);
}

/* Open micrograph --------------------------------------------------------- */
void Micrograph::open_micrograph(const FileName &_fn_micrograph)
{
    clear();
    // Micrograph name
    fn_micrograph = _fn_micrograph;
    // Look for micrograph dimensions
    auxI->read(fn_micrograph, HEADER);

    auxI->getDimensions(Xdim, Ydim, Zdim, Ndim);
    if ((Zdim > 1) || (Ndim > 1))
        REPORT_ERROR(
            ERR_MULTIDIM_DIM,
            "Micrograph::open_micrograph: Only files with a single micrograph may be processed. Error reading " + fn_micrograph);
    auxI->MDMainHeader.getValue(MDL_DATATYPE, datatype);
    __offset = 0;
    auxI->clear();
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "x,y,z,n, datatype : "
    << Xdim << " "
    << Ydim << " "
    << Zdim << " "
    << Ndim << " "
    << datatype << " "
    <<std::endl;
#endif
#undef DEBUG
    //3,4,6
    // Open micrograph and map
    int result;
    switch (datatype)
    {
    case DT_UChar:
        IUChar = new (Image<unsigned char> );
        result = IUChar->readMapped(fn_micrograph, FIRST_IMAGE);
        pixelDesvFilter(IUChar->data, stdevFilter);
        break;
    case DT_UShort:
        IUShort = new (Image<unsigned short> );
        result = IUShort->readMapped(fn_micrograph, FIRST_IMAGE);
        pixelDesvFilter(IUShort->data, stdevFilter);
        break;
    case DT_Short:
        IShort = new (Image<short> );
        result = IShort->readMapped(fn_micrograph, FIRST_IMAGE);
        pixelDesvFilter(IShort->data, stdevFilter);
        break;
    case DT_Int:
        IInt = new (Image<int> );
        result = IInt->readMapped(fn_micrograph, FIRST_IMAGE);
        pixelDesvFilter(IInt->data, stdevFilter);
        break;
    case DT_UInt:
        IUInt = new (Image<unsigned int> );
        result = IUInt->readMapped(fn_micrograph, FIRST_IMAGE);
        pixelDesvFilter(IUChar->data, stdevFilter);
        break;
    case DT_Float:
        IFloat = new (Image<float> );
        result = IFloat->readMapped(fn_micrograph, FIRST_IMAGE);
        pixelDesvFilter(IFloat->data, stdevFilter);
        break;
    default:
        std::cerr << "Micrograph::open_micrograph: Unknown datatype "
        << datatype << std::endl;
        REPORT_ERROR(ERR_TYPE_INCORRECT, "ERROR");
        break;
    }
    //fh_micrograph = open(fn_micrograph.c_str(), O_RDWR, S_IREAD | S_IWRITE);
    if (result < 0)
        REPORT_ERROR(
            ERR_IO_NOTEXIST,
            (std::string)"Micrograph::open_micrograph: There is a "
            "problem opening " + fn_micrograph + "\nCheck that the file has write permission");
}

/* Close micrograph -------------------------------------------------------- */
void Micrograph::close_micrograph()
{
    switch (datatype)
    {
    case DT_UChar:
        delete (IUChar);
        IUChar = NULL;
        break;
    case DT_UShort:
        delete (IUShort);
        IUShort = NULL;
        break;
    case DT_Short:
        delete (IShort);
        IShort = NULL;
        break;
    case DT_Int:
        delete (IInt);
        IInt = NULL;
        break;
    case DT_UInt:
        delete (IUInt);
        IUInt = NULL;
        break;
    case DT_Float:
        delete (IFloat);
        IFloat = NULL;
        break;
    default:
        std::cerr << "Micrograph::close_micrograph: Unknown datatype "
        << datatype << std::endl;
        REPORT_ERROR(ERR_TYPE_INCORRECT, "ERROR");
        break;
    }
}
/* Get datatype detph -------------------------------------------------------- */
int Micrograph::getDatatypeDetph() const
{
    switch (datatype)
    {
    case DT_UChar:
        return (8 * sizeof(unsigned char));
    case DT_UShort:
        return (8 * sizeof(unsigned short));
    case DT_Short:
        return (8 * sizeof(short));
    case DT_UInt:
        return (8 * sizeof(unsigned int));
    case DT_Int:
        return (8 * sizeof(int));
    case DT_Float:
        return (8 * sizeof(float));
    default:
        std::cerr << "Micrograph::getDatatypeDetph: Unknown datatype "
        << datatype << std::endl;
        REPORT_ERROR(ERR_TYPE_INCORRECT, "ERROR");
        break;
    }
}

/* Save coordinates to disk ------------------------------------------------ */
void Micrograph::write_coordinates(int label, double minCost,
                                   const FileName &_fn_coords)
{
    std::ofstream fh;
    if (_fn_coords != "")
        fn_coords = _fn_coords;

    MetaData MD;
    MD.setComment((std::string) "Selected Coordinates for file " + fn_coords);
    int imax = coords.size();
    size_t id;
    for (int i = 0; i < imax; i++)
    {
        if (coords[i].valid && coords[i].cost > minCost
                && coords[i].label == label)
        {
            id = MD.addObject();
            MD.setValue(MDL_XCOOR, coords[i].X, id);
            MD.setValue(MDL_YCOOR, coords[i].Y, id);
        }
    }
    MD.write(fn_coords);
}

/* Read coordinates from disk ---------------------------------------------- */
void Micrograph::read_coordinates(int label, const FileName &_fn_coords)
{
    std::ifstream fh;
    int line_no = 0;
    std::string line;

    fn_coords = _fn_coords;

    MetaData MD;
    MD.read(fn_coords);
    line_no = MD.size();

    // Resize coordinate list and read
    coords.reserve(line_no);
    struct Particle_coords aux;
    aux.valid = true;
    aux.label = label;
    aux.cost = 1;

    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        MD.getValue(MDL_XCOOR, aux.X, __iter.objId); //aux.X=x;
        MD.getValue(MDL_YCOOR, aux.Y, __iter.objId); //aux.Y=y;
        coords.push_back(aux);
    }
}

/* Transform all coordinates ---------------------------------------------- */
void Micrograph::transform_coordinates(const Matrix2D<double> &M)
{
    Matrix1D<double> m(3);
    SPEED_UP_temps012;

    int imax = coords.size();
    for (int i = 0; i < imax; i++)
    {
        if (coords[i].valid)
        {
            VECTOR_R3(m, coords[i].X, coords[i].Y, 1);
            M3x3_BY_V3x1(m, M, m);
            coords[i].X = (int) XX(m);
            coords[i].Y = (int) YY(m);
        }
    }
}

/* Multiply coordinates by a constant -------------------------------------- */
void Micrograph::scale_coordinates(const double &c)
{
    int imax = coords.size();
    for (int i = 0; i < imax; i++)
    {
        if (coords[i].valid)
        {
            coords[i].X = (int) (coords[i].X * c);
            coords[i].Y = (int) (coords[i].Y * c);
        }
    }

}

/* Scissor ----------------------------------------------------------------- */
int Micrograph::scissor(const Particle_coords &P, MultidimArray<double> &result,
                        double Dmin, double Dmax, double scaleX, double scaleY,
                        bool only_check)
{
    if (X_window_size == -1 || Y_window_size == -1)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     "Micrograph::scissor: window size not set");
    if (datatype == DT_UChar)
        return templateScissor(*IUChar, P, result, Dmin, Dmax, scaleX, scaleY,
                               only_check);
    else if (datatype == DT_UShort)
        return templateScissor(*IUShort, P, result, Dmin, Dmax, scaleX, scaleY,
                               only_check);
    else if (datatype == DT_Short)
        return templateScissor(*IShort, P, result, Dmin, Dmax, scaleX, scaleY,
                               only_check);
    else if (datatype == DT_UInt)
        return templateScissor(*IUInt, P, result, Dmin, Dmax, scaleX, scaleY,
                               only_check);
    else if (datatype == DT_Int)
        return templateScissor(*IInt, P, result, Dmin, Dmax, scaleX, scaleY,
                               only_check);
    else if (datatype == DT_Float)
        return templateScissor(*IFloat, P, result, Dmin, Dmax, scaleX, scaleY,
                               only_check);
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT,
                     "Micrograph::scissor: unknown datatype");

}

/* Produce all images ------------------------------------------------------ */
void Micrograph::produce_all_images(int label, double minCost,
                                    const FileName &fn_rootIn, const FileName &fn_image, double ang,
                                    double tilt, double psi, bool rmStack)
{
    MetaData SF;
    Image<double> I;
    Micrograph *M;

    // Set Source image
    if (fn_image == "")
        M = this;
    else
    {
        M = new Micrograph;
        M->open_micrograph(fn_image/*, swapbyte*/);
        M->set_window_size(X_window_size, Y_window_size);
        M->set_transmitance_flag(compute_transmitance);
        M->set_inverse_flag(compute_inverse);
    }

    // Set scale for particles
    int MXdim, MYdim, thisXdim, thisYdim;
    M->size(MXdim, MYdim);
    this->size(thisXdim, thisYdim);
    double scaleX = (double) MXdim / thisXdim;
    double scaleY = (double) MYdim / thisYdim;

    // Compute max and minimum if compute_transmitance
    // or compute_inverse flags are ON
    double Dmax=0., Dmin=0.;
    if (compute_transmitance || compute_inverse)
    {
        (*this).computeDoubleMinMax(Dmin, Dmax);

        if (compute_transmitance)
        {
            if (Dmin > 1)
                Dmin = log10(Dmin);
            if (Dmax > 1)
                Dmax = log10(Dmax);
        }
    }
    // Scissor all particles
    if (ang != 0)
        std::cout << "Angle from Y axis to tilt axis " << ang << std::endl
        << "   applying appropriate rotation\n";
    int nmax = ParticleNo();
    FileName fn_aux;
    FileName _ext = fn_rootIn.getFileFormat();
    FileName fn_out;
    FileName fn_root = fn_rootIn.removeFileFormat().removeLastExtension();
    if (fn_rootIn.hasStackExtension())
        fn_out=fn_root.addExtension(_ext);
    else
    	fn_out=fn_rootIn.addExtension("stk");

    if (rmStack)
        fn_out.deleteFile();
    size_t ii = 0;
    size_t id;
    for (int n = 0; n < nmax; n++)
        if (coords[n].valid && coords[n].cost > minCost && coords[n].label == label)
        {
            fn_aux.compose(++ii, fn_out);
            id = SF.addObject();
            SF.setValue(MDL_IMAGE, fn_aux, id);
            SF.setValue(MDL_MICROGRAPH, M->fn_micrograph, id);
            SF.setValue(MDL_XCOOR, coords[n].X, id);
            SF.setValue(MDL_YCOOR, coords[n].Y, id);
            bool t = M->scissor(coords[n], I(), Dmin, Dmax, scaleX, scaleY);
            if (!t)
            {
                std::cout << "Particle " << fn_aux
                << " is very near the border, "
                << "corresponding image is set to blank\n";
                SF.setValue(MDL_ENABLED, -1, id);
            }
            else
                SF.setValue(MDL_ENABLED, 1, id);
            //  if (ang!=0) I().rotate(-ang);
            I.write(fn_out, ii, true, WRITE_APPEND);
        }
    SF.write(fn_out.withoutExtension() + ".xmd");


    // Free source image??
    if (fn_image != "")
    {
        M->close_micrograph();
        delete M;
    }
}

/* Search coordinate near a position --------------------------------------- */
int Micrograph::search_coord_near(int x, int y, int prec) const
{
    int imax = coords.size();
    int prec2 = prec * prec;
    for (int i = 0; i < imax; i++)
        if ((coords[i].X - x) * (coords[i].X - x)
                + (coords[i].Y - y) * (coords[i].Y - y) < prec2
                && coords[i].valid)
            return i;
    return -1;
}

/* Invalidate a coordinate ------------------------------------------------- */
void Micrograph::invalidate_coord(int n)
{
    if (n < 0 || n >= ParticleNo())
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,
                     "Micrograph::invalidate_coord: Index out of range");
    coords[n].valid = false;
}

/* Add coordinate ---------------------------------------------------------- */
int Micrograph::add_coord(int x, int y, int label, double cost)
{
    struct Particle_coords aux;
    aux.valid = true;
    aux.X = x;
    aux.Y = y;
    aux.label = label;
    aux.cost = cost;
    coords.push_back(aux);
    return coords.size() - 1;
}

/* Move last coordinate ---------------------------------------------------- */
void Micrograph::move_last_coord_to(int x, int y)
{
    if (coords.size() > 0)
    {
        coords.back().X = x;
        coords.back().Y = y;
    }
}

void Micrograph::resize(int Xdim, int Ydim, const FileName &filename)
{
    this->Xdim = Xdim;
    this->Ydim = Ydim;
    this->Zdim = 1;
    this->Ndim = 1;
    if (datatype == DT_UChar)
    {
        if (IUChar == NULL)
            IUChar = new Image<unsigned char>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IUChar->data.setMmap(true);
            IUChar->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == DT_UShort)
    {
        if (IUShort == NULL)
            IUShort = new Image<unsigned short int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IUShort->data.setMmap(true);
            IUShort->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == DT_Short)
    {
        if (IShort == NULL)
            IShort = new Image<short int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IShort->data.setMmap(true);
            IShort->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == DT_UInt)
    {
        if (IUInt == NULL)
            IUInt = new Image<unsigned int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IUInt->data.setMmap(true);
            IUInt->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == DT_Int)
    {
        if (IInt == NULL)
            IInt = new Image<int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IInt->data.setMmap(true);
            IInt->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == DT_Float)
    {
        if (IFloat == NULL)
            IFloat = new Image<float>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IFloat->data.setMmap(true);
            IFloat->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT, "Unknown datatype");

}
void Micrograph::write(const FileName &fileName, CastWriteMode castMode)
{
    if (datatype == DT_UChar)
    {
        IUChar->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == DT_UShort)
    {
        IUShort->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == DT_Short)
    {
        IShort->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == DT_UInt)
    {
        IUInt->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == DT_Int)
    {
        IInt->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == DT_Float)
    {
        IFloat->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT,
                     "Micrograph::set_val::(): unknown datatype");

}

/* Tilt pair aligner ------------------------------------------------------ */
TiltPairAligner::TiltPairAligner()
{
    clear();
}

void TiltPairAligner::clear()
{
    coordU.clear();
    coordT.clear();
    Au.clear();
    Bt.clear();
    Put.clear();
    Ptu.clear();
    Au.initZeros(3, 3);
    Bt.initZeros(3, 3);
    Nu = 0;
    m.resizeNoCopy(3);
}

void TiltPairAligner::addCoordinatePair(int _muX, int _muY, int _mtX,
                                        int _mtY)
{
    coordU.push_back(_muX);
    coordU.push_back(_muY);
    coordT.push_back(_mtX);
    coordT.push_back(_mtY);
    Nu++; // Number of particles

#ifdef _DEBUG

    std::cout << "Adding point U(" << U.X << "," << U.Y << ") T(" << T.X << ","
    << T.Y << ")\n";
    std::cout << "A at input" << Au << "B at input" << Bt;
#endif
    // Adjust untilted dependent matrix
    MAT_ELEM(Au,0, 0) += _muX * _muX;
    MAT_ELEM(Au,0, 1) += _muX * _muY;
    MAT_ELEM(Au,0, 2) += _muX;
    MAT_ELEM(Au,1, 0) = MAT_ELEM(Au,0, 1);
    MAT_ELEM(Au,1, 1) += _muY * _muY;
    MAT_ELEM(Au,1, 2) += _muY;
    MAT_ELEM(Au,2, 0) = MAT_ELEM(Au,0, 2);
    MAT_ELEM(Au,2, 1) = MAT_ELEM(Au,1, 2);
    MAT_ELEM(Au,2, 2) = Nu;

    // Adjust tilted dependent matrix
    MAT_ELEM(Bt,0, 0) += _mtX * _muX;
    MAT_ELEM(Bt,0, 1) += _mtY * _muX;
    MAT_ELEM(Bt,0, 2) = MAT_ELEM(Au,0, 2);
    MAT_ELEM(Bt,1, 0) += _mtX * _muY;
    MAT_ELEM(Bt,1, 1) += _mtY * _muY;
    MAT_ELEM(Bt,1, 2) = MAT_ELEM(Au,1, 2);
    MAT_ELEM(Bt,2, 0) += _mtX;
    MAT_ELEM(Bt,2, 1) += _mtY;
    MAT_ELEM(Bt,2, 2) = MAT_ELEM(Au,2, 2);

#ifdef _DEBUG

    std::cout << "A at output" << Au << "B at output" << Bt;
#endif


}

/* Adjust passing matrix --------------------------------------------------- */
void TiltPairAligner::adjustPassingMatrix(int _muX, int _muY, int _mtX,
        int _mtY)
{
    addCoordinatePair(_muX, _muY, _mtX, _mtY);
    if (Nu > 3)
    {
        solve(Au, Bt, Put);
        Put = Put.transpose();
        Ptu = Put.inv();
    }
}

/* Passing to tilted ------------------------------------------------------- */
void TiltPairAligner::passToTilted(int _muX, int _muY, int &_mtX, int &_mtY)
{


    if (Nu > 3)
    {
        SPEED_UP_temps012;
        VECTOR_R3(m, _muX, _muY, 1);
        M3x3_BY_V3x1(m, Put, m);

        _mtX = (int) XX(m);
        _mtY = (int) YY(m);
    }
    else
    {
        _mtX = _muX;
        _mtY = _muY;
    }

}

/* Passing to tilted ------------------------------------------------------- */
void TiltPairAligner::passToUntilted(int _mtX, int _mtY, int &_muX, int &_muY)
{
    if (Nu > 3)
    {
        SPEED_UP_temps012;

        VECTOR_R3(m, _mtX, _mtY, 1);
        M3x3_BY_V3x1(m, Ptu, m);
        _muX = (int) XX(m);
        _muY = (int) YY(m);
    }
    else
    {
        _muX = _mtX;
        _muY = _mtY;
    }
}

/* Compute tilting angle --------------------------------------------------- */
void TiltPairAligner::computeGamma()
{
#define TRIANGLE_NO 15000
#define MIN_AREA       15
#define MAX_AREA   250000
    gamma = 0;
    Matrix1D<int> iju(2), iku(2), ijt(2), ikt(2); // From i to j in untilted
    // From i to k in untilted
    // From i to j in tilted
    // From i to k in tilted
    int triang = 0; // Number of triangles considered
    int i, j, k, counter1;
    counter1 = 0;
    randomize_random_generator();
    long noCombinations;
    noCombinations = Nu * (Nu - 1) * (Nu - 2) / 6;
    while (triang < TRIANGLE_NO && counter1 < noCombinations)
    {
        counter1++;
        i = (int)round(rnd_unif(0, Nu - 1));
        j = (int)round(rnd_unif(0, Nu - 1));
        k = (int)round(rnd_unif(0, Nu - 1));

        // Compute area of triangle in untilted micrograph
        VECTOR_R2(iju, coordU[j] - coordU[i], coordU[j+1] - coordU[i+1]);
        VECTOR_R2(iku, coordU[k] - coordU[i], coordU[k+1] - coordU[i+1]);
        double untilted_area = fabs(dotProduct(iju, iku)/*/2*/);
        if (untilted_area < MIN_AREA
           )
            continue; // For numerical stability

        // Compute area of the same triangle in the tilted micrograph
        VECTOR_R2(ijt, coordT[j] - coordT[i], coordT[j+1] - coordT[i+1]);
        VECTOR_R2(ikt, coordT[k] - coordT[i], coordT[k+1] - coordT[i+1]);
        double tilted_area = fabs(dotProduct(ijt, ikt)/*/2*/);
        if (tilted_area < MIN_AREA
           )
            continue; // For numerical stability
        if (tilted_area > MAX_AREA
           )
            continue; // micrograph are not perfect
        // sheets so avoid
        // very far away particles

        // Now we know that tilted_area=untilted_area*cos(gamma)
        if (tilted_area > untilted_area)
            continue; // There are user errors
        // In the point selection
        gamma += acos(tilted_area / untilted_area);
        triang++;
    }
    gamma /= triang;
    gamma = RAD2DEG(gamma);
    if (triang < 100)
        std::cout << "Not many particles, tilt angle may not be accurate"
        << std::endl;
}

/* Compute alphas ---------------------------------------------------------- */
double matrix_fitness(double *p, void *prm)
{
    TiltPairAligner *aligner = (TiltPairAligner *) prm;
    Euler_angles2matrix(-p[1], p[3], p[2], aligner->pair_E);
    double retval = 0;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            double error = fabs(
                               MAT_ELEM(aligner->pair_E,i, j)
                               - MAT_ELEM(aligner->Put, i, j));
            retval += error * error;
        }
    return retval;
}


void TiltPairAligner::computeAngles(double &ualpha, double &talpha, double &ogamma)
{
    alpha_u = alpha_t = 0;
    Matrix1D<double> angles(3);
    angles.initZeros();
    double fitness;
    int iter;

    // Coarse search
    double *aux = angles.adaptForNumericalRecipes();
    double best_alpha_u = 0, best_alpha_t = 0, best_fit = 1e8;
    aux[3] = gamma;
    for (aux[1] = 0; aux[1] < 180; aux[1] += 5)
        for (aux[2] = 0; aux[2] < 180; aux[2] += 5)
        {
            double fit = matrix_fitness(aux, this);
            if (fit < best_fit)
            {
                best_fit = fit;
                best_alpha_u = aux[1];
                best_alpha_t = aux[2];
            }
        }
    angles.killAdaptationForNumericalRecipes(aux);
    angles(0) = best_alpha_u;
    angles(1) = best_alpha_t;
    angles(2) = gamma;

    // Fine search
    Matrix1D<double> steps(3);
    steps.initConstant(1);
    powellOptimizer(angles, 1, 3, &matrix_fitness, this, 0.001, fitness, iter,
                    steps, false);
    alpha_u = angles(0);
    alpha_t = angles(1);
    gamma = angles(2);

    ualpha = alpha_u;
    talpha = alpha_t;
    ogamma = gamma;
}
