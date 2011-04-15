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
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#ifdef LINUX
#include <unistd.h>
#endif
Micrograph::Micrograph()
{
    auxI = new(Image<char>);
    IUChar=NULL;
    IShort=NULL;
    IUShort=NULL;
    IInt=NULL;
    IUInt=NULL;
    IFloat=NULL;
    stdevFilter=-1;
}
Micrograph::~Micrograph()
{
    delete(auxI);
    delete(IUChar);
    delete(IShort);
    delete(IUShort);
    delete(IInt);
    delete(IUInt);
    delete(IFloat);
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
    delete(IUChar);
    delete(IShort);
    delete(IUShort);
    delete(IInt);
    delete(IUInt);
    delete(IFloat);
}

/* Open micrograph --------------------------------------------------------- */
void Micrograph::open_micrograph(const FileName &_fn_micrograph)
{
    clear();
    struct stat info;
    // Micrograph name
    fn_micrograph = _fn_micrograph;
    // Look for micrograph dimensions
    auxI->read(fn_micrograph, HEADER);

    auxI->getDimensions(Xdim,Ydim, Zdim, Ndim);
    if((Zdim >1 )|| (Ndim >1))
        REPORT_ERROR(ERR_MULTIDIM_DIM,"Micrograph::open_micrograph: Only files with a single micrograph may be processed. Error reading " + fn_micrograph );
    auxI->MDMainHeader.getValue(MDL_DATATYPE,datatype);
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
    case UChar:
        IUChar = new(Image<unsigned char>);
        result=IUChar->read(fn_micrograph,DATA, FIRST_IMAGE, true);
        pixelDesvFilter(IUChar->data, stdevFilter);
        break;
    case UShort:
        IUShort = new(Image<unsigned short>);
        result=IUShort->read(fn_micrograph,DATA, FIRST_IMAGE, true);
        pixelDesvFilter(IUShort->data, stdevFilter);
        break;
    case Short:
        IShort = new(Image< short>);
        result=IShort->read(fn_micrograph,DATA, FIRST_IMAGE, true);
        pixelDesvFilter(IShort->data, stdevFilter);
        break;
    case Int:
        IInt = new(Image< int>);
        result=IInt->read(fn_micrograph,DATA, FIRST_IMAGE, true);
        pixelDesvFilter(IInt->data, stdevFilter);
        break;
    case UInt:
        IUInt = new(Image< unsigned int>);
        result=IUInt->read(fn_micrograph,DATA, FIRST_IMAGE, true);
        pixelDesvFilter(IUChar->data, stdevFilter);
        break;
    case Float:
        IFloat = new(Image<float>);
        result=IFloat->read(fn_micrograph,DATA, FIRST_IMAGE, true);
        pixelDesvFilter(IFloat->data, stdevFilter);
        break;
    default:
        std::cerr << "Micrograph::open_micrograph: Unknown datatype " << datatype <<std::endl;
        REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR");
        break;
    }
    //fh_micrograph = open(fn_micrograph.c_str(), O_RDWR, S_IREAD | S_IWRITE);
    if (result < 0)
        REPORT_ERROR(ERR_IO_NOTEXIST, (std::string)"Micrograph::open_micrograph: There is a "
                     "problem opening " + fn_micrograph +
                     "\nCheck that the file has write permission");
}


/* Close micrograph -------------------------------------------------------- */
void Micrograph::close_micrograph()
{
    switch (datatype)
    {
    case UChar:
        delete(IUChar);
        IUChar=NULL;
        break;
    case UShort:
        delete(IUShort);
        IUShort=NULL;
        break;
    case Short:
        delete(IShort);
        IShort=NULL;
        break;
    case Int:
        delete(IInt);
        IInt=NULL;
        break;
    case UInt:
        delete(IUInt);
        IUInt=NULL;
        break;
    case Float:
        delete(IFloat);
        IFloat=NULL;
        break;
    default:
        std::cerr << "Micrograph::close_micrograph: Unknown datatype " << datatype <<std::endl;
        REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR");
        break;
    }
}
/* Get datatype detph -------------------------------------------------------- */
int Micrograph::getDatatypeDetph() const
{
    switch (datatype)
    {
    case UChar:
        return (8*sizeof(unsigned char));
    case UShort:
        return (8*sizeof(unsigned short));
    case Short:
        return (8*sizeof(short));
    case UInt:
        return (8*sizeof(unsigned int));
    case Int:
        return (8*sizeof(int));
    case Float:
        return (8*sizeof(float));
    default:
        std::cerr << "Micrograph::getDatatypeDetph: Unknown datatype " << datatype <<std::endl;
        REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR");
        break;
    }
}

/* Save coordinates to disk ------------------------------------------------ */
void Micrograph::write_coordinates(int label, double minCost, const FileName &_fn_coords)
{
    std::ofstream fh;
    if (_fn_coords != "")
        fn_coords = _fn_coords;

    MetaData MD;
    MD.setComment((std::string)"Selected Coordinates for file " + fn_coords);
    int imax = coords.size();
    size_t id;
    for (int i = 0; i < imax; i++)
    {
        if (coords[i].valid && coords[i].cost>minCost && coords[i].label==label)
        {
            id = MD.addObject();
            MD.setValue(MDL_XINT,coords[i].X,id);
            MD.setValue(MDL_YINT,coords[i].Y,id);
        }
    }
    MD.write(fn_coords);
}

/* Read coordinates from disk ---------------------------------------------- */
void Micrograph::read_coordinates(int label, const FileName &_fn_coords)
{
    std::ifstream  fh;
    int            line_no = 0;
    std::string    line;

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
        int x,y;
        MD.getValue(MDL_XINT,aux.X,__iter.objId); //aux.X=x;
        MD.getValue(MDL_YINT,aux.Y,__iter.objId); //aux.Y=y;
        coords.push_back(aux);
    }
}

/* Transform all coordinates ---------------------------------------------- */
void Micrograph::transform_coordinates(const Matrix2D<double> &M)
{
    Matrix1D<double> m(3);
    SPEED_UP_temps;

    int imax = coords.size();
    for (int i = 0; i < imax; i++)
    {
        if (coords[i].valid)
        {
            VECTOR_R3(m,coords[i].X, coords[i].Y,1);
            M3x3_BY_V3x1(m, M, m);
            coords[i].X=(int)XX(m);
            coords[i].Y=(int)YY(m);
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
            coords[i].X = (int)coords[i].X*c;
            coords[i].Y = (int)coords[i].Y*c;
        }
    }

}

/* Scissor ----------------------------------------------------------------- */
int Micrograph::scissor(const Particle_coords &P, MultidimArray<double> &result,
                        double Dmin, double Dmax, double scaleX, double scaleY,
                        bool only_check)
{
    if (X_window_size == -1 || Y_window_size == -1)
        REPORT_ERROR(ERR_MULTIDIM_SIZE, "Micrograph::scissor: window size not set");

    if (datatype == UChar)
        return templateScissor(*IUChar,P,result,Dmin,Dmax,scaleX,scaleY,only_check);
    else if (datatype == UShort)
        return templateScissor(*IUShort,P,result,Dmin,Dmax,scaleX,scaleY,only_check);
    else if (datatype == Short)
        return templateScissor(*IShort,P,result,Dmin,Dmax,scaleX,scaleY,only_check);
    else if (datatype == UInt)
        return templateScissor(*IUInt,P,result,Dmin,Dmax,scaleX,scaleY,only_check);
    else if (datatype == Int)
        return templateScissor(*IInt,P,result,Dmin,Dmax,scaleX,scaleY,only_check);
    else if (datatype == Float)
        return templateScissor(*IFloat,P,result,Dmin,Dmax,scaleX,scaleY,only_check);
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT, "Micrograph::scissor: unknown datatype");

}

/* Produce all images ------------------------------------------------------ */
void Micrograph::produce_all_images(int label, double minCost, const FileName &fn_root,
                                    int starting_index, const FileName &fn_image,
                                    double ang, double tilt, double psi,
                                    bool rmStack)
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
    double scaleX = (double)MXdim / thisXdim;
    double scaleY = (double)MYdim / thisYdim;

    // Compute max and minimum if compute_transmitance
    // or compute_inverse flags are ON
    double Dmax, Dmin;
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
    int i = starting_index;
    int nmax = ParticleNo();
    FileName fn_aux;
    FileName fn_out = fn_root+".stk";
    if(exists(fn_out) && rmStack)
    {
        unlink(fn_out.c_str());
    }
    size_t ii=0;
    size_t id;
    for (int n = 0; n < nmax; n++)
        if (coords[n].valid && coords[n].cost>minCost && coords[n].label == label)
        {
            fn_aux.compose(++ii,fn_out);
            id = SF.addObject();
            SF.setValue( MDL_IMAGE, fn_aux, id);
            SF.setValue( MDL_MICROGRAPH, M->fn_micrograph, id);
            SF.setValue( MDL_XINT, coords[n].X, id);
            SF.setValue( MDL_YINT, coords[n].Y, id);
            bool t = M->scissor(coords[n], I(), Dmin, Dmax, scaleX, scaleY);
            if (!t)
            {
                std::cout << "Particle " << fn_aux << " is very near the border, "
                << "corresponding image is set to blank\n";
                SF.setValue( MDL_ENABLED, -1, id);
            }
            else
                SF.setValue( MDL_ENABLED, 1, id);
            //  if (ang!=0) I().rotate(-ang);
            I.write(fn_out,ii,true,WRITE_APPEND);
        }
    SF.write(fn_root + ".sel");

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
        if ((coords[i].X - x)*(coords[i].X - x) + (coords[i].Y - y)*(coords[i].Y - y) < prec2
            && coords[i].valid)
            return i;
    return -1;
}

/* Invalidate a coordinate ------------------------------------------------- */
void Micrograph::invalidate_coord(int n)
{
    if (n < 0 || n >= ParticleNo())
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "Micrograph::invalidate_coord: Index out of range");
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
    if (datatype == UChar)
    {
        if (IUChar == NULL)
            IUChar = new Image<unsigned char>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IUChar->data.setMmap(true);
            IUChar->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == UShort)
    {
        if (IUShort == NULL)
            IUShort = new Image<unsigned short int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IUShort->data.setMmap(true);
            IUShort->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == Short)
    {
        if (IShort == NULL)
            IShort = new Image<short int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IShort->data.setMmap(true);
            IShort->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == UInt)
    {
        if (IUInt == NULL)
            IUInt = new Image<unsigned int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IUInt->data.setMmap(true);
            IUInt->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == Int)
    {
        if (IInt == NULL)
            IInt = new Image<int>(Xdim, Ydim, 1, 1, filename);
        else
        {
            IInt->data.setMmap(true);
            IInt->data.resize(1, 1, Ydim, Xdim);
        }
    }
    else if (datatype == Float)
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
    if (datatype == UChar)
    {
        IUChar->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == UShort)
    {
        IUShort->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == Short)
    {
        IShort->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == UInt)
    {
        IUInt->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == Int)
    {
        IInt->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else if (datatype == Float)
    {
        IFloat->write(fileName, FIRST_IMAGE, false, WRITE_OVERWRITE, castMode);
    }
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT, "Micrograph::set_val::(): unknown datatype");

}
