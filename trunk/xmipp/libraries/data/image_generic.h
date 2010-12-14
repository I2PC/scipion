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

#ifndef IMAGE_GENERIC_H_
#define IMAGE_GENERIC_H_

#include "image.h"
#include "multidim_array_generic.h"

/// @addtogroup Images

//@{

/**
 *  ImageGeneric class to handle images with independence of data type
 **/
class ImageGeneric
{
public:
    DataType datatype;
    MultidimArrayGeneric * data;

private:
    ImageBase* image;

public:

    /** Empty constructor.
     *
     * No internal image class is declared.
     */
    ImageGeneric()
    {
        init();
    }
    /** Constructor passing the data type of the image.
     *
     * Defines the data type of the image and then declares the internal image class.
     */
    ImageGeneric(DataType _datatype)
    {
        datatype = _datatype;
        setImageBase();
    }
    /** Destructor.
     */
    ~ImageGeneric();

    /** Initialize the parameters.
     */
    void init()
    {
        image = NULL;
        datatype = Unknown_Type;
    }

    /** Clear the parameters and initialize them.
     */
    void clear()
    {
        if (image != NULL)
        {
            image->clear();
            delete image;
            init();
        }
    }

    /** Get Image dimensions
    */
    void getDimensions(int &Xdim, int &Ydim, int &Zdim, unsigned long &Ndim) const
    {
        image->getDimensions(Xdim, Ydim, Zdim, Ndim);
    }
    void getDimensions(int &Xdim, int &Ydim, int &Zdim) const
    {
        unsigned long Ndim;
        image->getDimensions(Xdim, Ydim, Zdim, Ndim);
    }

    /** Get Euler angles from image header
    */
    void getEulerAngles(double &rot, double &tilt, double &psi,
                        long int n = 0)
    {
        image->getEulerAngles(rot, tilt, psi, n);
    }

    /** Set the data type for the generic image
     */
    void setDatatype(DataType _datatype)
    {
        datatype = _datatype;
        setImageBase();
    }

    /** Read image from file.
     */
    int read(const FileName &name, bool readdata=true, int select_img = -1,
             bool apply_geo = false, bool only_apply_shifts = false,
             MDRow * row = NULL, bool mapData = false);
    /** Write image to file.
    */
    void write(const FileName &name="", int select_img=-1, bool isStack=false,
               int mode=WRITE_OVERWRITE,bool adjust=false)
    {
        image->write(name,select_img,isStack,mode,adjust);
    }

    /* Create an empty image file of format given by filename and map it to memory.
     */
    void newMappedFile(int Xdim, int Ydim, int Zdim, int Ndim, FileName _filename)
    {
        image->newMappedFile(Xdim,Ydim,Zdim,Ndim,_filename);
    }

    /* MultidimArrayGeneric data access
     */
    MultidimArrayGeneric operator()()
    {
        return *data;
    }

    const MultidimArrayGeneric* operator()() const
    {
        return data;
    }

private:
    /* Declare the internal image class with the stored data type
     */
    void setImageBase();

}
;

//@}

#endif /* IMAGE_GENERIC_H_ */
