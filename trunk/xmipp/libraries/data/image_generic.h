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

#include "image_base.h"
#include "multidim_array_generic.h"

/// @addtogroup Images

//@{

/**
 *  ImageGeneric class to handle images with independence of data type
 **/
class ImageGeneric
{

public:
    ImageBase* image;
    DataType datatype;
    MultidimArrayGeneric * data;


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
    ImageGeneric(DataType _datatype);

    /** Destructor.
     */
    ~ImageGeneric();

    /** Initialize the parameters.
     */
    void init();

    /** Clear the parameters and initialize them.
     */
    void clear();

    /** Clear the image header
     */
    void clearHeader()
    {
        image->clearHeader();
    }

    /** Init geometry transformation with defaults values
     */
    void initGeometry(const size_t n = 0)
    {
        image->initGeometry(n);
    }

    /** Return geometry row
     */
    MDRow& getGeometry(const size_t n = 0)
    {
        return image->MD[n];
    }

    /** Get Image dimensions
    */
    void getDimensions(int &Xdim, int &Ydim, int &Zdim, size_t &Ndim) const
    {
        image->getDimensions(Xdim, Ydim, Zdim, Ndim);
    }
    void getDimensions(int &Xdim, int &Ydim, int &Zdim) const
    {
        size_t Ndim;
        image->getDimensions(Xdim, Ydim, Zdim, Ndim);
    }

    /** Get number of elements in image
     */
    size_t getSize() const
    {
        return NZYXSIZE(*(data->im));
    }

    /** Get Euler angles from image header
    */
    void getEulerAngles(double &rot, double &tilt, double &psi,
                        size_t n = 0)
    {
        image->getEulerAngles(rot, tilt, psi, n);
    }

    /** Get Tilt angle from image header
    */
    double tilt(const size_t n = 0) const
    {
        return image->tilt(n);
    }

    /** Set the data type for the generic image
     */
    void setDatatype(DataType _datatype);

    /** Set image dataMode */
    void setDataMode(DataMode mode)
    {
        image->setDataMode(mode);
    }

    /** Get the data type
     */
    DataType getDatatype()const
    {
        return datatype;
    }

    /** Check if image is mapped on file
      */
    inline bool isMapped()
    {
        return image->isMapped();
    }

    /** Read image from file.
     */
    int read(const FileName &name, DataMode datamode = DATA, size_t select_img = ALL_IMAGES,
             bool mapData = false);

    /** Read image from file with a header applied.
     */
    int readApplyGeo(const FileName &name, const MDRow &row, bool only_apply_shifts = false, DataMode datamode = DATA, size_t select_img = ALL_IMAGES);

    /** Read image from file.
     */
    int readApplyGeo(const FileName &name, const MetaData &md, size_t objId, bool only_apply_shifts = false,
        DataMode datamode = DATA, size_t select_img = ALL_IMAGES);

    /** Read an image from metadata, filename is taken from MDL_IMAGE */
    int readApplyGeo(const MetaData &md, size_t objId, bool only_apply_shifts = false,
        DataMode datamode = DATA, size_t select_img = ALL_IMAGES );

    /** Read image mapped from file.
     */
    inline int readMapped(const FileName &name, int select_img = 0)
    {
        read(name, DATA, select_img, true);
    }

    /* Read an image with a lower resolution as a preview image.
    * If Zdim parameter is not passed, then all slices are rescaled.
    */
    int readPreview(const FileName &name, int Xdim, int Ydim, int Zdim = -1, size_t select_img = FIRST_IMAGE)
    {
        image->readPreview(name, Xdim, Ydim, Zdim, select_img);
    }

    /** Write image to file.
    */
    inline void write(const FileName &name="", size_t select_img = ALL_IMAGES, bool isStack=false,
                      int mode=WRITE_OVERWRITE,bool adjust=false)
    {
        image->write(name,select_img,isStack,mode,adjust);
    }

    /* Create an empty image file of format given by filename and map it to memory.
     */
    inline  void newMappedFile(int Xdim, int Ydim, int Zdim, int Ndim, FileName _filename,
                               bool createTempFile=false)
    {
        image->newMappedFile(Xdim,Ydim,Zdim,Ndim,_filename,createTempFile);
    }

    /* MultidimArrayGeneric data access
     */
    inline MultidimArrayGeneric& operator()()
    {
        return *data;
    }

    inline const MultidimArrayGeneric& operator()() const
    {
        return *data;
    }

    /** Get pixel value
     */
    inline double getPixel(unsigned long n, int k, int i, int j) const
    {
#define GETVALUE(type) return NZYX_ELEM(*(MultidimArray<type>*)data->im,n,k,i,j)
        SWITCHDATATYPE(datatype,GETVALUE);
#undef GETVALUE

    }

    /** Set pixel value
     */
    inline void setPixel(unsigned long n, int k, int i, int j, double value) const
    {
#define GETVALUE(type) NZYX_ELEM(*(MultidimArray<type>*)data->im,n,k,i,j) = value
        SWITCHDATATYPE(datatype,GETVALUE);
#undef GETVALUE

    }

    void print() const;

private:
    /* Declare the internal image class with the stored data type
     */
    DataType getImageType(const FileName &imgName);
}
;

//@}

#endif /* IMAGE_GENERIC_H_ */
