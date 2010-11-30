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


//// abstract base class
//class TFunctorImage
//{
//public:
//
//    // two possible functions to call member function. virtual cause derived
//    // classes will use a pointer to an object and a pointer to a member function
//    // to make the function call
//    virtual int read(const FileName &name, bool readdata=true, int select_img = -1,
//                     bool apply_geo = false, bool only_apply_shifts = false,
//                     MDRow * row = NULL, bool mapData = false)
//    {
//
//    }   ;  // call using operator
//    //    virtual void write(const FileName &name="", int select_img=-1, bool isStack=false,
//    //                       int mode=WRITE_OVERWRITE,bool adjust=false);        // call using function
//};


//// derived template class
//template <typename TClass>
//class TSpecificFunctorImage : public TFunctorImage
//{
//private:
//    TClass* pt2Object;                  // pointer to object
//
//public:
//
//    // constructor - takes pointer to an object and pointer to a member and stores
//    // them in two private variables
//    TSpecificFunctorImage(TClass* _pt2Object)
//    {
//        pt2Object = _pt2Object;
//    };
//
//    // override operator "()"
//    virtual int read(const FileName &name, bool readdata=true, int select_img = -1,
//                     bool apply_geo = false, bool only_apply_shifts = false,
//                     MDRow * row = NULL, bool mapData = false)
//    {
//        pt2Object->read(name,readdata,select_img,apply_geo, only_apply_shifts,row,
//                        mapData);
//    };
//    // execute member function
//
//    // override function "Call"
//    //    void write(const FileName &name="", int select_img=-1, bool isStack=false,
//    //                       int mode=WRITE_OVERWRITE,bool adjust=false)
//    //    {
//    //        pt2Object->write(name,select_img,isStack,mode,adjust);
//    //    }
//    //    ;             // execute member function
//};


class ImageGeneric
{
public:
    DataType datatype;
    MultidimArrayGeneric * data;

private:

    ImageBase* image;

public:

    ImageGeneric()
    {
        init();
    }

    ImageGeneric(DataType _datatype)
    {
        datatype = _datatype;
        setImageBase();
    }

    ~ImageGeneric()
    {}

    void init()
    {
        image = NULL;
        datatype = Unknown_Type;
    }

    void clear()
    {
        if (image != NULL)
        {
            image->clear();
            delete image;
            init();
        }
    }

    void setDatatype(DataType _datatype)
    {
        datatype = _datatype;
        setImageBase();
    }

    int read(const FileName &name, bool readdata=true, int select_img = -1,
             bool apply_geo = false, bool only_apply_shifts = false,
             MDRow * row = NULL, bool mapData = false);

    void write(const FileName &name="", int select_img=-1, bool isStack=false,
               int mode=WRITE_OVERWRITE,bool adjust=false)
    {
        image->write(name,select_img,isStack,mode,adjust);
    }


    void newMappedFile(int Xdim, int Ydim, int Zdim, int Ndim, FileName _filename)
    {
        image->newMappedFile(Xdim,Ydim,Zdim,Ndim,_filename);
    }


private:

    void setImageBase();

}
;


#endif /* IMAGE_GENERIC_H_ */
