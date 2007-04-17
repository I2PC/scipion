/***************************************************************************
 *
 * Authors: Lorenzo Zampighi and Nelson Tang
 *
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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

#ifndef XMIPPIMAGIC_HH
#define XMIPPIMAGIC_HH

#include <vector>

#include "xmippImages.hh"

#define GCC_VERSION (__GNUC__ * 10000 \
    + __GNUC_MINOR__ * 100 \
    + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
    #include <sstream>
#else
    #include <strstream.h>
#endif

/// @defgroup Imagic IMAGIC Images

/** Extensions for the Imagic header and image files
 * @ingroup Imagic
 */
static const char* IMAGIC_HEADER_EXT = "hed";
static const char* IMAGIC_IMAGE_EXT  = "img";

/** Length of the tag
 * @ingroup Imagic
 */
static const short unsigned IMAGIC_TAG_LEN = 7;

/** Separator character
 * @ingroup Imagic
 */
static const char IMAGIC_TAG_SEP = ':';

/** Types of Imagic files
 * @ingroup Imagic
 * 
 * Valid types are IMAGIC_REAL, IMAGIC_INTG, IMAGIC_PACK, IMAGIC_COMP.
 */
enum ImageImagicType { IMAGIC_REAL, IMAGIC_INTG, IMAGIC_PACK, IMAGIC_COMP };

/** Structure that holds Imagic image information
 * @ingroup Imagic
 */
struct ImageImagicInfo
{
    unsigned int num_img;
    unsigned int xsize, ysize;
    vector< ImageImagicType > img_types;
};

/** Imagic Image class
 * @ingroup Imagic
 */
template<typename T>
class ImageImagicT: public ImageT< T >
{
public:
    /** Empty constructor
     */
    ImageImagicT() : ImageT< T >()
    {
        name_parsed = false;
    }
    
    /** Constructor with size
     */
    ImageImagicT(int Ydim, int Xdim) : ImageT< T >(Ydim, Xdim)
    {
        name_parsed = false;
    }
    
    /** Constructor with image name
     */
    ImageImagicT(FileName _name) : ImageT< T >(_name)
    {
        name_parsed = false;
    }

    /** Copy constructor
     */
    ImageImagicT(const ImageImagicT& I) : ImageT< T >(I)
    {
        if ((name_parsed = I.name_parsed))
        {
            hedfname = I.hedfname;
            imgfname = I.imgfname;
            imgnum = I.imgnum;
        }
    }

    /** Rename
     */
    virtual void rename(FileName newName)
    {
        if (newName != ImageT< T >::fn_img)
        {
            ImageT< T >::rename (newName);
            name_parsed = false;
        }
    }

    /** Clear
     */
    virtual void clear()
    {
        ImageT< T >::clear();
        name_parsed = false;
    }

    /** Assignment operator
     */
    ImageImagicT< T >& operator=(const ImageImagicT< T >& I)
    {
        if (this != &I)
        {
            this->ImageT< T >::operator=(I);
            if ((name_parsed = I.name_parsed))
            {
                hedfname = I.hedfname;
                imgfname = I.imgfname;
                imgnum = I.imgnum;
            }
        }
        return (*this);
    }

    /** Assignment operator
     */
    ImageImagicT< T >& operator=(const ImageT< T >& I)
    {
        if (this != &I)
        {
            this->ImageT< T >::operator=(I);
            name_parsed = false;
        }
        return (*this);
    }

    /** Assignment operator
     */
    template<typename T1>
    ImageImagicT< T >& operator=(const matrix2D< T1 >& m)
    {
        if (&ImageT< T >::img != (matrix2D< T >*)& m)
        {
            this->ImageT< T >::operator=(m);
            name_parsed = false;
        }
        return (*this);
    }

    /** Show
     */
    friend ostream& operator<<(ostream& out, const ImageImagicT< T >& I)
    {
        out << (ImageT< T >&) I;
        
        /* COSS: These functions are not defined
           out << "IMAGIC header fname: " << get_hedfname() << endl;
           out << "IMAGIC image fname:  " << get_imgfname() << endl;
           out << "image number:        " << get_imgnum() << endl;
        */
        
        return (out);
    }

    /** Read file
     */
    bool read(const FileName& name)
    {
        rename(name);
        ImageImagicInfo img_info = ImagicGetImgInfo(getHedFname());

        FileName img_fname = getImgFname();
        if (img_fname == "")
            REPORT_ERROR (1501, "ImageImagic::read: File " + name +
                          " doesn't seem fit Imagic format");
        
        FILE* img_fh;
        if ((img_fh = fopen(img_fname.c_str(), "rb")) == NULL)
            REPORT_ERROR (1501, "ImageImagic::read: IMAGIC file " + img_fname +
                          " not found");

        const int imgnum = getImgNum();
        const size_t img_offset = img_info.xsize * img_info.ysize;
        
        // Read the image data
        ImageT< T >::img.resize(img_info.ysize, img_info.xsize);
        const bool reversed = false;
        switch (img_info.img_types[imgnum])
        {
        case IMAGIC_REAL:
            {
                float data;
                const unsigned size = 4;
                fseek (img_fh, imgnum * size * img_offset, SEEK_SET);
                FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(ImageT< T >::img)
                {
                    FREAD(&data, size, 1, img_fh, reversed);
                    MULTIDIM_ELEM(ImageT< T >::img, i) = data;
                }
                break;
            }
        case IMAGIC_INTG:
            {
                short int data;
                const unsigned size = 2;
                fseek(img_fh, imgnum * size * img_offset, SEEK_SET);
                FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(ImageT< T >::img)
                {
                    FREAD(&data, size, 1, img_fh, reversed);
                    MULTIDIM_ELEM(ImageT< T >::img, i) = data;
                }
                break;
            }
        case IMAGIC_PACK:
            {
                unsigned char data;
                const unsigned size = 1;
                fseek(img_fh, imgnum * size * img_offset, SEEK_SET);
                FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(ImageT< T >::img)
                {
                    FREAD(&data, size, 1, img_fh, reversed);
                    MULTIDIM_ELEM(ImageT< T >::img, i) = data;
                }
                break;
            }
        default:
            REPORT_ERROR(1501,
                "ImageImagicType not supported for this imgtype!");
            break;
        }
        fclose (img_fh);
        return (true);
    }

    /** Write file
     * 
     * FIXME Not implemented
     */
    void write(const FileName& name="", bool reversed=FALSE,
                Image_Type image_type=IBYTE)
    {
        REPORT_ERROR (1503, "ImageImagic::write: can't directly save");
    }

    /** Low level write
     * 
     * FIXME Not implemented
     */
    void write(FILE*& fh, bool reversed, Image_Type image_type)
    {
        REPORT_ERROR (1503, "ImageImagic::write: can't directly save");
    }

    /** Get Header name
     */
    const FileName& getHedFname()
    {
        parseFname();
        return (hedfname);
    }
    
    /** Get Image name
     */
    const FileName& getImgFname()
    {
        parseFname();
        return (imgfname);
    }
    
    /** Get image number
     */
    int getImgNum()
    {
        parseFname();
        return (imgnum);
    }

    virtual void parseFname()
    {
        // Not meant to be called directly, but needs to
        // be kept public
        if (!name_parsed)
        {
            hedfname = "";
            imgfname = "";
            imgnum = -1;
            
            // Look for special IMAGIC format: 'imagic:<hedfile>:<imgnum>'
            if (ImageT< T >::fn_img.find(IMAGIC_TAG) != string::npos)
            {
                const string::size_type imgnumpos = ImageT< T >::fn_img.rfind(
                    IMAGIC_TAG_SEP);
                if (imgnumpos > IMAGIC_TAG_LEN)
                {
                    hedfname = ImageT< T >::fn_img.substr(IMAGIC_TAG_LEN,
                        imgnumpos - IMAGIC_TAG_LEN);
                    imgfname = hedfname.substitute_extension(IMAGIC_HEADER_EXT,
                               IMAGIC_IMAGE_EXT);
                    imgnum = atoi((
                        ImageT< T >::fn_img.substr(imgnumpos+1)).c_str());
                }
            }
            name_parsed = true;
        }
    }


protected:
    bool name_parsed;
    FileName hedfname, imgfname;
    int imgnum;
};

/** Alias for ImageImagic type
 * @ingroup Imagic
 */
typedef ImageImagicT< double > ImageImagic;

/** Creates a string of the format 'imagic:hedfname:num' suitable for use as an
 * image name
 * @ingroup Imagic
 */
inline string ImagicMakeName(const char* hed_fname, unsigned int imgnum)
{
#if GCC_VERSION < 30300
    char aux[15];
    ostrstream ss(aux, sizeof(aux));
#else
    ostringstream ss;
#endif

    ss << IMAGIC_TAG << hed_fname << IMAGIC_TAG_SEP << imgnum;
    return (ss.str());
}

/** Looks at an IMAGIC header file for information about the images it
 * contains
 * @ingroup Imagic
 */
const ImageImagicInfo ImagicGetImgInfo(const FileName& hed_fname);

/** Creates a new Imagic header/image file pair, filling it with the data
 * pointed to by the vector parameter
 * @ingroup Imagic
 */
template<typename T>
bool ImagicWriteImagicFile(const FileName& hed_fname,
                            const vector< ImageT< T >* >& imgs,
                            ImageImagicType img_type = IMAGIC_REAL);

#endif
