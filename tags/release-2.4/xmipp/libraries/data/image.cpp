/***************************************************************************
 *
 * Authors: Pedro Antonio de Alarc�n (pedro@cnb.csic.es)
 *          Carlos Oscar S. Sorzano
 *          Alberto Pascual Montano
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medistd::cine
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "image.h"
#include "selfile.h"
#include "volume.h"

/* Is Xmipp image? --------------------------------------------------------- */
int Is_ImageXmipp(const FileName &fn, bool skip_type_check,
                  bool reversed)
{
    FILE *fp;
    int result;
    headerXmipp header(headerXmipp::IMG_XMIPP);

    // Open file
    if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "Is_ImageXmipp: File " + fn + " not found");

    // Read header
    result = header.read(fp, skip_type_check, reversed);

    fclose(fp);

    return result;
}

/* Is Fourier Xmipp image? ------------------------------------------------- */
int Is_FourierImageXmipp(const FileName &fn, bool skip_type_check,
                         bool reversed)
{
    FILE *fp;
    int result;
    headerXmipp header(headerXmipp::IMG_FOURIER);

    // Open file
    if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "Is_FourierImageXmipp: File " + fn + " not found");

    // Read header
    result = header.read(fp, skip_type_check, reversed);

    fclose(fp);

    return result;
}

/* Get Image size ---------------------------------------------------------- */
void GetXmippImageSize(const FileName &fn, int &Ydim, int &Xdim)
{
    FILE *fp;
    int result;
    headerXmipp header(headerXmipp::IMG_XMIPP);

    // Open file
    if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "Is_ImageXmipp: File " + fn + " not found");

    // Read header
    result = header.read(fp, false, false);

    fclose(fp);
    if (result)
    {
        Ydim = header.iYdim();
        Xdim = header.iXdim();
    }
    else
    {
        Ydim = Xdim = -1;
    }
}

/* Convert a Xmipp Image into  a Fourier Xmipp Image -------------------------------*/
// A simple copy of real numbers is done
void ImageXmipp_to_FourierImageXmipp(ImageXmipp &I, FourierImageXmipp &F)
{
    // Adjust the size of the Fourier Image
    F().resize(I().rowNumber(), I().colNumber());
    // And copy
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
    {
        F(i, j) = I(i, j);
    }
}

/* Convert a Fourier Xmipp Image  into a Xmipp Image -------------------------------*/
// A simple copy of real parts is done
void FourierImageXmipp_to_ImageXmipp(FourierImageXmipp &F, ImageXmipp &I)
{
    // Adjust the size of the Fourier Image
    I().resize(F().rowNumber(), F().colNumber());
    // And copy
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I())
    {
        I(i, j) = F(i, j).real();
    }
}

// Read stack from stack file ----------------------------------------------
bool ImageXmippStack::readFromStack(const FileName& name,
    bool reversed, bool apply_geo, bool only_apply_shifts)
{
    headerXmipp header;
    FILE* fp;
    bool ret;

    fn=name;
    if ((fp = fopen(name.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, (std::string) "ImageXmippStack::read: File " +
                     name + " not found");

    // Read header
    if (!header.read(fp, false, reversed))
        REPORT_ERROR(1502, "ImageXmippStack::read: File " + name +
                     " is not a valid Xmipp file");

    // Compute the number of images in the stack.
    long int headerSize = header.get_header_size();
    fseek(fp, 0, SEEK_END);
    long int fileSize = ftell(fp);
    long int imageSize = header.Ydim() * header.Xdim() * 4;
    int Nimgs = (fileSize-headerSize)/(imageSize+headerSize);    

    // Read all images
    fseek(fp,headerSize,SEEK_SET);
    for (int i=0; i<Nimgs; i++) 
    {
        ImageXmipp I;
        if (!I.getHeader().read(fp, true, reversed, true))
            REPORT_ERROR(1502,
                (std::string)"ImageXmippStack::readFromStack: Cannot read image "+
                integerToString(i));
        if (!I.readImageContent(fp, apply_geo, only_apply_shifts))
            REPORT_ERROR(1502,
                (std::string)"ImageXmippStack::readFromStack: Cannot read (2) image "
                +integerToString(i));
        I.rename(name.insert_before_extension(integerToString(i,5)));
        stack.push_back(I);
    }

    // Close file
    fclose(fp);
    return (ret);
}
    
// Write as stack file.-----------------------------------------------------
void ImageXmippStack::writeAsStack(const FileName& name, bool force_reversed)
{
    FILE* fp;
    if (name != "")
        fn=name;

    if ((fp = fopen(fn.c_str(), "wb")) == NULL)
        REPORT_ERROR(1503, (std::string) "ImageXmippStack::write: File " + fn +
                     " cannot be written");

    getImage(0).adjust_header();
    bool reversed = (force_reversed) ?
        !getImage(0).getHeader().reversed():getImage(0).getHeader().reversed();
    headerXmipp overallHeader;
    overallHeader=getImage(0).getHeader();
    overallHeader.fIangle()=0;
    overallHeader.fStack()=1;
    overallHeader.fMaxim()=ImgNo();
    overallHeader.fLastidx()=0;
    overallHeader.write(fp, force_reversed);
    for (int i=0; i<ImgNo(); i++)
    {
        getImage(i).getHeader().fImgnum()=i+1;
        getImage(i).write(fp, reversed);
    }
    
    fclose(fp);
}
    
// Write as selfile --------------------------------------------------------
void ImageXmippStack::writeAsSelFile(const FileName& rootname,
    bool force_reversed)
{
    SelFile SF;
    for (int i=0; i<ImgNo(); i++)
    {
        ImageXmipp I=getImage(i);
        I.rename(rootname+integerToString(i+1,5)+".xmp");
        I.write();
        SF.insert(I.name());
    }
    SF.write(rootname+".sel");
}

// Write as volume ---------------------------------------------------------
void ImageXmippStack::writeAsVolume(const FileName& name, bool force_reversed)
{
    VolumeXmipp V;
    V().resize(ImgNo(), YSIZE(getImage(0)()), XSIZE(getImage(0)()));
    for (int i=0; i<ImgNo(); i++)
    {
        V().setSlice(i,getImage(i)());
    }
    V.write(name,force_reversed);
}

// Initialise an oversampled image (ready for work) ------------------------
void ImageOver::init(int _vmin, int _vmax, int _vistep,
                     int _umin, int _umax, int _uistep)
{
    overvmin = _vmin;
    overumin = _umin;
    overvmax = _vmax;
    overumax = _umax;
    vistep = _vistep;
    uistep = _uistep;
    //   img.initZeros((_vmax-_vmin+1)*_vistep,(_umax-_umin+1)*_uistep);
    img.initZeros((_vmax - _vmin)*_vistep + 1, (_umax - _umin)*_uistep + 1);
    STARTINGY(img) = 0;
    STARTINGX(img) = 0;
    //   STARTINGY(img)=_vmin*_vistep - (_vistep-1)/2;
    //   STARTINGX(img)=_umin*_uistep - (_uistep-1)/2;
}

// Window ------------------------------------------------------------------
void ImageOver::window(int _v0, int _u0, int _vF, int _uF)
{
    overvmin = _v0;
    overumin = _u0;
    overvmax = _vF;
    overumax = _uF;

    int newYdim = (_vF - _v0) * vistep + 1;
    int newXdim = (_uF - _u0) * uistep + 1;
    img.setXmippOrigin();
    img.window(FIRST_XMIPP_INDEX(newYdim), FIRST_XMIPP_INDEX(newXdim),
               LAST_XMIPP_INDEX(newYdim), LAST_XMIPP_INDEX(newXdim));
    STARTINGY(img) = 0;
    STARTINGX(img) = 0;
}

// Clear -------------------------------------------------------------------
void ImageOver::clear()
{
    overvmin = overvmax = 0;
    overumin = overumax = 0;
    vistep = uistep = 0;
    Image::clear();
}

// Generate the normal image by averaging ----------------------------------
void ImageOver::downsample(ImageT< double > *I) const
{
    IMGMATRIX(*I).resize(overvmax - overvmin + 1, overumax - overumin + 1);
    for (int i = overvmin; i <= overvmax; i++)
        for (int j = overumin; j <= overumax; j++)
        {
            IMGPIXEL(*I, i, j) = 0;
            for (int v = (i - overvmin) * vistep; v < (i + 1 - overvmin)*vistep; v++)
                for (int u = (j - overumin) * uistep; u < (j + 1 - overumin)*uistep; u++)
                {
                    IMGPIXEL(*I, i, j) += IMGPIXEL(*this, u, v);
                }
            IMGPIXEL(*I, i, j) /= vistep * uistep;
        }
}

// Generate the oversample image by interpolation --------------------------
void ImageOver::oversample(ImageT < double > *I) const
    {}

/****************************************************************************/
/* IMAGIC IMAGES                                                            */
/****************************************************************************/

template <>
bool ImageImagicT<std::complex<double> >::read(const FileName &name)
{
    rename(name);
    ImageImagicinfo img_info = ImagicGetImgInfo(getHedFname());

    FileName img_fname = getImgFname();
    if (img_fname == "")
        REPORT_ERROR(1501, "ImageImagic::read: File " + name +
                     " doesn't seem fit Imagic format");
    FILE *img_fh;
    if ((img_fh = fopen(img_fname.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501, "ImageImagic::read: IMAGIC file " + img_fname +
                     " not found");

    const int imgnum = getImgNum();
    const size_t img_offset = img_info.xsize * img_info.ysize;
    // Read the image data
    const bool reversed = false;
    switch (img_info.img_types[imgnum])
    {
    case IMAGIC_COMP:
    {
        float a, b;
        const unsigned size = 4;
        fseek(img_fh, imgnum*(size*2)*img_offset, SEEK_SET);
    	std::complex<double>* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(img,n,ptr)
        {
            // read real part of a complex number
            FREAD(&a, size, 1, img_fh, reversed);
            // read imaginary part of a complex number
            FREAD(&b, size, 1, img_fh, reversed);
            // Assign the number
            std::complex<double> c(a, b);
            *ptr = c;
        }
        break;
    }
    default:
        REPORT_ERROR(1501, "ImageImagicType not supported for this imgtype!");
        break;
    }
    fclose(img_fh);
    return (true);
}

const ImageImagicinfo ImagicGetImgInfo(const FileName &hed_fname)
{
    ImageImagicinfo info;

    FILE *fp;
    if ((fp = fopen(hed_fname.c_str(), "rb")) != NULL)
    {
        // Read how many images (IFOL) and img size (NPIXEL)
        fseek(fp, IMAGIC_IFOL_OFFSET, SEEK_SET);
        unsigned int num_img;
        fread(&num_img, IMAGIC_WORD_LEN, 1, fp);
        num_img++; // Don't forget to include the current image
        info.num_img = num_img;
        fseek(fp, IMAGIC_IXLP_OFFSET, SEEK_SET);
        fread(&info.xsize, IMAGIC_WORD_LEN, 1, fp);
        fseek(fp, IMAGIC_IYLP_OFFSET, SEEK_SET);
        fread(&info.ysize, IMAGIC_WORD_LEN, 1, fp);
        for (unsigned int i = 0; i < num_img; i++)
        {
            // Determine what data type for this image
            char typeval[IMAGIC_WORD_LEN+1];
            fseek(fp, i*IMAGIC_RECORD_LEN + IMAGIC_TYPE_OFFSET, SEEK_SET);
            fread(&typeval, IMAGIC_WORD_LEN, 1, fp);
            typeval[IMAGIC_WORD_LEN] = '\0';
            // fIform values are defined in headerXmipp::set_header
            if (strncmp(typeval, "REAL", 4) == 0)
                info.img_types.push_back(IMAGIC_REAL);
            else if (strncmp(typeval, "INTG", 4) == 0)
                info.img_types.push_back(IMAGIC_INTG);
            else if (strncmp(typeval, "PACK", 4) == 0)
                info.img_types.push_back(IMAGIC_PACK);
            else if (strncmp(typeval, "COMP", 4) == 0)
                info.img_types.push_back(IMAGIC_COMP);
            else
                ; // throw an error
        }
        fclose(fp);
    }
    return (info);
}
