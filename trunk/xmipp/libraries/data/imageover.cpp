#include "imageover.h"


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
    //   data.initZeros((_vmax-_vmin+1)*_vistep,(_umax-_umin+1)*_uistep);
    data.initZeros((_vmax - _vmin)*_vistep + 1, (_umax - _umin)*_uistep + 1);
    STARTINGY(data) = 0;
    STARTINGX(data) = 0;
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
    data.setXmippOrigin();
    data.window(FIRST_XMIPP_INDEX(newYdim), FIRST_XMIPP_INDEX(newXdim),
                LAST_XMIPP_INDEX(newYdim), LAST_XMIPP_INDEX(newXdim));
    STARTINGY(data) = 0;
    STARTINGX(data) = 0;
}

// Clear -------------------------------------------------------------------
void ImageOver::clear()
{
    overvmin = overvmax = 0;
    overumin = overumax = 0;
    vistep = uistep = 0;
    Image<double>::clear();
}

// Generate the normal image by averaging ----------------------------------
void ImageOver::downsample(Image< double > *I) const
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
void ImageOver::oversample(Image < double > *I) const
    {}

