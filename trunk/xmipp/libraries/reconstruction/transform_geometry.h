#ifndef TRANSFORMGEOMETRY_H
#define TRANSFORMGEOMETRY_H

#include <data/transformations.h>
#include <data/xmipp_image_generic.h>
#include <data/metadata_extension.h>
#include <data/xmipp_fftw.h>
#include <data/xmipp_program.h>
#include <data/matrix2d.h>


typedef enum { SCALE_NONE, SCALE_FACTOR, SCALE_DIM, SCALE_FOURIER, SCALE_PYRAMID_EXPAND, SCALE_PYRAMID_REDUCE } ScaleType;
#define INTERP_FOURIER -1

class ProgTransformGeometry: public XmippMetadataProgram
{
public:
  /** Constructor and destructor, just to avoid vtable undefined references errors */
  ProgTransformGeometry();
  ~ProgTransformGeometry();

protected:
    ScaleType scale_type;

    int             zdim, ydim, xdim, splineDegree, dim, pyramid_level, fourier_threads;
    size_t          ndim;
    bool            applyTransform, inverse, wrap, isVol, flip;
    Matrix2D<double> R, T, S, A, B;
    Matrix1D<double>          shiftV, rotV, scaleV;
    MDRow            input, transformation;
    ImageGeneric img, imgOut;
    bool          output_is_stack;

    void defineParams();
    void readParams();
    void preProcess();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId);

};
#endif TRANSFORMGEOMETRY_H
