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

    int             splineDegree, dim, pyramid_level, fourier_threads;
    bool            applyTransform, inverse, wrap, isVol, flip, disableMetadata, temporaryOutput;
    Matrix2D<double> R, T, S, A, B;
    Matrix1D<double>          shiftV, rotV, scaleV;
    ImageGeneric img, imgOut;

    void defineParams();
    void readParams();
    void preProcess();
    void postProcess();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

};
#endif TRANSFORMGEOMETRY_H
