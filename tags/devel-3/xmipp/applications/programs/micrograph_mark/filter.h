/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 *
 *****************************************************************************/

#ifndef __QT_FILTER_HH__
#define __QT_FILTER_HH__

#include <data/micrograph.h>
#include <data/image.h>
#include <reconstruction/fourier_filter.h>

#include <qstring.h>

/* Filters generic class --------------------------------------------------- */
class QtFilter
{
public:
    static QString name;

protected:
    bool              __active;
    const Micrograph *__M;

public:
    // Constructor
    QtFilter(const Micrograph *_M);

    // Apply the filter
    virtual void apply(Image *_img) = 0;
};


/* Invert contrast --------------------------------------------------------- */
class QtInvertContrastFilter : public QtFilter
{
public:
    static QString name;
    QtInvertContrastFilter(const Micrograph *_M): QtFilter(_M)
    {}
    virtual void apply(Image *_img);
};

/* Enhance contrast -------------------------------------------------------- */
class QtEnhanceContrastFilter : public QtFilter
{
public:
    QtEnhanceContrastFilter(const Micrograph *_M): QtFilter(_M)
    {}
    static QString name;
    virtual void apply(Image *_img);
};

/* Substract background ---------------------------------------------------- */
class QtSubstractBackgroundFilter : public QtFilter
{
public:
    QtSubstractBackgroundFilter(const Micrograph *_M): QtFilter(_M)
    {}
    static QString name;
    virtual void apply(Image *_img);
};

/* Remove Outliers --------------------------------------------------------- */
class QtRemoveOutlierFilter : public QtFilter
{
public:
    QtRemoveOutlierFilter(const Micrograph *_M): QtFilter(_M)
    {}
    static QString name;
    virtual void apply(Image *_img);
};

/* Low Pass ---------------------------------------------------------------- */
class QtLowPassFilter : public QtFilter
{
public:
    QtLowPassFilter(const Micrograph *_M);
    static QString name;
    FourierMask filter;
    virtual void apply(Image *_img);
};

/* High Pass --------------------------------------------------------------- */
class QtHighPassFilter : public QtFilter
{
public:
    QtHighPassFilter(const Micrograph *_M);
    static QString name;
    FourierMask filter;
    virtual void apply(Image *_img);
};

#endif
