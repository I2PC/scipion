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

#ifndef __QT_IMAGE_CONVERTER_HH__
#define __QT_IMAGE_CONVERTER_HH__

#include <data/image.h>

#include <qimage.h>

/*  --------------------------------------------------- */
class QtImageConverter
{
public:
    // Constructor
    QtImageConverter();

    // qt to xmipp converter
    Image<double> *qt2xmipp(QImage *qtImg);

    // xmipp to qt converter
    QImage *xmipp2qt(Image<double> *xmippImg);
};

#endif
