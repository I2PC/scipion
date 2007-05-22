/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#ifndef __QT_COLORLABEL_HH__
#define __QT_COLORLABEL_HH__

/* Includes ---------------------------------------------------------------- */

/* Forward declarations ---------------------------------------------------- */
class QColor;

/* Set of colors for the labels -------------------------------------------- */
class QtColorLabel
{

private:
    static const int  __numColors = 15;
    QColor           *__col;

public:
    // Constructor
    QtColorLabel();
    ~QtColorLabel();

    // Get the color for the i family
    QColor col(int i);
};

#endif
