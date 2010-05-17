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

#ifndef __QT_FILTERS_CONTROLLER_HH__
#define __QT_FILTERS_CONTROLLER_HH__

#include <vector>
#include <qwidget.h>

#include <data/micrograph.h>

/* Forward declarations ---------------------------------------------------- */
class QImage;
class QtFilter;
class QDialog;
#ifdef QT3_SUPPORT
class Q3ListBox;
#else
class QListBox;
#endif

/* Filters generic class --------------------------------------------------- */
class QtFiltersController : public QWidget
{
    Q_OBJECT

private:
    std::vector<QtFilter*>  __filterList;
    QDialog           *__addFilterDialog;
#ifdef QT3_SUPPORT
    Q3ListBox          *__listFilters;
#else
    QListBox          *__listFilters;
#endif
    const Micrograph  *__M;

    enum filters
    {
        invertContrastFilter,
        enhanceContrastFilter,
        substractBackgroundFilter,
        removeOutlierFilter,
        lowpassFilter,
        highpassFilter
    };

public:
    // Constructor
    QtFiltersController(QWidget * _parent, const Micrograph *_M);
    ~QtFiltersController();

    // Apply the filters list
    void applyFilters(QImage *_img);

public slots:
    void slotAddFilter();
    void slotAddFilter(int _f);
    void slotCleanFilters();
};

#endif
