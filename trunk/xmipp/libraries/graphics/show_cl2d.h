/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef SHOWCL2D_H
#define SHOWCL2D_H

#include "show_selfile.h"

/*@defgroup ShowCL2D Show CL2D images
   @ingroup GraphicsLibrary */
//@{
/** Class to show CL2D images. */
class ShowCL2D: public ShowSel
{
    Q_OBJECT;
public:
    // Apply geometric transformation when showing average & stddev of original images
    // and to compute the error image
    bool        apply_geo;
    // Filter suffix
    std::string filterSuffix;
protected:
    // Array of selfiles, a selfile for each codevector
    SelFile *SFcv;
    // Histogram of images assigned
    std::string  *hisAssigned;

    /* Initialization.
       Sets SFcv = NULL; and then calls to ShowSel::init() */
    virtual void init();
    /* Clear everything */
    virtual void clear();
    /* Read CL2D files and Sel file. */
    virtual void readFile(const FileName &_fn,
                          double _minGray = 0, double _maxGray = 0);

    /* Extract the represented images to a SelFile.
       The selfile is not cleared at the beginning. */
    virtual void extractRepresented(SelFile &SF_represented);

    /* Initialize right click menubar */
    virtual void initRightclickMenubar();
    /* Cell label depending on labeltype */
    virtual const char* cellLabel(int i) const;
    /* Some menu items for the options */
    int mi_hisAsLabels, mi_nameAsLabels;
#define Histogram_LABEL 2
#define Name_LABEL      3
    /* Change the labels */
    virtual void changeLabel(int _clicked_mi);

protected slots:
    /* Save assigned images in a new selfile */
    virtual void saveAssigned();
    /* Different kind of labels */
    void changeLabelsToName();
    void changeLabelsToHis();
    /* Particular shows */
    void showAssigned();
    void showThisImage();
    void showThisClass();
    void showThisPCA();
    void showThisOutliers();

public:
    /** Init with SOM files.
        The given parameter is the rootname and the following files are needed:
    */
    void initWithFile(const FileName &_fn_root);
};
//@}

#endif
