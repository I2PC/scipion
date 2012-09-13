/***************************************************************************
 *
 * Authors:      Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef SHOWASSIGNCTF_H
#define SHOWASSIGNCTF_H

#include "show_2d.h"
#include "show_tools.h"

#include <reconstruction/ctf_estimate_from_micrograph.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QPaintEvent>
#endif

/**@defgroup AssignCTFViewer Assign CTF Viewer
   @ingroup GraphicsLibrary */
//@{

/** Assign CTF Viewer.
    This class shows a psd and allows to select some parameters for
    adjusting the CTF.
*/
class AssignCTFViewer: public ImageViewer
{
    Q_OBJECT
public:
    /// Constructor with a PSD image
    AssignCTFViewer(const FileName &_fn_psd, Prog_assign_CTF_prm &_assign_ctf_prm);

    /// Draw first zero
    void drawFirstZero(std::vector<float> &prm);

    /// Update mask
    void updateMask(std::vector<float> &prm);
public:
    // Assign CTF parameters
    Prog_assign_CTF_prm assign_ctf_prm;

    // PSD filename
    FileName fn_psd;

    // Scroll parameters
    ScrollParam *select_prm;

    // Current vector of parameters
    // 0: Min freq
    // 1: Max freq
    // 2: Defocus U
    // 3: Defocus V
    // 4: Angle
    std::vector<float> current_prm;

    // Backup of the original image
    MultidimArray<double> xmippImage_backup;
public slots:
    // Set a new set of parameters
    void set_prm(std::vector<float> new_prm);
    // cancel
    void cancel();
    // proceed
    void okToProceed();
protected slots:
    // For repainting
    virtual void paintEvent(QPaintEvent *);
};
//@}

#endif
