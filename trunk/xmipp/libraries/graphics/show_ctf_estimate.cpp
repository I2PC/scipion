/***************************************************************************
 *
 * Authors:      Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "show_ctf_estimate.h"

#include <qpainter.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QPaintEvent>
#endif

// Constructor -------------------------------------------------------------
AssignCTFViewer::AssignCTFViewer(const FileName &_fn_psd,
                                 Prog_assign_CTF_prm &_assign_ctf_prm):
        ImageViewer(_fn_psd.c_str(), false)
{
    // Get input parameters .................................................
    assign_ctf_prm = _assign_ctf_prm;
    fn_psd = _fn_psd;
    if (assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.K == 0)
        assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.K = 1;

    // Open a window for the scroll parameters ..............................
    std::vector<float> min, max, initial_value;
    std::vector<char *> prm_name;

    // Set min_freq
    min.push_back(0);
    max.push_back(50);
    prm_name.push_back("Minimum freq. (x100)");
    initial_value.push_back(assign_ctf_prm.adjust_CTF_prm.min_freq*100);
    // Set max_freq
    min.push_back(0);
    max.push_back(50);
    prm_name.push_back("Maximum freq. (x100)");
    initial_value.push_back(assign_ctf_prm.adjust_CTF_prm.max_freq*100);
    // Set DefocusU
    min.push_back(-100);
    max.push_back(-1);
    prm_name.push_back("DefocusU (/1000)");
    initial_value.push_back(
        assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.DeltafU / 1000);
    // Set DefocusV
    min.push_back(-100);
    max.push_back(-1);
    prm_name.push_back("DefocusV (/1000)");
    initial_value.push_back(
        assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.DeltafV / 1000);
    // Set Angle
    min.push_back(0);
    max.push_back(360);
    prm_name.push_back("Angle");
    initial_value.push_back(
        assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.azimuthal_angle);
    current_prm = initial_value;

    // Open window for scrolls
    select_prm = new ScrollParam(min, max, initial_value, prm_name,
                                 "Recompute CTF model", 0, "new window", Qt::WDestructiveClose, 0);
    connect(select_prm, SIGNAL(new_value(std::vector<float>)),
            this,       SLOT(set_prm(std::vector<float>)));
    connect(select_prm, SIGNAL(signal_ok_clicked()),
            this,       SLOT(okToProceed()));
    connect(select_prm, SIGNAL(signal_close_clicked()),
            this,       SLOT(cancel()));
    select_prm->setFixedSize(200, 175);
    select_prm->show();

    // Show the PSD .........................................................
    apply_geo = false;
    loadImage(fn_psd.c_str(), 0, 0, ImageViewer::PSD_mode);
    xmippImage_backup = xmippImage();
    updateMask(current_prm);
    show();
}

// Set parameters ----------------------------------------------------------
void AssignCTFViewer::set_prm(std::vector<float> new_prm)
{
    // If there is a change of min, max freq
    if (current_prm[0] != new_prm[0] || current_prm[1] != new_prm[1])
        updateMask(new_prm);

    // If there is a change of defocus or angle
    if (current_prm[2] != new_prm[2] || current_prm[3] != new_prm[3] ||
        current_prm[4] != new_prm[4])
    {
        drawFirstZero(current_prm); // To remove current ellipse
        drawFirstZero(new_prm);     // To draw a new one
    }
    current_prm = new_prm;
}

// Update mask -------------------------------------------------------------
void AssignCTFViewer::updateMask(std::vector<float> &prm)
{
    xmippImage() = xmippImage_backup;
    double min_freq = prm[0] / 100;
    double max_freq = prm[1] / 100;
    double r2_min = min_freq * XSIZE(xmippImage_backup);
    r2_min *= r2_min;
    double r2_max = max_freq * XSIZE(xmippImage_backup);
    r2_max *= r2_max;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(xmippImage_backup)
    {
        double r2 = i * i + j * j;
        if (r2 < r2_min || r2 > r2_max) xmippImage(i, j) = 0;
    }
    xmipp2Qt(xmippImage);
    showImage();
    repaint();
    updateStatus();
}

// Draw first zero ---------------------------------------------------------
void AssignCTFViewer::drawFirstZero(std::vector<float> &prm)
{
    // Setup a CTF model with the known parameters
    XmippCTF ctfmodel;
    ctfmodel = assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel;
    ctfmodel.DeltafU        = prm[2] * 1000;
    ctfmodel.DeltafV        = prm[3] * 1000;
    ctfmodel.azimuthal_angle = prm[4];
    ctfmodel.Produce_Side_Info();

    double fy = (double)(height() - status->height()) / (double)image.height();
    double fx = (double)(width()) / (double)image.width();

    QPainter painter(this);
    QBrush brush(Qt::NoBrush);
    QPen myPen(Qt::red, 3);
    painter.setPen(myPen);
    painter.setBrush(brush);
    /*
     * painter.setRasterOp(XorROP);
     */

    Matrix1D<int> previous_point(2);
    for (double ang = 0; ang <= 360; ang += 1)
    {
        Matrix1D<double> freq(2), dir(2);
        Matrix1D<int>    current_point(2);

        // Compute first zero in the given direction
        VECTOR_R2(dir, COSD(ang), SIND(ang));
        ctfmodel.zero(1, dir, freq);
        contfreq2digfreq(freq, freq, ctfmodel.Tm);
        XX(current_point) = (int)(XX(freq) * XSIZE(xmippImage_backup))
                            - STARTINGX(xmippImage_backup);
        YY(current_point) = (int)(YY(freq) * YSIZE(xmippImage_backup))
                            - STARTINGY(xmippImage_backup);
        if (ang != 0)
            painter.drawLine(
                ROUND(fx*XX(current_point)),  ROUND(fy*YY(current_point)),
                ROUND(fx*XX(previous_point)), ROUND(fy*YY(previous_point)));
        previous_point = current_point;
    }
}

// Cancel ------------------------------------------------------------------
void AssignCTFViewer::cancel()
{
    close();
}

// Ok to Proceed -----------------------------------------------------------
void AssignCTFViewer::okToProceed()
{
    // Create adjust parameters file
    FileName fn_random;
    fn_random.init_random(15);
    fn_random = (std::string)"PPP" + fn_random + ".txt";
    std::ofstream fh_adjust_param;
    fh_adjust_param.open(fn_random.c_str());
    if (!fh_adjust_param)
        REPORT_ERROR(1, "ShowSel::recomputeCTFmodel: Cannot open "
                     "file for output");

    // Write adjust_CTF parameters
    fh_adjust_param
    << "ctf=         " << fn_psd             << std::endl
    << "min_freq=    " << current_prm[0] / 100 << std::endl
    << "max_freq=    " << current_prm[1] / 100 << std::endl
    ;
    if (!assign_ctf_prm.adjust_CTF_prm.modelSimplification)
        fh_adjust_param << "radial_noise=yes\n";
    fh_adjust_param << "defocus_range=4000\n";
    fh_adjust_param << "show_optimization=yes\n";
    fh_adjust_param << std::endl;

    // Write CTF parameters
    fh_adjust_param
    << "voltage=             " << assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.kV << std::endl
    << "spherical_aberration=" << assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.Cs << std::endl
    << "sampling_rate=       " << assign_ctf_prm.adjust_CTF_prm.initial_ctfmodel.Tm << std::endl
    << "defocusU=            " << current_prm[2]*1000              << std::endl
    << "defocusV=            " << current_prm[3]*1000              << std::endl
    << "azimuthal_angle=     " << current_prm[4]                   << std::endl
    << "ctfmodelSize=        " << assign_ctf_prm.adjust_CTF_prm.ctfmodelSize      << std::endl
    ;

    fh_adjust_param.close();

    // Recompute the model
    system(((std::string)"( xmipp_ctf_estimate_from_psd -i " + fn_random + " ; rm " + fn_random + " ) &").c_str());

    // Close this window
    close();
}

// Paint event ------------------------------------------------------------
void AssignCTFViewer::paintEvent(QPaintEvent *e)
{
    ImageViewer::paintEvent(e);
    drawFirstZero(current_prm);
}
