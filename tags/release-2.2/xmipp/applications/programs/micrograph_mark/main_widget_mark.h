/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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

#ifndef __QT_MAIN_WIDGET_MARK_HH__
#define __QT_MAIN_WIDGET_MARK_HH__

#include <qwidget.h>
#include <qlayout.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <Q3GridLayout>
#endif

#include <data/micrograph.h>

/* Forward declarations ---------------------------------------------------- */
class QtWidgetMicrograph;
class QtDialogFamilies;
class QtFiltersController;
#ifdef QT3_SUPPORT
class Q3Accel;
#else
class QAccel;
#endif

/* Main mark widget -------------------------------------------------------- */
class QtMainWidgetMark : public QWidget
{
    // For accepting signals and slots
    Q_OBJECT

private:
    QtWidgetMicrograph  *__mWidget;
    QtWidgetMicrograph  *__mTiltedWidget;
    QtDialogFamilies    *__familyDialog;
#ifdef QT3_SUPPORT
    Q3GridLayout         *__gridLayout;
    Q3Accel              *__ctrlPlus;
    Q3Accel              *__ctrlMinus;
    Q3Accel              *__otherCtrl;
#else
    QGridLayout         *__gridLayout;
    QAccel              *__ctrlPlus;
    QAccel              *__ctrlMinus;
    QAccel              *__otherCtrl;
#endif
    QtFiltersController *__filtersController;

    // For tilted-untilted correspondance
    Matrix2D<double>    __Au;     // Untilted "positions"
    Matrix2D<double>    __Bt;     // Tilted   "positions"
    Matrix2D<double>    __Put;    // Passing matrix from untilted to tilted
    Matrix2D<double>    __Ptu;    // Passing matrix from tilted to untilted
    int                 __Nu;     // Number of points in matrices
    double              __gamma;  // Tilting angle in radians
    double              __alpha_u;// Angle of axis with X axis in radians
    double              __alpha_t;// Anfle of axis with X axis in radians

    // To balance the selfiles of a pair
    bool                __untilted_generated;
    bool                __tilted_generated;


public:
    // Constructor
    QtMainWidgetMark(Micrograph *_m, Micrograph *_mTilted = NULL);
    ~QtMainWidgetMark();

    // Access to WidgetMicrographs
    const QtWidgetMicrograph * untilted_widget()
    {
        return __mWidget;
    }
    const QtWidgetMicrograph * tilted_widget()
    {
        return __mTiltedWidget;
    }

    // Add point to least squares matrices
    void add_point(const Particle_coords &U, const Particle_coords &T);

    // Add point to least squares matrices and solve for the Passing matrix
    void adjust_passing_matrix(const Particle_coords &U,
                               const Particle_coords &T);

    /* Recalculate Passing matrix.
       This is only needed when a particle is moved */
    void recalculate_passing_matrix();

    // Pass from untilted to tilted
    void pass_to_tilted(int _muX, int _muY, int &_mtX, int &_mtY,
                        bool _update_passing_matrix);

    // Pass from tilted to untilted
    void pass_to_untilted(int _mtX, int _mtY, int &_muX, int &_muY);

    // Can use passing matrix?
    bool can_use_passing_matrix()
    {
        return __Nu > 3;
    }

    // Compute gamma
    void compute_gamma();

    // Compute alphas
    void compute_alphas();

    // Write angles in file ".ang"
    void write_angles();

    // Draw tilt axes in both micrographs
    void draw_axes();

    // Return alpha for untilted
    double alpha_u() const
    {
        return __alpha_u;
    }

    // Return alpha for tilted
    double alpha_t() const
    {
        return __alpha_t;
    }

    // Return tilt for tilted
    double gamma_t() const
    {
        return __gamma;
    }

    // True if there is a tilted micrograph
    bool there_is_tilted() const
    {
        return __mTiltedWidget != NULL;
    }

    // The main widget is informed when each of the micrographs
    // has generated its set of values. When both micrographs
    // generated the images, the selfiles are compared to remove
    // those images that were discarded in only one of them
    void generated(bool _this_is_tilted, const std::string &_label);
public slots:
    void slotAddCoordTilted(int _mX, int _mY, int _f);
    void slotAddCoordUntilted(int _mX, int _mY, int _f);
    void slotRecalculateTiltMatrix()
    {
        recalculate_passing_matrix();
    }
    void slotActualizeTiltedOverview(int _mX, int _mY);
    void slotSaveCoords();
};

#endif
