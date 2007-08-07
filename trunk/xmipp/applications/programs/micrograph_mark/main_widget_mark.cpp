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

#include "main_widget_mark.h"
#include "widget_micrograph.h"
#include "image_overview_micrograph.h"
#include "dialog_families.h"
#include "filters_controller.h"

#include <data/matrix2d.h>
#include <data/geometry.h>
#include <data/selfile.h>

#include <qaccel.h>

/* Constructor ------------------------------------------------------------- */
QtMainWidgetMark::QtMainWidgetMark(Micrograph *_m, Micrograph *_mTilted)
{
    int mXdim, mYdim;
    _m->size(mXdim, mYdim);
    float maspect_ratio = (float)mYdim / mXdim;
    if (_mTilted == NULL)
    {
        int suggested_X = 600;
        int suggested_Y = (int)(suggested_X * maspect_ratio);
        if (suggested_Y > 600)
        {
            suggested_Y = 600;
            suggested_X = (int)(suggested_Y / maspect_ratio);
        }
        setMinimumSize(suggested_X, suggested_Y);
        QRect G = geometry();
        G.setWidth(suggested_X);
        G.setHeight(suggested_Y);
        setGeometry(G);
    }
    else
    {
        int tXdim, tYdim;
        _m->size(tXdim, tYdim);
        float taspect_ratio = (float)tYdim / tXdim;
        int suggested_Y = 400;
        int suggested_X = (int)(suggested_Y * (1 / maspect_ratio + 1 / taspect_ratio));
        if (suggested_X > 800)
        {
            int suggested_YU = (int)(400 * maspect_ratio);
            int suggested_YT = (int)(400 * taspect_ratio);
            suggested_Y = MAX(suggested_YU, suggested_YT);
        }
        setMinimumSize(suggested_X, suggested_Y);
        QRect G = geometry();
        G.setWidth(suggested_X);
        G.setHeight(suggested_Y);
        setGeometry(G);
    }
    __tilted_generated = false;
    __untilted_generated = false;

    __gridLayout         = new QGridLayout(this, 1, 2, 0);
    __filtersController  = new QtFiltersController(this , _m);

    __mWidget            = new QtWidgetMicrograph(this, __filtersController);
    if (_mTilted == NULL) __mTiltedWidget = NULL;
    else
    {
        __mTiltedWidget = new QtWidgetMicrograph(this, __filtersController);
        __mTiltedWidget->setTilted();
    }

    __familyDialog       = new QtDialogFamilies(this);
    __familyDialog->setCaption("Families");
    __familyDialog->show();

    __ctrlPlus    = new QAccel(this);
    __ctrlMinus   = new QAccel(this);
    __otherCtrl   = new QAccel(this);

    connect(__familyDialog, SIGNAL(signalActiveFamily(int)),
            __mWidget, SLOT(slotActiveFamily(int)));
    connect(__mWidget, SIGNAL(signalAddFamily(const char*)),
            __familyDialog, SLOT(slotAddFamily(const char*)));

    __ctrlPlus->connectItem(__ctrlPlus->insertItem(Key_Plus + CTRL, 200),
                            __mWidget->image(), SLOT(slotZoomIn(void)));
    __ctrlMinus->connectItem(__ctrlMinus->insertItem(Key_Minus + CTRL, 201),
                             __mWidget->image(), SLOT(slotZoomOut(void)));
    __otherCtrl->connectItem(__otherCtrl->insertItem(Key_Q + CTRL, 200),
                             __mWidget, SLOT(slotQuit(void)));
    __otherCtrl->connectItem(__otherCtrl->insertItem(Key_S + CTRL, 201),
                             this, SLOT(slotSaveCoords(void)));

    if (_mTilted != NULL)
    {
        connect(__familyDialog, SIGNAL(signalActiveFamily(int)),
                __mTiltedWidget, SLOT(slotActiveFamily(int)));
        connect(__mTiltedWidget, SIGNAL(signalAddFamily(const char*)),
                __familyDialog, SLOT(slotAddFamily(const char*)));

        __ctrlPlus->connectItem(200, __mTiltedWidget->image(),
                                SLOT(slotZoomIn(void)));
        __ctrlMinus->connectItem(201, __mTiltedWidget->image(),
                                 SLOT(slotZoomOut(void)));

        connect((QObject*)__mWidget->overview(),
                SIGNAL(signalActualizeOtherOverview(int, int)),
                this,
                SLOT(slotActualizeTiltedOverview(int, int)));
        /* Disconnect the tilted->untilted connection
        connect( (QObject*)__mTiltedWidget->overview(),
                 SIGNAL(signalActualizeOtherOverview(int,int)),
                 (QObject*)__mWidget->overview(),
                 SLOT(slotActualizeOtherOverview(int,int)) );
        */

        connect((QObject*)__mTiltedWidget->image(),
                SIGNAL(signalRepaint(void)),
                (QObject*)__mWidget->overview(),
                SLOT(slotRepaint(void)));

        connect((QObject*)__mWidget->image(),
                SIGNAL(signalAddCoordOther(int, int, int)),
                this,
                SLOT(slotAddCoordTilted(int, int, int)));
        /* connect( (QObject*)__mTiltedWidget->image(),
                 SIGNAL(signalAddCoordOther(int,int,int)),
                 this,
                 SLOT(slotAddCoordUntilted(int,int,int)) );
        */

        connect((QObject*)__mWidget->image(),
                SIGNAL(signalDeleteMarkOther(int)),
                (QObject*)__mTiltedWidget,
                SLOT(slotDeleteMarkOther(int)));
        connect((QObject*)__mTiltedWidget->image(),
                SIGNAL(signalDeleteMarkOther(int)),
                (QObject*)__mWidget,
                SLOT(slotDeleteMarkOther(int)));

        connect((QObject*)__mWidget->image(),
                SIGNAL(signalChangeFamilyOther(int, int)),
                (QObject*)__mTiltedWidget,
                SLOT(slotChangeFamilyOther(int, int)));
        connect((QObject*)__mTiltedWidget->image(),
                SIGNAL(signalChangeFamilyOther(int, int)),
                (QObject*)__mWidget,
                SLOT(slotChangeFamilyOther(int, int)));

        connect((QObject*)__mWidget->image(),
                SIGNAL(signalRecalculateTiltMatrix()),
                this,
                SLOT(slotRecalculateTiltMatrix()));
        connect((QObject*)__mTiltedWidget->image(),
                SIGNAL(signalRecalculateTiltMatrix()),
                this,
                SLOT(slotRecalculateTiltMatrix()));
    }

    __mWidget->setMicrograph(_m);
    __familyDialog->setMicrograph(_m);
    if (_mTilted != NULL)
    {
        __mTiltedWidget->setMicrograph(_mTilted);
        __familyDialog->setTiltedMicrograph(_mTilted);
    }

    if (_mTilted == NULL)
    {
        __gridLayout->addMultiCellWidget(__mWidget, 0, 0, 0, 1);
    }
    else
    {
        __gridLayout->addWidget(__mWidget, 0, 0);
        __gridLayout->addWidget(__mTiltedWidget, 0, 1);
    }

    // Passing matrix initialization
    __Au.initZeros(3, 3);
    __Bt.initZeros(3, 3);
    __Nu = 0;
    __gamma = __alpha_u = __alpha_t = 0;
}

QtMainWidgetMark::~QtMainWidgetMark()
{
    delete __ctrlPlus;
    delete __ctrlMinus;
    delete __gridLayout;
    delete __familyDialog;
    delete __filtersController;
    delete __mWidget;
    if (__mTiltedWidget != NULL) delete __mTiltedWidget;
}

//#define _DEBUG
/* Add point to least squares matrices ------------------------------------- */
void QtMainWidgetMark::add_point(const Particle_coords &U,
                                 const Particle_coords &T)
{
    __Nu++; // Number of particles

#ifdef _DEBUG
    cout << "Adding point U(" << U.X << "," << U.Y << ") T(" << T.X << ","
    << T.Y << ")\n";
    cout << "A at input" << __Au << "B at input" << __Bt;
#endif
    // Adjust untilted dependent matrix
    __Au(0, 0) += U.X * U.X  ;
    __Au(0, 1) += U.X * U.Y  ;
    __Au(0, 2) += U.X;
    __Au(1, 0) = __Au(0, 1);
    __Au(1, 1) += U.Y * U.Y  ;
    __Au(1, 2) += U.Y;
    __Au(2, 0) = __Au(0, 2);
    __Au(2, 1) = __Au(1, 2);
    __Au(2, 2) = __Nu;

    // Adjust tilted dependent matrix
    __Bt(0, 0) += T.X * U.X  ;
    __Bt(0, 1) += T.Y * U.X  ;
    __Bt(0, 2) = __Au(0, 2);
    __Bt(1, 0) += T.X * U.Y  ;
    __Bt(1, 1) += T.Y * U.Y  ;
    __Bt(1, 2) = __Au(1, 2);
    __Bt(2, 0) += T.X      ;
    __Bt(2, 1) += T.Y      ;
    __Bt(2, 2) = __Au(2, 2);

#ifdef _DEBUG
    cout << "A at output" << __Au << "B at output" << __Bt;
#endif
}

/* Adjust passing matrix --------------------------------------------------- */
void QtMainWidgetMark::adjust_passing_matrix(const Particle_coords &U,
        const Particle_coords &T)
{
    add_point(U, T);
    if (can_use_passing_matrix())
    {
        solve(__Au, __Bt, __Put);
        __Put = __Put.transpose();
#ifdef _DEBUG
        cout << "Solved " << __Put;
#endif
        __Ptu = __Put.inv();
    }
}

#define Mu __mWidget->getMicrograph()
#define Mt __mTiltedWidget->getMicrograph()
/* Recalculate passing matrix ---------------------------------------------- */
void QtMainWidgetMark::recalculate_passing_matrix()
{
    __Au.initZeros();
    __Bt.initZeros();
    __Nu = 0;
    for (int i = 0; i < Mu->ParticleNo(); i++)
        add_point(Mu->coord(i), Mt->coord(i));
    if (can_use_passing_matrix())
    {
        solve(__Au, __Bt, __Put);
        __Put = __Put.transpose();
#ifdef _DEBUG
        cout << "Solved " << __Put;
#endif
        __Ptu = __Put.inv();
    }
}

/* Passing to tilted ------------------------------------------------------- */
void QtMainWidgetMark::pass_to_tilted(int _muX, int _muY,
                                      int &_mtX, int &_mtY, bool _update_passing_matrix)
{
    if (can_use_passing_matrix())
    {
        Matrix1D<double> m(3);
        SPEED_UP_temps;

        VECTOR_R3(m, _muX, _muY, 1);
#ifdef _DEBUG
        cout << "Input=" << m.transpose() << endl;
        cout << "Matrix=\n" << __Put;
#endif
        M3x3_BY_V3x1(m, __Put, m);
#ifdef _DEBUG
        cout << "Output=" << m.transpose() << endl;
#endif

        _mtX = (int)XX(m);
        _mtY = (int)YY(m);
    }
    else
    {
        _mtX = _muX;
        _mtY = _muY;
    }
    if (_update_passing_matrix)
    {
        Particle_coords T, U;
        U.X = _muX;
        U.Y = _muY;
        T.X = _mtX;
        T.Y = _mtY;
        adjust_passing_matrix(U, T);
    }
}
#undef DEBUG

/* Passing to tilted ------------------------------------------------------- */
void QtMainWidgetMark::pass_to_untilted(int _mtX, int _mtY, int &_muX,
                                        int &_muY)
{
    if (can_use_passing_matrix())
    {
        Matrix1D<double> m(3);
        SPEED_UP_temps;

        VECTOR_R3(m, _mtX, _mtY, 1);
        M3x3_BY_V3x1(m, __Ptu, m);
        _muX = (int)XX(m);
        _muY = (int)YY(m);
    }
    else
    {
        _muX = _mtX;
        _muY = _mtY;
    }
    Particle_coords T, U;
    U.X = _muX;
    U.Y = _muY;
    T.X = _mtX;
    T.Y = _mtY;
    adjust_passing_matrix(U, T);
}

/* Compute tilting angle --------------------------------------------------- */
void QtMainWidgetMark::compute_gamma()
{
#define TRIANGLE_NO 15000
#define MIN_AREA       15
#define MAX_AREA   250000
    __gamma = 0;
    // If there is no tilted image there is nothing to do
    if (__mTiltedWidget == NULL) return;

    int step = CEIL(pow((double)__Nu * __Nu * __Nu / TRIANGLE_NO, 1.0 / 3));
    Matrix1D<int> iju(2), iku(2), ijt(2), ikt(2); // From i to j in untilted
    // From i to k in untilted
    // From i to j in tilted
    // From i to k in tilted
    int triang = 0; // Number of triangles considered
    int i, j, k, counter1;
    counter1 = 0;
    randomize_random_generator();
    long noCombinations;
    noCombinations = __Nu * (__Nu - 1) * (__Nu - 2) / 6;
    while (triang < TRIANGLE_NO && counter1 < noCombinations)
    {
        counter1++;
        i = ROUND(rnd_unif(0, __Nu - 1));
        j = ROUND(rnd_unif(0, __Nu - 1));
        k = ROUND(rnd_unif(0, __Nu - 1));

        // Compute area of triangle in untilted micrograph
        VECTOR_R2(iju, Mu->coord(j).X - Mu->coord(i).X,
                  Mu->coord(j).Y - Mu->coord(i).Y);
        VECTOR_R2(iku, Mu->coord(k).X - Mu->coord(i).X,
                  Mu->coord(k).Y - Mu->coord(i).Y);
        double untilted_area = ABS(dotProduct(iju, iku)/*/2*/);
        if (untilted_area < MIN_AREA) continue; // For numerical stability

        // Compute area of the same triangle in the tilted micrograph
        VECTOR_R2(ijt, Mt->coord(j).X - Mt->coord(i).X,
                  Mu->coord(j).Y - Mu->coord(i).Y);
        VECTOR_R2(ikt, Mt->coord(k).X - Mt->coord(i).X,
                  Mt->coord(k).Y - Mt->coord(i).Y);
        double tilted_area = ABS(dotProduct(ijt, ikt)/*/2*/);
        if (tilted_area < MIN_AREA) continue; // For numerical stability
        if (tilted_area > MAX_AREA) continue; // micrograph are not perfect
        // sheets so avoid
        // very far away particles

        // Now we know that tilted_area=untilted_area*cos(gamma)
        if (tilted_area > untilted_area) continue; // There are user errors
        // In the point selection
        __gamma += acos(tilted_area / untilted_area);
        triang++;
    }
    __gamma /= triang;
    __gamma = RAD2DEG(__gamma);
    if (triang < 100)
        cout << "Not many particles, tilt angle may not be accurate" << endl;
}

/* Compute alphas ---------------------------------------------------------- */
Matrix2D<double> *pair_Put;
Matrix2D<double> pair_E;

double matrix_fitness(double *p)
{
    Euler_angles2matrix(p[1], p[3], -p[2], pair_E);
    double retval = 0;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            double error = ABS(pair_E(i, j) - (*pair_Put)(i, j));
            retval += error * error;
        }
    return retval;
}

void QtMainWidgetMark::compute_alphas()
{
    __alpha_u = __alpha_t = 0;
    // If there is no tilted image there is nothing to do
    if (__mTiltedWidget == NULL) return;

    Matrix1D<double> angles(3);
    angles.initZeros();
    double fitness;
    int iter;
    pair_Put = &__Put;

    // Coarse search
    double *aux = angles.adaptForNumericalRecipes();
    double best_alpha_u = 0, best_alpha_t = 0, best_fit = 1e8;
    aux[3] = __gamma;
    for (aux[1] = 0; aux[1] < 180; aux[1] += 10)
        for (aux[2] = 0; aux[2] < 180; aux[2] += 10)
        {
            double fit = matrix_fitness(aux);
            if (fit < best_fit)
            {
                best_fit = fit;
                best_alpha_u = aux[1];
                best_alpha_t = aux[2];
            }
        }
    angles.killAdaptationForNumericalRecipes(aux);
    angles(0) = best_alpha_u;
    angles(1) = best_alpha_t;
    angles(2) = __gamma;

    // Fine search
    Matrix1D<double> steps(3);
    steps.init_constant(1);
    powellOptimizer(angles, 1, 3, &matrix_fitness,
                     0.001, fitness, iter, steps, false);
    __alpha_u = angles(0);
    __alpha_t = angles(1);
    __gamma = angles(2);
}

/* Write angles in file ---------------------------------------------------- */
void QtMainWidgetMark::write_angles()
{
    if (!there_is_tilted()) return;
    recalculate_passing_matrix();
    compute_gamma();
    compute_alphas();
    draw_axes();

    FileName fn = __mWidget->getMicrograph()->micrograph_name();
    fn = fn.without_extension();
    fn = fn.add_extension("ang");

    ofstream out;
    out.open(fn.c_str(), ios::out);
    if (!out)
        REPORT_ERROR(1, (string)"QtMainWidgetMark::write_angles: Cannot open "
                     + fn + " for output\n");
    out << "# alpha_u alpha_t gamma\n"
    << __alpha_u << " " << __alpha_t << " " << __gamma << endl;
    out.close();
}

/* Draw tilt axes ---------------------------------------------------------- */
void QtMainWidgetMark::draw_axes()
{
    __mWidget->draw_axis(__alpha_u);
    __mTiltedWidget->draw_axis(__alpha_t);
}

/* Get informed that one of the micrographs generated images --------------- */
void QtMainWidgetMark::generated(bool _this_is_tilted,
                                 const string &_label)
{
    cerr << "Balancing both selfiles ...\n";
    if (_this_is_tilted)   __tilted_generated = true;
    else                 __untilted_generated = true;
    if (__untilted_generated && __tilted_generated)
    {
        FileName fn_untilted =
            __mWidget->getMicrograph()->micrograph_name() + "." + _label + ".sel";
        FileName fn_tilted =
            __mTiltedWidget->getMicrograph()->micrograph_name() + "." + _label + ".sel";
        SelFile SFUntilted(fn_untilted);
        SelFile SFTilted(fn_tilted);
        SFUntilted.go_beginning();
        SFTilted.go_beginning();
        while (!SFUntilted.eof() && !SFTilted.eof())
        {
            if (SFUntilted.Is_DISCARDED() || SFTilted.Is_DISCARDED())
            {
                SFUntilted.set_current(SelLine::DISCARDED);
                SFTilted.set_current(SelLine::DISCARDED);
                cerr << "Images " << SFUntilted.get_current_file() << " and "
                << SFTilted.get_current_file() << " are discarded\n";
            }
            SFUntilted.next();
            SFTilted.next();
        }
        SFUntilted.write(fn_untilted);
        SFTilted.write(fn_tilted);
    }
}

void QtMainWidgetMark::slotAddCoordTilted(int _muX, int _muY, int _f)
{
    int mtX, mtY;

    pass_to_tilted(_muX, _muY, mtX, mtY, true);
    cout << "     --> in the tilted image (" << mtX << "," << mtY << ")\n";
    __mTiltedWidget->getMicrograph()->add_coord(mtX, mtY, _f);

    __mTiltedWidget->slotDrawEllipse(mtX, mtY, _f);
    slotActualizeTiltedOverview(_muX, _muY);
// CO:   __mWidget->slotDrawEllipse( _muX, _muY, _f);
    __mTiltedWidget->slotDrawLastEllipse(mtX, mtY, _f);
}

void QtMainWidgetMark::slotAddCoordUntilted(int _mtX, int _mtY, int _f)
{
    int muX, muY;

    pass_to_untilted(_mtX, _mtY, muX, muY);
    __mWidget->getMicrograph()->add_coord(muX, muY, _f);

// CO:   __mTiltedWidget->slotDrawEllipse(_mtX, _mtY, _f);
    __mWidget->slotDrawEllipse(muX, muY, _f);
}

void QtMainWidgetMark::slotActualizeTiltedOverview(int _muX, int _muY)
{
    int mtX, mtY;
    pass_to_tilted(_muX, _muY, mtX, mtY, false);
    int mMaxX, mMaxY;
    __mTiltedWidget->getMicrograph()->size(mMaxX, mMaxY);
    mtX = CLIP(mtX, 0, mMaxX - 1);
    mtY = CLIP(mtY, 0, mMaxY - 1);
    __mTiltedWidget->overview()->slotActualizeOtherOverview(mtX, mtY);
}

void QtMainWidgetMark::slotSaveCoords()
{
    __mWidget->file_menu()->slotSaveCoords();
    if (__mTiltedWidget != NULL)
    {
        __mTiltedWidget->file_menu()->slotSaveCoords();
        __mWidget->file_menu()->slotSaveAngles();
    }
}
