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

#ifndef __QT_MAIN_WIDGET_MARK_HH__
#define __QT_MAIN_WIDGET_MARK_HH__

/* Includes ---------------------------------------------------------------- */
#include "qwidget.h"
#include "qlayout.h"
#include "XmippData/xmippMicrograph.hh"

/* Forward declarations ---------------------------------------------------- */ 
class QtWidgetMicrograph;
class QtDialogFamilies;
class QtFiltersController;
class QAccel;

/* Main mark widget -------------------------------------------------------- */
class QtMainWidgetMark : public QWidget {
   // For accepting signals and slots
   Q_OBJECT   
   
private:
   QtWidgetMicrograph  *__mWidget;
   QtWidgetMicrograph  *__mTiltedWidget;
   QtDialogFamilies    *__familyDialog;
   QGridLayout         *__gridLayout;
   QAccel              *__ctrlPlus;
   QAccel              *__ctrlMinus;
   QtFiltersController *__filtersController;
   
   // For tilted-untilted correspondance
   matrix2D<double>    __Au;     // Untilted "positions"
   matrix2D<double>    __Bt;     // Tilted   "positions"
   matrix2D<double>    __Put;    // Passing matrix from untilted to tilted
   matrix2D<double>    __Ptu;    // Passing matrix from tilted to untilted
   int                 __Nu;     // Number of points in matrices
   double              __gamma;  // Tilting angle in radians
   double              __alpha_u;// Angle of axis with X axis in radians
   double              __alpha_t;// Anfle of axis with X axis in radians

public:
   // Constructor
   QtMainWidgetMark( Micrograph *_m, Micrograph *_mTilted=NULL );
   ~QtMainWidgetMark();

   // Access to WidgetMicrographs
   const QtWidgetMicrograph * untilted_widget() {return __mWidget;}
   const QtWidgetMicrograph * tilted_widget() {return __mTiltedWidget;}

   // Add point to least squares matrices
   void add_point(const Particle_coords &U, const Particle_coords &T);
   
   // Add point to least squares matrices and solve for the Passing matrix
   void adjust_passing_matrix(const Particle_coords &U,
      const Particle_coords &T);

   /* Recalculate Passing matrix.
      This is only needed when a particle is moved */
   void recalculate_passing_matrix();
   
   // Pass from untilted to tilted
   void pass_to_tilted(int _muX, int _muY, int &_mtX, int &_mtY);
   
   // Pass from tilted to untilted
   void pass_to_untilted(int _mtX, int _mtY, int &_muX, int &_muY);
   
   // Can use passing matrix?
   bool can_use_passing_matrix() {return __Nu>3;}

   // Compute gamma
   void compute_gamma();

   // Compute alphas
   void compute_alphas();
   
   // Write angles in file ".ang"
   void write_angles() _THROW;
   
   // Draw tilt axes in both micrographs
   void draw_axes();

   // Return alpha for untilted
   double alpha_u() const {return __alpha_u;}

   // Return alpha for tilted
   double alpha_t() const {return __alpha_t;}

   // True if there is a tilted micrograph
   bool there_is_tilted() const {return __mTiltedWidget!=NULL;}
   
public slots:
   void slotAddCoordTilted( int _mX, int _mY, int _f );
   void slotAddCoordUntilted( int _mX, int _mY, int _f );
   void slotRecalculateTiltMatrix() { recalculate_passing_matrix(); }
   void slotActualizeTiltedOverview( int _mX, int _mY);
};

#endif
