/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
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

#ifndef __QT_WIDGET_MICROGRAPH_HH__
#define __QT_WIDGET_MICROGRAPH_HH__

/* Includes ---------------------------------------------------------------- */
#include <qwidget.h>
#include <qpainter.h>
#include <qlayout.h>
#include <qmenubar.h>
#include <qaccel.h>
#include <qscrollbar.h>
#include <qlabel.h>
#include "QtImageMicrograph.hh"
#include "QtImageOverviewMicrograph.hh"
#include "QtFileMenu.hh"
#include <XmippData/xmippMasks.hh>
#include <vector>

/* Forward declarations ---------------------------------------------------- */ 
class QtMainWidgetMark;
class QtImageMicrograph;
class QtPopupMenuMark;
class QtFiltersController;
class Micrograph;

/* Particle ---------------------------------------------------------------- */
class Particle {
public:
   int x,y;          	   // position in micrograph
   int idx;                // Index of this particle within the micrograph
                           // list of coordinates
   char status;      	   // rejected=0, selected=1 or moved=2
   matrix1D<double> vec;   // vector of that particle
   double dist;      	   // distance from the avg vector

   // Print
   friend ostream & operator << (ostream &_out, const Particle &_p);
   
   // Read
   void read(istream &_in, int _vec_size);
};

/* Classification model ---------------------------------------------------- */
class Classification_model {
public:
   // Average of the model
   matrix1D<double>           __avg;
   // Covariance of the model
   matrix2D<double>           __sigma;
   // Inverse of sigma
   matrix2D<double>           __sigma_inv;
   // Example vectors
   vector< Particle >         __training_particle;
   // largest distance in the example set
   double                     __largest_distance;
   // Well posed
   bool                       __well_posed;
public:
   // Clear
   void clear();

   // Reserve space for a number of example particles
   void reserve_examples(int my_N) {__training_particle.reserve(my_N);}

   // Add example particle to model
   void add_particle(const Particle &p) {__training_particle.push_back(p);}

   // Import particles from another model
   void import_particles(const Classification_model &_model);

   // Build average and sigma
   void build_model();
   
   // Well posed. The model is well posed if the determinant of the covariance
   // matrix is different from zero
   bool well_posed() {
      __well_posed=(XSIZE(__sigma)>0)?(__sigma.det()>1e-1):false;
      return __well_posed;
   }
   
   // Compute the largest distance within the example set
   void compute_largest_distance();
   
   // Distance between two vectors
   double distance(const matrix1D<double> &my_X, const matrix1D<double> &my_Y);

   // Distance between a vector and the average
   double distance_to_average(const matrix1D<double> &my_X)
      {return distance(__avg,my_X);}
   
   // Euclidean distance between a vector and the average
   double euclidean_distance_to_average(const matrix1D<double> &my_X);
   
   // Print
   friend ostream & operator << (ostream &_out, const Classification_model &_m);
   
   // Read
   friend istream & operator >> (istream &_in, Classification_model &_m);
   
   // Print model
   void print_model(ostream &_out);
};

/* Widget for the micrograph ----------------------------------------------- */
class QtWidgetMicrograph : public QWidget {   
   Q_OBJECT
   
private:
   Micrograph                *__m;
   QtFiltersController       *__filtersController;
   int                        __activeFamily;   
   QMenuBar                  *__menuBar;      
   QtImageMicrograph         *__mImage;
   QtImageOverviewMicrograph *__mImageOverview;
   QVBoxLayout               *__gridLayout;
   QtFileMenu                *__file_menu;
   bool                       __tilted;
   int                        __mingray;
   int                        __maxgray;
   float                      __gamma;
   float                      __ellipse_radius;

   bool                       __learn_particles_done;
   bool                       __autoselection_done;
   Mask_Params                __mask;
   vector < matrix2D<int> * > __mask_classification;
   bool                       __use_background;
   vector < matrix1D<int> * > __radial_val;
   Classification_model       __training_model;
   Classification_model       __training_loaded_model;
   Classification_model       __auto_model;
   Classification_model       __auto_loaded_model;
   Classification_model       __selection_model;
   vector<Classification_model> __error_model;
   int                        __Nerror_models;
   bool                       __use_euclidean_distance_for_errors;
   int                        __auto_label;
   vector<int>                __error_index;
   matrix2D<double>           __piece;
   int                        __gray_bins;
   int                        __radial_bins;
   double                     __keep;
   double                     __highpass_cutoff;
   int                        __piece_xsize;
   int                        __piece_ysize;
   int                        __particle_radius;
   int                        __mask_size;
   int                        __min_distance_between_particles;
   int                        __output_scale;
   int                        __reduction; // Of the piece with respect
                                           // to the micrograph
   int                        __piece_overlap;
   int                        __particle_overlap;
   int                        __numin;
   int                        __numax;
   double                     __th1; // 0.6
   double                     __th2; // 0.9
   int                        __th3; // 50
   double                     __th4; // 0.8
   int                        __th5; // 10

public:
   // Constructor
   QtWidgetMicrograph( QtMainWidgetMark *_mainWidget, 
                       QtFiltersController *_f,
                       Micrograph *_m = NULL );
   ~QtWidgetMicrograph();
   
   // Set Micrograph
   void setMicrograph( Micrograph *_m );
   
   // Get Micrograph
   Micrograph *getMicrograph() { return( __m ); }
      
   // Set this as tilted micrograph
   void setTilted() {__tilted=TRUE; __mImage->setTilted();}

   // Is tilted?
   bool isTilted() {return __tilted;}

   // Get filters controller
   QtFiltersController *getFiltersController() { return(__filtersController); }
   
   // Get active family
   int activeFamily() { return( __activeFamily); }
   
   // Get overview
   QtImageOverviewMicrograph *overview() { return( __mImageOverview ); }
   
   // Get Image
   QtImageMicrograph *image() { return( __mImage ); }   
   
   // Get Filemenu
   QtFileMenu *file_menu() {return __file_menu;}
   
   // Add menu item
   void addMenuItem( const char *_msg, const QtPopupMenuMark *_item ) {
       __menuBar->insertItem( _msg, (QPopupMenu*)_item );
   }
   
   // Draw axis
   void draw_axis(double _ang)
      {__mImageOverview->enableAxis(); __mImageOverview->draw_axis(_ang);}
   
   // Open menu.
   // Add your menus to this function
   void openMenus();
   
   // Change contrast
   void changeContrast(int _mingray, int _maxgray, float _gamma);
   
   // Change contrast
   void changeCircleRadius(float _circle_radius);

   // Repaint 
   void repaint( int t=FALSE );
   
   // Learn particles
   void learnParticles();
   
   // Build Selection model. The selection model is made up of the training
   // and the automatically selected particles.
   void buildSelectionModel();

   // Create mask for learning particles
   void createMask();

   // Classify mask
   void classifyMask();

   // Build vectors
   void buildVectors(vector<int> &_idx, Classification_model &_model);
   
   // Build classfication vector
   // x,y are in the coordinate system of the piece (that might be
   // a reduced version of a piece in the micrograph)
   // (0,0) is the top-left corner
   // Returns true if the vector is successfully built
   bool build_vector(int _x, int _y, matrix1D<double> &_result);

   // Get a piece of the micrograph centered at position x,y (if possible)
   // the position of (x,y) in the piece is returned in (posx, posy)
   void get_centered_piece(int _x, int _y, int &_posx, int &_posy);

   // Get a piece whose top-left corner is at the desired position (if possible)
   // Returns true if the piece could be taken, and false if the whole
   // micrograph has been scanned. To scan the full micrograph,
   // Top and left should be initialized to 0,0 and this function should
   // succesive times. It returns the next top and left coordinates to
   // get the next piece along with the skips. The skips indicate which
   // part of the piece has been already scanned. This happens towards the
   // right and bottom boundaries of the micrograph, where the pieces
   // need to be shifted in order to fit with the required size.
   bool get_corner_piece(int _top, int _left, int _skip_y,
      int &_next_skip_x,int &_next_skip_y,int &_next_top,int &_next_left);
   
   // Denoise, reject outliers and equalize histogram
   // Returns true, if successful. False if unsuccessful (skip this piece)
   // Usually, it is unsuccessful if the denoising fails to work because
   // some "weird" features of the piece
   bool prepare_piece();

   //To get the neighbours of the particle at position (x,y) in the micrograph
   // (with actual coordinates in the piece posx,posy)
   // and their positions in the piece image
   void QtWidgetMicrograph::find_nbr(vector<int> &_idx, int _index,
      int _x, int _y,
      int _posx, int _posy, matrix1D<char> &_visited,
      vector< matrix1D<int> > &_nbr);

   // Automatically Select Particles
   void automaticallySelectParticles();

   // Given a current scanning position, this function returns
   // the next scanning position whithin the current piece.
   // The skips are given by get_corner_piece.
   // It returns true, if the next scanning position can be computed.
   // Otherwise, if the piece has been completely scanned, it returns false
   // Initialize _x,_y to 0,0 to scan the full piece (even if there are skips)
   bool get_next_scanning_pos(
      int &_x,int &_y, int _skip_x, int _skip_y);

   // Run over the list sorted by distances. If two particles are within
   // a given distance then either reject both or the one with largest distance
   // depending upon _reject_both. This function returns the number of particles
   // that are still candidates.
   int reject_within_distance(
      vector<Particle> &_Input, double _min_dist, bool _reject_both);

   // Reject those automatically selected particles that are very close
   // to manually selected ones.
   void reject_previously_selected(const Classification_model &_model,
      vector<Particle> &_candidate_vec);

   // Refine the position of a particle within the current piece
   void refine_center(Particle &my_P);

   // Add family.
   // The family label is returned
   int add_family(vector<Particle> &_list,
      const string &_family_name);

   // Move particle.
   // The input index is the index of the moved particle in the micrograph list
   void move_particle(int _idx);

   // Delete particle.
   // The input index is the index of the moved particle in the micrograph list
   void delete_particle(int _idx);
   
   // Rebuild the automatic vectors that have been moved
   void rebuild_moved_automatic_vectors();

   // load models
   void loadModels();
   
   // Save models
   void saveModels();
   
   // Classify errors
   void classify_errors();
   
   // Write all important information for particle selecting to a file
   void write() _THROW;
   
   // Configure auto
   void configure_auto();

public slots:
   void slotActiveFamily( int _f );
   void slotAddFamily( const char *_familyName );
   void slotDeleteMarkOther( int _coord );
   void slotChangeFamilyOther( int _coord, int _f );
   void slotRepaint() { repaint( FALSE ); }
   void slotDrawEllipse(int _x, int _y, int _f);
   void slotDrawLastEllipse(int _x, int _y, int _f);
   void slotQuit();
   void slotChangeContrast();
   void slotChangeCircleRadius();
signals:
   void signalActiveFamily( int _f );
   void signalAddFamily( const char *_familyName );
};

/** Class to adjust contrast
*/
class AdjustContrastWidget : public QWidget {
   Q_OBJECT
public:
   /** Constructor */
   AdjustContrastWidget(int min, int max, float gamma, 
      QtWidgetMicrograph *_qtwidgetmicrograph,
      QWidget *parent=0, const char *name=0, int wflags=0);
private:
   QtWidgetMicrograph *__qtwidgetmicrograph;
   QScrollBar 	      *__scroll_min;
   QScrollBar 	      *__scroll_max;
   QScrollBar 	      *__scroll_gamma;
   QLabel     	      *__label_min;
   QLabel     	      *__label_max;
   QLabel     	      *__label_gamma;
private slots:
   void scrollValueChanged(int);  
};

/** Class to adjust contrast
*/
class AdjustCircleRadiustWidget : public QWidget {
   Q_OBJECT
public:
   /** Constructor */
   AdjustCircleRadiustWidget(int min, int max, int start_with, 
      QtWidgetMicrograph *_qtwidgetmicrograph,
      QWidget *parent=0, const char *name=0, int wflags=0);
private:
   QtWidgetMicrograph *__qtwidgetmicrograph;
   QScrollBar 	      *__scroll_radius;
   QLabel     	      *__label_radius;   

private slots:
   void scrollValueChanged(int);  
};

#endif
