/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
 *              Enrique Recarte Llorens (erecallo@hotmail.com)
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

#include "widget_micrograph.h"
#include "filter_menu.h"
#include "auto_menu.h"
#include "image_micrograph.h"
#include "image_overview_micrograph.h"

#include <data/micrograph.h>
#include <data/args.h>
#include <data/filters.h>
#include <data/rotational_spectrum.h>
#include <data/denoise.h>
#include <data/fft.h>
#include <reconstruction/fourier_filter.h>
#include <graphics/show_tools.h>

#include <algorithm>
#include <sys/time.h>

#include <qinputdialog.h>
#include <qmessagebox.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qpushbutton.h>

#ifdef QT3_SUPPORT
#include <q3grid.h>
//Added by qt3to4:
#include <Q3GridLayout>
#include <Q3VBoxLayout>
#else
#include <qgrid.h>
#endif

//#define DEBUG_AUTO
//#define DEBUG_MORE_AUTO
//#define DEBUG_PREPARE
//#define DEBUG_CLASSIFY
//#define DEBUG_BUILDVECTOR
//#define DEBUG_IMG_BUILDVECTOR

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const Particle &_p)
{
    _out << _p.x      << " " << _p.y << " "
         << _p.idx    << " "
         << (int)_p.status << " "
         << _p.cost   << " "
         << _p.vec.transpose()
         << std::endl
    ;
    return _out;
}

/* Read--------------------------------------------------------------------- */
void Particle::read(std::istream &_in, int _vec_size)
{
    _in >> x >> y
        >> idx
        >> status
        >> cost;
    status -= '0';
    vec.resize(_vec_size);
    _in >> vec;
}

/* Clear ------------------------------------------------------------------- */
void Classification_model::clear()
{
    __particles_picked.clear();
    __micrographs_number = 0;
    __falsePositives.clear();
    
    int imax = __training_particles.size();
    for (int i = 0; i < imax; i++)
        __training_particles[i].clear();
}

/* Import model ------------------------------------------------------------ */
void Classification_model::import_model(const Classification_model &_model)
{
    // Import particles
    int jmax = _model.__classNo;
    for (int j = 0; j < jmax; j++)
    {
        int imax = _model.__training_particles[j].size();
	    for (int i = 0; i < imax; i++)
            addParticleTraining(_model.__training_particles[j][i], j);
    }

    // Import parameters
    int imax = _model.__micrographs_number;
    for (int i = 0; i < imax; i++)
    {
        addParticlePicked(_model.__particles_picked[i]);    
        addMicrographScanned(_model.__micrographs_scanned[i]);
        addFalsePositives(_model.__falsePositives[i]);
    }

    // Update the number of micrographs
    __micrographs_number += _model.__micrographs_number;
}

/* Initialize -------------------------------------------------------------- */
void Classification_model::initNaiveBayes(
    const std::vector < Matrix2D<double> > &features,
    const Matrix1D<double> &probs, int discreteLevels,
    double penalization)
{
    __bayesNet=new NaiveBayes(features, probs, discreteLevels);
    int K=features.size();
    Matrix2D<double> cost(K,K);
    cost.initConstant(1);
    for (int i=0; i<XSIZE(cost); i++)
       cost(i,i)=0;
    cost(0,K-1)=penalization;
    __bayesNet->setCostMatrix(cost);
    __bayesClassifier=0;
}

void Classification_model::initNaiveBayesEnsemble(
    const std::vector < Matrix2D<double> > &features,
    const Matrix1D<double> &probs, int discreteLevels,
    double penalization, int numberOfClassifiers,
    double samplingFeatures, double samplingIndividuals,
    const std::string &newJudgeCombination)
{
    __bayesEnsembleNet=new EnsembleNaiveBayes(features, probs, discreteLevels,
        numberOfClassifiers, samplingFeatures, samplingIndividuals,
        newJudgeCombination);
    int K=features.size();
    Matrix2D<double> cost(K,K);
    cost.initConstant(1);
    for (int i=0; i<XSIZE(cost); i++)
       cost(i,i)=0;
    cost(0,K-1)=penalization;
    __bayesEnsembleNet->setCostMatrix(cost);
    __bayesClassifier=1;
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const Classification_model &_m)
{
    _out << "#Already_processed_parameters...\n";
    _out << "#Micrographs_processed= " << _m.__micrographs_number << std::endl;
    _out << "#Micrographs_processed_points= ";
    
    for (int i = 0; i < _m.__micrographs_number; i++)
        _out << _m.__micrographs_scanned[i] << " ";
    _out << "\n";

    _out << "#Particles_picked_per_micrograph= ";
    for (int i = 0; i < _m.__micrographs_number; i++)
        _out << _m.__particles_picked[i] << " ";
    _out << "\n";
    
    _out << "#FalsePos_Picked= ";
    for (int i = 0; i < _m.__micrographs_number; i++)
        _out << _m.__falsePositives[i] << " ";
    _out << "\n";
       
    _out << "#Model_parameters..." << std::endl;
    _out << "#Vector_size= " << XSIZE(_m.__training_particles[0][0].vec) << std::endl;
    _out << "#Class_Number= " << _m.__classNo << std::endl;
    
    for (int i = 0; i < _m.__classNo; i++)
    {
        int particlesNo = _m.__training_particles[i].size();
	    _out << "#Particles_No= " << particlesNo << std::endl;
	    for (int j = 0; j < particlesNo; j++)
	        _out << _m.__training_particles[i][j];
    }
    
    return _out;
}

void Classification_model::printShape() const
{
    for (int i=0; i<__training_particles.size(); i++)
        std::cout << "Class " << i << " has " << __training_particles[i].size()
                  << " elements" << std::endl;
}

/* Read -------------------------------------------------------------------- */
std::istream & operator >> (std::istream &_in, Classification_model &_m)
{
    std::string dummy;
    int classNo, vec_size;
    _in >> dummy;
    _in >> dummy >> _m.__micrographs_number;
    _in >> dummy;

    _m.__micrographs_scanned.resize(_m.__micrographs_number);
    for(int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__micrographs_scanned[i];
    _in >> dummy;

    _m.__particles_picked.resize(_m.__micrographs_number);
    for(int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__particles_picked[i];
    _in >> dummy;
    
    _m.__falsePositives.resize(_m.__micrographs_number);
    for(int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__falsePositives[i];
    _in >> dummy;

    _in >> dummy >> vec_size;
    _in >> dummy >> _m.__classNo;
    _m.__training_particles.resize(_m.__classNo);
    for (int i = 0; i < _m.__classNo; i++)
    {
        int particlesNo;
	    _in >> dummy >> particlesNo;
	    for (int j = 0; j < particlesNo; j++)
	    {
	        Particle *P = new Particle;
                P->read(_in, vec_size);
                _m.addParticleTraining(*P, i);
	    }
    }
    return _in;
}

/* Constructor ------------------------------------------------------------ */
QtWidgetMicrograph::QtWidgetMicrograph(QtMainWidgetMark *_mainWidget,
                                       QtFiltersController *_f,
                                       Micrograph *_m) :
                                       QWidget((QWidget*) _mainWidget)
{
    __filtersController              = _f;
    __m                              = NULL;
    __activeFamily                   = -1;
    __tilted                         = false;
    __learn_particles_done           = false;
    __autoselection_done             = false;
    __modelRootName                  = "";
    __numThreads                     = 1;
    __auto_label                     = -1;
    __gray_bins                      = 8;
    __radial_bins                    = 16;
    __piece_xsize                    = 512;
    __output_scale                   = 1;
    __highpass_cutoff                = 0.02;
    __penalization                   = 10;
    __minCost                        = 0;
    __reduction                      = (int)std::pow(2.0, __output_scale);
    __particle_radius                = 0;
    __min_distance_between_particles = 0.5 * __particle_radius;
    __mask_size                      = 2 * __particle_radius;
    __piece_overlap                  = 2 * __particle_radius;
    __scan_overlap                   = (int)(2 * __particle_radius * 0.9);
    __learn_overlap                  = (int)(2 * __particle_radius * 0.5);
    __classNo                        = 3;
    __is_model_loaded                = false;
    
#ifdef QT3_SUPPORT
    __gridLayout                     = new Q3VBoxLayout(this);
#else
    __gridLayout                     = new QVBoxLayout(this);
#endif
    __menuBar                        = new QMenuBar(this);
    __menuBar->setSeparator(QMenuBar::InWindowsStyle);
    __mImage                         = new QtImageMicrograph(0);
    __mImage->setWidgetMicrograph(this);
    __mImageOverview                 = new QtImageOverviewMicrograph(this);
    __mImageOverview->setWidgetMicrograph(this);
    __file_menu                      = NULL;
    __ellipse_radius                 = 5;
    __ellipse_type                   = MARK_CIRCLE;

    __mImage->setFiltersController(_f);
    __mImageOverview->setFiltersController(_f);

    connect(__mImageOverview, SIGNAL(signalSetCoords(int, int)),
            __mImage, SLOT(slotSetCoords(int, int)));
    connect(__mImage, SIGNAL(signalSetCoords(int, int)),
            __mImageOverview, SLOT(slotSetCoords(int, int)));
    connect(__mImage, SIGNAL(signalSetWidthHeight(int, int)),
            __mImageOverview, SLOT(slotSetWidthHeight(int, int)));
    connect(__mImage, SIGNAL(signalRepaint(void)),
            __mImageOverview, SLOT(repaint(void)));
    connect(__mImageOverview, SIGNAL(signalRepaint(void)),
            __mImage, SLOT(repaint(void)));
    connect(__mImage, SIGNAL(signalRepaint(void)),
            __mImage, SLOT(repaint(void)));
    connect(__mImageOverview, SIGNAL(signalRepaint(void)),
            __mImageOverview, SLOT(repaint(void)));
    connect(__mImage, SIGNAL(signalAddCoordOther(int, int, int)),
            this, SLOT(slotDrawEllipse(int, int, int)));

    connect(this, SIGNAL(signalActiveFamily(int)),
            __mImage, SLOT(slotActiveFamily(int)));
    connect(this, SIGNAL(signalActiveFamily(int)),
            __mImageOverview, SLOT(slotActiveFamily(int)));
    connect(this, SIGNAL(signalRepaint()),
            this, SLOT(slotRepaint()));

#ifdef QT3_SUPPORT
    Q3Accel *ctrl = new Q3Accel(this);
#else
    QAccel *ctrl  = new QAccel(this);
#endif
    ctrl->connectItem(ctrl->insertItem(Qt::Key_G + Qt::CTRL, 200),
                      this, SLOT(slotChangeContrast(void)));
#ifdef QT3_SUPPORT
    Q3Accel *ctrl2 = new Q3Accel(this);
#else
    QAccel *ctrl2  = new QAccel(this);
#endif
    ctrl2->connectItem(ctrl2->insertItem(Qt::Key_R + Qt::CTRL, 200),
                       this, SLOT(slotChangeCircleRadius(void)));

    setMicrograph(_m);

    __gridLayout->setMenuBar(__menuBar);
    __gridLayout->addWidget(__mImageOverview);

    openMenus();
}

/* Destructor--------------------------------------------------------------- */
QtWidgetMicrograph::~QtWidgetMicrograph()
{
    delete __mImage;
    delete __mImageOverview;
    delete __menuBar;
    delete __gridLayout;
}

/* Open all windows -------------------------------------------------------- */
void QtWidgetMicrograph::openAllWindows()
{
    __mImage->show();
    __mImageOverview->show();
}

/* Set Micrograph ---------------------------------------------------------- */
void QtWidgetMicrograph::setMicrograph(Micrograph *_m)
{
    if (_m != NULL)
    {
        __m = _m;
        __mImage->setMicrograph(_m);
        __mImageOverview->setMicrograph(_m);
    }
}

/* Open menus --------------------------------------------------------------*/
void QtWidgetMicrograph::openMenus()
{
    __file_menu = new QtFileMenu(this);
    connect(__mImage, SIGNAL(signalAddCoordOther(int, int, int)),
            __file_menu, SLOT(slotCoordChange()));

    QtFilterMenu *filterMenu = new QtFilterMenu(this);

    QtAutoMenu *autoMenu = new QtAutoMenu(this);

    addMenuItem("&File", (QtPopupMenuMark *)(__file_menu));
    addMenuItem("F&ilters", (QtPopupMenuMark *)(filterMenu));
    addMenuItem("&AutoSelection", (QtPopupMenuMark *)(autoMenu));

    connect(__file_menu, SIGNAL(signalAddFamily(const char *)),
            this, SLOT(slotAddFamily(const char*)));
    connect((QObject*)filterMenu, SIGNAL(signalAdjustContrast()),
            this, SLOT(slotChangeContrast(void)));
    connect((QObject*)filterMenu, SIGNAL(signalCrop()),
            this, SLOT(slotChangeCrop(void)));
    connect((QObject*)filterMenu, SIGNAL(signalAddFilter()),
            (QObject*)__filtersController, SLOT(slotAddFilter()));
    connect((QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            (QObject*)__filtersController, SLOT(slotCleanFilters()));
    connect((QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            this, SLOT(repaint()));

    // *** Add your own menus
}

/* Set the number of threads ----------------------------------------------- */
void QtWidgetMicrograph::setNumThreads(int _numThreads)
{
    __numThreads=_numThreads;
}

/* Learn particles --------------------------------------------------------- */
void QtWidgetMicrograph::learnParticles()
{
    if (__particle_radius==0) {
        QMessageBox helpmsg("Information",
            "You cannot learn particles without configuring the model",
            QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, false);
        helpmsg.exec();
        return;
    }
        
    std::cerr << "\n------------------Learning Phase-----------------------\n";
    createMask();
    
    // Find in the particle list the images to learn
    std::vector<int> all_idx;
    int num_part = __m->ParticleNo();
    for (int i = 0; i < num_part; i++)
        if (__m->coord(i).valid && __m->coord(i).label != __auto_label)
            all_idx.push_back(i);
    
    // If there is nothing to learn, return
    if (all_idx.size() == 0 && !__is_model_loaded)
    {
        std::cerr << "No valid particles marked." << std::endl;
	return;
    }
    
    // If we have already learnt or autoselected, delete the training vectors
    if (__learn_particles_done || __autoselection_done)
    {
        __training_model.clear();
        __training_model=__training_loaded_model;
    }
    
    // Actually learn
    buildVectors(all_idx, __training_model);
    getAutoTruePositives(__training_model);
    buildNegativeVectors(__training_model,false);
    getAutoFalsePositives(__training_model);
    buildNegativeVectors(__training_model,true);

    // Count the number of points scanned in a micrograph
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();
    int top = 0, left = 0, next_top = 0, next_left = 0;
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    int Nscanned=0;
    Matrix2D<double> piece;
    while (get_corner_piece(piece, top, left, skip_y,
                            next_skip_x, next_skip_y, next_top,
			    next_left, __piece_overlap, true))
    {
        top = next_top;
        left = next_left;
        skip_x /= __reduction;
        skip_y /= __reduction;
        int posx = 0, next_posx = 0, posy = 0, next_posy = 0;
        next_posx = posx = skip_x + XSIZE(mask) / 2;
        next_posy = posy = skip_y + YSIZE(mask) / 2;

        while (get_next_scanning_pos(piece, next_posx, next_posy, skip_x, 
            skip_y, __scan_overlap))
        {
            Nscanned++;
            posx = next_posx;
            posy = next_posy;
        }

        top = next_top;
        left = next_left;
        skip_x = next_skip_x;
        skip_y = next_skip_y;
    }
    __training_model.addMicrographScanned(Nscanned);
    __selection_model = __training_model;

    __learn_particles_done = true;    
    std::cerr << "Learning process finished..." << std::endl;
}

/* Build Selection Model---------------------------------------------------- */
void QtWidgetMicrograph::buildSelectionModel()
{
    __training_model.import_model(__training_loaded_model);
    __selection_model = __training_model;
}

/* produceFeatures --------------------------------------------------------- */
void QtWidgetMicrograph::produceFeatures(
    const Classification_model &_model,
    std::vector < Matrix2D<double> > &_features)
{
    _features.clear();
    std::vector < std::vector<Particle> > particles=
        _model.__training_particles;
    int vec_size=0;
    for (int i=0; i<particles.size(); i++)
        if (particles[i].size()!=0) vec_size=XSIZE(particles[i][0].vec);
    
    for(int j = 0; j < __classNo; j++)
    {
        int imax = particles[j].size();
        //if we do not have false-positives, then we only use 2 classes.
        if (imax == 0)
        {
            _features.resize(__classNo - 1);
            break;
        }
        _features.push_back(*(new Matrix2D<double>));
        _features[j].initZeros(imax, vec_size);
        for(int i = 0; i < imax; i++)
        {
            particles[j][i].vec.setRow();
            _features[j].setRow(i, particles[j][i].vec);	  
        }
    }
}

/* produceClassesProbabilities --------------------------------------------- */
void QtWidgetMicrograph::produceClassesProbabilities(
    const Classification_model &_model, Matrix1D<double> &probabilities)
{
    double micrographsScanned = 0.0;
    double particlesMarked = 0.0;
    if (_model.__training_particles.size()>0)
        particlesMarked = _model.__training_particles[0].size();
    double falsePositives = 0.0;
    if (_model.__training_particles.size()>2)
        falsePositives = _model.__training_particles[2].size();
    probabilities.initZeros(__classNo);
    
    for (int i = 0; i < _model.__micrographs_number; i++)
        if (i<_model.__micrographs_scanned.size())
            micrographsScanned += _model.__micrographs_scanned[i];
    micrographsScanned=XMIPP_MAX(micrographsScanned,
        2*(particlesMarked+falsePositives));
    int NtrainingVectors=0;
    for (int i=0; i<_model.__training_particles.size(); i++)
        NtrainingVectors+=_model.__training_particles[i].size();
    micrographsScanned=XMIPP_MAX(micrographsScanned, NtrainingVectors);
    
    probabilities(0) = particlesMarked / micrographsScanned;    
    probabilities(1) = (micrographsScanned - particlesMarked - falsePositives)/
                       micrographsScanned;
    probabilities(2) = falsePositives / micrographsScanned;
}

/* Automatic phase ----------------------------------------------------------*/
static pthread_mutex_t particleAdditionMutex = PTHREAD_MUTEX_INITIALIZER;

void * automaticallySelectParticlesThread(void * args)
{
    // Pick the input parameters
    AutomaticallySelectThreadParams * prm=
        (AutomaticallySelectThreadParams *)args;
    QtWidgetMicrograph *widgetmicrograph=prm->widgetmicrograph;
    int idThread=prm->idThread;

    int thVotes=1;
    if (widgetmicrograph->__selection_model.__training_particles.size()==2 ||
        (widgetmicrograph->__selection_model.__training_particles.size()==3 &&
         widgetmicrograph->__selection_model.__training_particles[2].size()==0))
        thVotes=8;
     //top,left corner of the piece
    int top = 0, left = 0, next_top = 0, next_left = 0;
    // If the piece available is small then include the scanned part
    // because we need bigger image for denoising but for scanning
    // particles we skip the already scanned part
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    Matrix1D<double> v;
    int N = 1, particle_idx = 0, Nscanned=0;
    Matrix2D<double> piece, original_piece;
    const Matrix2D<int> &mask = widgetmicrograph->__mask.get_binary_mask2D();
    std::vector< Particle > threadCandidates;
    do
    {
        bool isMine=(N%widgetmicrograph->__numThreads==idThread);
        bool pieceOK=widgetmicrograph->get_corner_piece(piece, top, left,
            skip_y, next_skip_x, next_skip_y, next_top, next_left,
            widgetmicrograph->__piece_overlap,isMine);
        if (!pieceOK) break;
        if (isMine)
        {
            if (idThread==0) std::cerr << ".";
            #ifdef DEBUG_MORE_AUTO
                std::cerr << "thread " << idThread
                          << " processing piece " << N << "...\n";
                std::cerr << "    (top,left)=" << top << "," << left
                          << " skip y,x=" << next_skip_y << "," << next_skip_x
                          << " next=" << next_top << "," << next_left << std::endl;
            #endif

            // Get a piece and prepare it
            if (!widgetmicrograph->prepare_piece(piece, original_piece))
            {
                top = next_top;
                left = next_left;
                std::cerr << "bad piece...skipping" << std::endl;
                N++;
                continue;
            }

            // Express the skip values in the reduced image
            skip_x /= widgetmicrograph->__reduction;
            skip_y /= widgetmicrograph->__reduction;
            #ifdef DEBUG_MORE_AUTO
                std::cerr << "Skip(y,x)=" << skip_y << "," << skip_x << std::endl;
            #endif

            // Scan this piece
            int posx = 0, next_posx = 0, posy = 0, next_posy = 0;
            next_posx = posx = skip_x + XSIZE(mask) / 2;
            next_posy = posy = skip_y + YSIZE(mask) / 2;

            while (widgetmicrograph->get_next_scanning_pos(
                piece, next_posx, next_posy, skip_x, skip_y,
                widgetmicrograph->__scan_overlap))
            {
                // COSS: Uncomment the next sentence for fast debugging
                // if (rnd_unif(0,1)<0.98) continue;
            
                #ifdef DEBUG_MORE_AUTO
                   std::cerr << "Pos(y,x)=" << posy << "," << posx
                             << " Micro(y,x)=" << posy*widgetmicrograph->__reduction + top
                             << "," << posx*widgetmicrograph->__reduction + left
                             << " Next pos(y,x)=" << next_posy << "," << next_posx
                             << std::endl;
                #endif

	        if (widgetmicrograph->build_vector(piece, original_piece,
                    posx, posy, v))
                {
                    double cost;
                    int votes=widgetmicrograph->__selection_model.isParticle(v,cost);
	            if (votes>thVotes)
	            {
                        #ifdef DEBUG_MORE_AUTO
		            std::cout << "Particle Found: "
                                << left + posx * widgetmicrograph->__reduction << ","
                                << top  + posy *widgetmicrograph->__reduction
                                << " votes=" << votes
                                << std::endl;
                            std::cout << "Press any key to continue...\n";
                            char c;
                            std::cin >> c;
                        #endif

                        // Build the Particle structure
                        Particle P;
                        P.x = left + posx * widgetmicrograph->__reduction;
                        P.y = top + posy * widgetmicrograph->__reduction;
                        P.idx = -1;
                        P.status = 1;
                        P.vec = v;
                        P.cost = cost;
                        threadCandidates.push_back(P);
                    }
                    Nscanned++;
                }

                // Go to next scanning position
                posx = next_posx;
                posy = next_posy;
            }
        }

        // Go to next piece in the micrograph
        top = next_top;
        left = next_left;
        skip_x = next_skip_x;
        skip_y = next_skip_y;
        N++;
    } while (true);

    int imax=threadCandidates.size();
    pthread_mutex_lock( &particleAdditionMutex );
    for (int i=0; i<imax; i++)
        widgetmicrograph->__auto_candidates.push_back(threadCandidates[i]);
    pthread_mutex_unlock( &particleAdditionMutex );
}

void QtWidgetMicrograph::automaticallySelectParticles()
{
    // Check that there is a valid model
    if (__selection_model.isEmpty())
    {
        std::cerr << "No model has been created." << std::endl;
        return;
    }

    std::cerr << "-----------------Automatic Phase--------------------------\n";

    // Initialize some variables
    __auto_candidates.resize(0);
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    // Get the training features and the a priori probabilities    
    std::vector < Matrix2D<double> > features;
    Matrix1D<double> probs;
    produceFeatures(__selection_model,features);
    produceClassesProbabilities(__selection_model,probs);

    // Initialize classifier
    __selection_model.initNaiveBayesEnsemble(features, probs, 8,
        __penalization, 10, 1, 1, "mm");
    #ifdef DEBUG_AUTO
        std::cout << "Probabilities of the classes:"
                  << probs.transpose() << std::endl;
    #endif

    // Automatically select particles with threads
    pthread_t * th_ids = new pthread_t[__numThreads];
    AutomaticallySelectThreadParams * th_args=
        new AutomaticallySelectThreadParams[__numThreads];
    for( int nt = 0 ; nt < __numThreads ; nt ++ )
    {
        th_args[nt].widgetmicrograph = this;
        th_args[nt].idThread = nt;
        pthread_create( &th_ids[nt], NULL,
            automaticallySelectParticlesThread, &th_args[nt]);
    }

    // Waiting for threads to finish
    for( int nt=0; nt<__numThreads; nt++)
        pthread_join(th_ids[nt], NULL);

    // Thread structures are not needed any more
    delete [] th_ids;
    delete [] th_args;

    // Sort particles by cost
    std::sort(__auto_candidates.begin(),__auto_candidates.end(),
        SAscendingParticleSort());

    #ifdef DEBUG_AUTO
       std::cerr << "Number of automatically selected particles = " 
                 << __auto_candidates.size() << std::endl;
    #endif

    // Reject the candidates that are pointing to the same particle
    int Nalive = reject_within_distance(__auto_candidates, __particle_radius,
        false);

    #ifdef DEBUG_AUTO
       std::cerr << "Number of automatically selected particles after distance rejection = " 
                 << __auto_candidates.size() << std::endl;
    #endif

    // Apply a second classifier for classifying between particle
    // and false positive. For that, remove the middle class (background)
    if (features.size()==3)
    {
        int imax;
        if (Nalive > 0)
        {
            __selection_model2.clear();
            __selection_model2.init(2);
            std::vector < Matrix2D<double> >::iterator featuresIterator=
                features.begin();
            featuresIterator++;
            features.erase(featuresIterator);
            probs(1)=probs(2);
            probs.resize(2);
            probs/=probs.sum();
            __selection_model2.initNaiveBayesEnsemble(features, probs, 8,
                __penalization,10,1,1,"mm");
            #ifdef DEBUG_AUTO
	        std::cout << "Second classification\n";
            #endif
            imax = __auto_candidates.size();
            Nalive = 0;
            for (int i = 0; i < imax; i++)
                if (__auto_candidates[i].status == 1)
                {
                    double p;
                    int votes=__selection_model2.isParticle(__auto_candidates[i].vec,p);
	            if (votes<8)
                    {
                        __auto_candidates[i].status=0;
                        __auto_candidates[i].cost=-1;
                        #ifdef DEBUG_AUTO
	                    std::cout << __auto_candidates[i].x << ", "
	                              << __auto_candidates[i].y
                                      << " is considered as a false positive\n";
                        #endif
                    }
                    else
                        Nalive++;
                }
        }

        // Apply a third classifier to distinguish between particles and
        // very tough particles
        if (Nalive > 0)
        {
            __selection_model3.clear();
            __selection_model3.init(2);
            // Remove from the error class, all those errors that
            // the previous classifier was able to classify correctly
            std::vector<int> toKeep;
            for (int i=0; i<YSIZE(features[1]); i++)
            {
                Matrix1D<double> trialFeatures;
                features[1].getRow(i,trialFeatures);
                double cost;
                int votes=__selection_model2.isParticle(trialFeatures,cost);
                if (votes<8)
                    toKeep.push_back(i);
            }

            if (toKeep.size()>0)
            {
                Matrix2D<double> difficultParticles;
                difficultParticles.initZeros(toKeep.size(),XSIZE(features[1]));
                FOR_ALL_ELEMENTS_IN_MATRIX2D(difficultParticles)
                    difficultParticles(i,j)=features[1](toKeep[i],j);
                int NErrors=YSIZE(features[1]);
                features.pop_back();
                features.push_back(difficultParticles);

                __selection_model3.initNaiveBayesEnsemble(features, probs, 8,
                    __penalization, 10, 1, 1, "mm");
                #ifdef DEBUG_AUTO
	            std::cout << "Third classification: " << YSIZE(difficultParticles)
                              << " difficult particles out of " << NErrors << "\n"
                              << "Before filtering there were " << Nalive << " particles\n";
                #endif
                imax = __auto_candidates.size();
                Nalive = 0;
                for (int i = 0; i < imax; i++)
                    if (__auto_candidates[i].status == 1)
                    {
                        double p;
                        int votes=__selection_model3.isParticle(__auto_candidates[i].vec,p);
	                if (votes<8)
                        {
                            __auto_candidates[i].status=0;
                            __auto_candidates[i].cost=-1;
                            #ifdef DEBUG_AUTO
	                        std::cout << __auto_candidates[i].x << ", "
	                                  << __auto_candidates[i].y
                                          << " is considered as a false positive 2\n";
                            #endif
                        }
                        else
                            Nalive++;
                    }
            }
        }
    }

    if (Nalive>0)
    {
        // Add a family to the micrograph is necessary
        bool addFamily=true;
        for (int i = 0; i < __m->LabelNo(); i++)
            if (__m->get_label(i)=="auto")
            {
                addFamily=false;
                break;
            }
        if (addFamily)
        {
            __m->add_label("auto");
            emit signalAddFamily("auto");
        }
        __auto_label=-1;
        for (int i = 0; i < __m->LabelNo(); i++)
            if (__m->get_label(i)=="auto")
            {
                __auto_label=i;
                break;
            }

        int imax = __auto_candidates.size();
        // Get the maximum and minimum cost
        bool first=true;
        double minCost, maxCost;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].status == 1)
            {
                if (!isinf(__auto_candidates[i].cost))
                {
                    if (first || __auto_candidates[i].cost<minCost)
                    {
                        minCost=__auto_candidates[i].cost;
                        first=false;
                    }
                    if (first || __auto_candidates[i].cost>maxCost)
                    {
                        maxCost=__auto_candidates[i].cost;
                        first=false;
                    }
                }
            }

        // Insert selected particles in the result
        int idxMicrograph=0;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].status == 1)
            {
                __auto_candidates[i].idx = idxMicrograph;
                #ifdef DEBUG_AUTO
	            std::cout << "Particle coords " << __auto_candidates[i].x << ", "
	                      << __auto_candidates[i].y << std::endl;
                #endif
                if (isinf(__auto_candidates[i].cost))
                    __auto_candidates[i].cost=1;
                else
                    __auto_candidates[i].cost=
                        (__auto_candidates[i].cost-maxCost)/
                        (minCost-maxCost);
                getMicrograph()->add_coord(__auto_candidates[i].x, 
	                                   __auto_candidates[i].y,
                                           __auto_label,
                                           __auto_candidates[i].cost);
                idxMicrograph++;
            }
            else
                __auto_candidates[i].status = 0;
    }

    emit signalRepaint();
    __autoselection_done = true;
    std::cout << "\nAutomatic process finished. Number of particles found: "
              << Nalive << std::endl;

    if (Nalive>0)
    {
        ScrollParam* param_window;
        param_window = new ScrollParam(0, 1, 0.01, "Quality",
            "Restrict selection", 0, "new window", Qt::WDestructiveClose);

        // Connect its output to my input (set_spacing)
        connect( param_window, SIGNAL(new_value(float)), this,
            SLOT(slotRestrictSelection(float)));

        // Show
        param_window->setFixedSize(300,150);
        param_window->show();
    }
    slotActiveFamily(0);
}

/* Restrict Selection ------------------------------------------------------ */
void QtWidgetMicrograph::slotRestrictSelection(float _cost)
{
    __minCost=__mImage->__minCost=__mImageOverview->__minCost=_cost;
    int Nremaining=0;
    int imax = __auto_candidates.size();
    for (int i = 0; i < imax; i++)
        if (__auto_candidates[i].status == 1 && __auto_candidates[i].cost>_cost)
            Nremaining++;
    std::cout << Nremaining << " particles remain" << std::endl;
    emit signalRepaint();
}

/* Add family -------------------------------------------------------------- */
int QtWidgetMicrograph::add_family(std::vector<Particle> &_list,
                                   const std::string &_family_name)
{
    int ilabel = __m->add_label(_family_name);
    int imax = _list.size();
    for (int i = 0; i < imax; i++)
    {
        int idx = __m->add_coord(_list.at(i).x, _list.at(i).y, ilabel,
            _list.at(i).cost);
        _list.at(i).idx = idx;
    }
    return ilabel;
}

/* Create mask for learning particles ---------------------------------------*/
void QtWidgetMicrograph::createMask()
{
    if (XSIZE(__mask.get_binary_mask2D()) != 0) return;

    int xsize = __mask_size / __reduction;
    int ysize = xsize;
    int radius = __particle_radius / __reduction;

    __mask.type = BINARY_CIRCULAR_MASK;
    __mask.mode = INNER_MASK;
    __mask.R1 = radius;
    __mask.resize(xsize, ysize);
    __mask.generate_2Dmask();
    __mask.get_binary_mask2D().setXmippOrigin();

    classifyMask();
}

/* Classify the mask pixels -------------------------------------------------*/
void QtWidgetMicrograph::classifyMask()
{
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();
    if (XSIZE(mask) == 0) return;
    double max_radius_particle = __particle_radius / __reduction;

    // Determine max_radius
    double max_radius = max_radius_particle;

    // Initialize some variables
    // 6 is the minimum radius to be informative
    double radial_step = (max_radius - 6) / __radial_bins;

    Matrix2D<int> *classif1 = new Matrix2D<int>;
    classif1->resize(mask);
    classif1->initConstant(-1);

    Matrix2D<int> *classif2 = new Matrix2D<int>;
    classif2->resize(mask);
    classif2->initConstant(-1);
    const double deltaAng=(PI/8.0);

    Matrix1D<int> Nrad(__radial_bins);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
    {
        double radius = sqrt((double)(i * i + j * j));
        double angle = atan2((double)i, (double)j);
        if (angle < 0) angle += 2 * PI;

        if (radius < max_radius)
        {
            // Classif1 is the classification for the radial mass distribution
            int radius_idx;
            if (radius > 6) radius_idx = XMIPP_MIN(__radial_bins - 1, 1 + 
	                                 FLOOR((radius - 6) / radial_step));
            else            radius_idx = 0;
            (*classif1)(i, j) = radius_idx;
            Nrad(radius_idx)++;

            // Classif2 is the classification by angles
            (*classif2)(i, j) = FLOOR(angle/deltaAng);
        }
    }
    __mask_classification.push_back(classif1);
    __mask_classification.push_back(classif2);

    // Create the holders for the radius values in classif1
    for (int i = 0; i < __radial_bins; i++)
    {
        Matrix1D<int> *aux1 = new Matrix1D<int>;
        aux1->initZeros(Nrad(i));
        __radial_val.push_back(aux1);
    }

    // Create the holders for the radius values in classif1 and classif2
    int angleBins=classif2->computeMax()+1;
    for (int i = 0; i < angleBins; i++)
    {
        Matrix1D<double> *aux2 = new Matrix1D<double>;
        aux2->initZeros(__radial_bins);
        __sector.push_back(aux2);
        
        Matrix1D<int> *iaux2 = new Matrix1D<int>;
        iaux2->initZeros(__radial_bins);
        __Nsector.push_back(iaux2);
    }

    // Compute the number of elements in each sector
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*classif1)
    {
        int idx1 = (*classif1)(i, j);
        if (idx1 != -1)
        {
            int idx2 = (*classif2)(i, j);
            if (idx2 != -1)
                (*__Nsector[idx2])(idx1)++;
        }
    }

    for (int i = 0; i < __radial_bins; i++)
    {
        Matrix1D<double> *aux3 = new Matrix1D<double>;
        aux3->initZeros(angleBins);
        __ring.push_back(aux3);
    }
        
    #ifdef DEBUG_CLASSIFY
        ImageXmipp save;
        typeCast(*classif1, save());
        save.write("PPPmaskClassification1.xmp");

        typeCast(*classif2, save());
        save.write("PPPmaskClassification2.xmp");
    #endif
}

/* Build training vectors ---------------------------------------------------*/
void QtWidgetMicrograph::buildVectors(std::vector<int> &_idx,
                                      Classification_model &_model)
{
    std::cerr << "Building " << _idx.size()
              << " particle vectors for this image. Please wait..." << std::endl;
    _model.addMicrographItem();
    int num_part = _idx.size();
    int width, height;
    __m->size(width, height);
    Matrix1D<double> v;
    int numParticles = 0;

    Matrix1D<char> visited(num_part);
    while (visited.sum() < num_part)
    {
        int part_i = 0;

        // get the first un-visited particle in the array
        while (visited(part_i) == 1) part_i++;
        if (part_i >= num_part) break;
        visited(part_i) = 1;

        // Get the piece containing that particle
        int part_idx = _idx.at(part_i);
        int x = __m->coord(part_idx).X;
        int y = __m->coord(part_idx).Y;
        int posx, posy;
        Matrix2D<double> piece, original_piece;
        get_centered_piece(piece, x, y, posx, posy);

        // Denoise, reduce, reject outliers and equalize histogram
        bool success = prepare_piece(piece, original_piece);
       
        if (!success) continue;
        posx = ROUND(posx / __reduction);
        posy = ROUND(posy / __reduction);
       
        //make vector from this particle
        success = build_vector(piece, original_piece, posx, posy, v);
        if (success)
        {
            Particle p;
            p.x = x;
            p.y = y;
            p.idx = part_idx;
            p.vec = v;
            p.status = 1;
            p.cost = 1.0;
            _model.addParticleTraining(p, 0);
            numParticles++;
        }

        // make vector from the neighbours
        std::vector< Matrix1D<int> > nbr;
        nbr.reserve(num_part);

        find_neighbour(piece, _idx, part_i, x, y, posx, posy, visited, nbr);

        for (int i = 0; i < nbr.size(); i++)
        {
            part_i = nbr.at(i)(0);
            part_idx = _idx.at(part_i);
            posx = nbr.at(i)(1);
            posy = nbr.at(i)(2);
            success = build_vector(piece, original_piece, posx, posy, v);
            visited(part_i) = 1;
            if (success)
            {
                Particle p;
                p.x = __m->coord(part_idx).X;
                p.y = __m->coord(part_idx).Y;
                p.idx = part_idx;
                p.vec = v;
                p.status = 1;
                p.cost = 1.0;
                _model.addParticleTraining(p, 0);
                numParticles++;
            }
        }
    }
    _model.addParticlePicked(numParticles);
}

/* Build vector from non particles------------------------------------------ */
void QtWidgetMicrograph::buildNegativeVectors(Classification_model &_model,
    bool checkForPalsePostives)
{
    if (checkForPalsePostives)
        std::cerr << "Building automatic false positives ..." << std::endl;
    else
        std::cerr << "Building non particles ..." << std::endl;
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    // Setup a classification model with the already known data
    std::vector < Matrix2D<double> > features;
    Matrix1D<double> probs;
    if (checkForPalsePostives)
    {
        // Gather all information for classification
        __selection_model = __training_model;

        // Prepare data to be classified
        produceFeatures(__selection_model,features);
        produceClassesProbabilities(__selection_model,probs);

        // Initialize classifier
        __selection_model.initNaiveBayesEnsemble(features, probs, 8,
            __penalization, 10, 1, 1, "mm");
    }

    // top,left corner of the piece
    int top = 0, left = 0, next_top = 0, next_left = 0;

    // If the piece available is small then include the scanned part
    // because we need bigger image for denoising but for scanning
    // particles we skip the already scanned part
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    Matrix1D<double> v;
    
    int N = 1, Nnonparticles=0, Nfalsepositives=0;
 
    // We do not want any overlap for this process,since it is only for 
    // counting the non particles and calculating their features. For 
    // the process of automatic selecting we will want an overlap so we
    // do not miss any particle.
    Matrix2D<double> piece, original_piece;
    while (get_corner_piece(piece, top, left, skip_y,
                            next_skip_x, next_skip_y, next_top, next_left, 0,
                            true))
    {
        // Get a piece and prepare it
        if (!prepare_piece(piece, original_piece))
        {
            top = next_top;
            left = next_left;
            std::cerr << "bad piece...skipping" << std::endl;
            N++;
            continue;
        }

        // Express the skip values in the reduced image
        skip_x /= __reduction;
        skip_y /= __reduction;

        // Scan this piece
        int posx = 0, next_posx = 0, posy = 0, next_posy = 0;
        next_posx = posx = skip_x + XSIZE(mask) / 2;
        next_posy = posy = skip_y + YSIZE(mask) / 2;

	// We do not want any overlap for this process, since it is only for
	// counting the non particles and calculating their features. For
	// the process of automatic selecting we will want an overlap so we
	// do not miss any particle.
        while (get_next_scanning_pos(piece, next_posx, next_posy, skip_x,
            skip_y, __learn_overlap))
        {
            // Check if there is any particle around
	    if (!anyParticle(left + posx  * __reduction,
	                      top + posy  * __reduction,
                              XSIZE(mask) * __reduction))
            {
                if (build_vector(piece, original_piece, posx, posy, v))
                {
                    // Build the Particle structure
                    Particle P;
                    P.x = left + posx * __reduction;
                    P.y = top + posy * __reduction;
                    P.idx = -1;
                    P.status = 1;
                    P.vec = v;
                    P.cost = -1;
                    double cost;
                    if (!checkForPalsePostives)
                    {
                        _model.addParticleTraining(P, 1);
                        Nnonparticles++;
                    }
                    else
                    {
                        int votes=__selection_model.isParticle(v,cost);
                        if (votes>5)
                        {
                            _model.addParticleTraining(P, 2);
                            Nfalsepositives++;
                        }
                    }
                }
            }
	    // Go to next scanning position
            posx = next_posx;
            posy = next_posy;
	}

        // Go to next piece in the micrograph
        top = next_top;
        left = next_left;
        skip_x = next_skip_x;
        skip_y = next_skip_y;
        N++;
    }
    if (checkForPalsePostives)
        std::cout << Nfalsepositives << " false positive chosen\n";
    else
        std::cout << Nnonparticles << " non particles randomly chosen\n";
}

bool QtWidgetMicrograph::anyParticle(int posx, int posy, int rect_size)
{
    int num_part = __m->ParticleNo();
    std::vector<Particle_coords> selected_particles = __m->Particles();
    for (int i = 0; i < num_part; i++)
    {
        if (__m->coord(i).valid && __m->coord(i).label != __auto_label)
	    {
	        int _x = selected_particles[i].X;
                int _y = selected_particles[i].Y;
    
                if((_x > posx - rect_size) && (_x < posx + rect_size))
	        {
	            if((_y > posy - rect_size) && (_y < posy + rect_size))
   	            return true;
	        }   
	    }
    }
    return false;
}

/* Build classification vector --------------------------------------------- */
bool QtWidgetMicrograph::build_vector(const Matrix2D<double> &piece,
    const Matrix2D<double> &original_piece, int _x, int _y,
    Matrix1D<double> &_result)
{
    #ifdef DEBUG_BUILDVECTOR
        std::cout << "build_vector(" << _x << "," << _y << "," << "_result)" << std::endl;
    #endif

    // Resize the output and make same aliases
    int angleBins=__sector.size();
    _result.initZeros(__radial_bins*(__gray_bins-1) // radial histograms
        +(angleBins-1)+(angleBins-1)*angleBins // sector correlations
        +(2*__radial_bins-19)*angleBins // ring correlations
        );
    const Matrix2D<int> &mask =         __mask.get_binary_mask2D();
    const Matrix2D<int> &classif1 =  (*(__mask_classification[0]));
    const Matrix2D<int> &classif2 =  (*(__mask_classification[1]));

    if (STARTINGX(mask) + _x < STARTINGX(piece)) return false;
    if (STARTINGY(mask) + _y < STARTINGY(piece)) return false;
    if (FINISHINGX(mask) + _x > FINISHINGX(piece)) return false;
    if (FINISHINGY(mask) + _y > FINISHINGY(piece)) return false;

    #ifdef DEBUG_IMG_BUILDVECTOR
        bool debug_go = false;
        ImageXmipp save, savefg, saveOrig;
        if (true)
        {
            save() = piece;
            save.write("PPP0.xmp");
            save().initZeros(YSIZE(mask), XSIZE(mask));
            STARTINGY(save()) = STARTINGY(mask);
            STARTINGX(save()) = STARTINGX(mask);
            savefg() = save();
            savefg().initConstant(-1);
            saveOrig() = savefg();
            debug_go = true;
        }
    #endif

    Matrix1D<int> radial_idx(__radial_bins);
    Matrix2D<double> particle;
    particle.initZeros(YSIZE(mask), XSIZE(mask));
    STARTINGY(particle) = STARTINGY(mask);
    STARTINGX(particle) = STARTINGX(mask);

    // Copy __radial_val, __sector, and __ring into local structures
    // so that threads are applicable
    std::vector< Matrix1D<int> > radial_val;
    for (int i=0; i<__radial_val.size(); i++)
    {
        Matrix1D<int> dummyi=*(__radial_val[i]);
        radial_val.push_back(dummyi);
    }

    std::vector< Matrix1D<double> > sector;
    Matrix1D<double> dummyd=*(__sector[0]);
    for (int i=0; i<__sector.size(); i++)
        sector.push_back(dummyd);

    std::vector< Matrix1D<double> > ring;
    dummyd=*(__ring[0]);
    for (int i=0; i<__ring.size(); i++)
        ring.push_back(dummyd);

    // Put the image values into the corresponding radial bins
    FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
    {
        int val = (int)piece(_y + i, _x + j);
        bool foreground = mask(i, j);

        int idx1 = classif1(i, j);
        if (idx1 != -1)
        {
            radial_val[idx1](radial_idx(idx1)++) = val;
            int idx2 = classif2(i, j);
            if (idx2 != -1)
                sector[idx2](idx1) += val;
        }

        // Get particle
        #ifdef DEBUG_IMG_BUILDVECTOR
            if (debug_go)
            {
                save(i, j) = val;
                if (foreground) 
	        {
	            savefg(i, j) = val;
	            saveOrig(i, j) =  original_piece(_y+i, _x+j);
	        }
            }
        #endif
    }

    // Compute the sector averages and reorganize the data in rings
    for (int j = 0; j < angleBins; j++)
        FOR_ALL_ELEMENTS_IN_MATRIX1D(sector[j])
        {
            sector[j](i)/=(*(__Nsector[j]))(i);
            ring[i](j)=sector[j](i);
        }

    // Compute the histogram of the radial bins and store them  
    int idx_result=0;
    for (int i = 0; i < __radial_bins; i++)
    {
        histogram1D hist;
        compute_hist(radial_val[i], hist, 0, __gray_bins - 1, __gray_bins);
        for (int j = 0; j < __gray_bins - 1; j++)
            _result(idx_result++) = hist(j);
    }

    // Compute the correlation of the sectors
    for (int step = 1; step<angleBins; step++)
    {
        Matrix1D<double> sectorCorr;
        sectorCorr.initZeros(angleBins);
        for (int i = 0; i<angleBins; i++)
        {
            sectorCorr(i)=correlation_index(
                sector[i],
                sector[intWRAP(i+step,0,angleBins-1)]);
        }
        _result(idx_result++) = sectorCorr.computeAvg();
        Matrix1D<double> sectorAutocorr;
        correlation_vector_no_Fourier(sectorCorr,sectorCorr,sectorAutocorr);
        for (int j = 0; j < XSIZE(sectorAutocorr); j++)
            _result(idx_result++) = sectorAutocorr(j);
    }

    // Compute the correlation in rings
    for (int step = 1; step<=2; step++)
         for (int i = 8; i<__radial_bins-step; i++)
         {
             Matrix1D<double> ringCorr;
             correlation_vector_no_Fourier(ring[i],ring[i+step],ringCorr);
             for (int j = 0; j < XSIZE(ringCorr); j++)
                 _result(idx_result++) = ringCorr(j);
         }

    #ifdef DEBUG_IMG_BUILDVECTOR
        if (debug_go)
        {
            save.write("PPP1.xmp");
            savefg.write("PPP2.xmp");
            saveOrig.write("PPP3.xmp");

            save().initZeros(savefg());
            FOR_ALL_ELEMENTS_IN_MATRIX2D(save())
            {
                int idx1 = classif1(i, j);
                if (idx1 == -1) continue;
                int idx2 = classif2(i, j);
                if (idx2 == -1) continue;
                save(i,j)=sector[idx2](idx1);
            }
            save.write("PPP4.xmp");

            save().initZeros();
            FOR_ALL_ELEMENTS_IN_MATRIX2D(save())
            {
                int idx1 = classif1(i, j);
                if (idx1 == -1) continue;
                int idx2 = classif2(i, j);
                if (idx2 == -1) continue;
                save(i,j)=ring[idx1](idx2);
            }
            save.write("PPP5.xmp");
        }
    #endif

    #ifdef DEBUG_BUILDVECTOR
        Matrix2D<double> A;
        A.resize(angleBins,XSIZE(sector[0]));
        FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
            A(i,j)=sector[i](j);
        A.write("PPPsector.txt");
        std::cout << _result.transpose() << std::endl;
    #endif
    
    #if defined(DEBUG_BUILDVECTOR) || defined(DEBUG_IMG_BUILDVECTOR)
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
    #endif
    return true;
}

/* Get piece --------------------------------------------------------------- */
// to get the piece containing (x,y) of size xsize,ysize
// return the position of x,y in the piece in posx,posy
void QtWidgetMicrograph::get_centered_piece(Matrix2D<double> &piece,
    int _x, int _y, int &_posx, int &_posy)
{
    piece.resize(__piece_xsize, __piece_xsize);
    int startx = _x - ROUND(__piece_xsize / 2);
    int endx   = _x + ROUND(__piece_xsize / 2);
    int starty = _y - ROUND(__piece_xsize / 2);
    int endy   = _y + ROUND(__piece_xsize / 2);
    int maxx, maxy;
    __m->size(maxx, maxy);
    _posx = ROUND(__piece_xsize / 2);
    _posy = ROUND(__piece_xsize / 2);
    
    // boundary adjustments
    if (startx < 0)
    {
        _posx += startx;
        startx = 0;
        endx = __piece_xsize - 1;
    }
    if (starty < 0)
    {
        _posy += starty;
        starty = 0;
        endy = __piece_xsize - 1;
    }
    if (endx > maxx - 1)
    {
        _posx += endx - (maxx - 1);
        endx = maxx - 1;
        startx = endx - __piece_xsize;
    }
    if (endy > maxy - 1)
    {
        _posy += endy - (maxy - 1);
        endy = maxy - 1;
        starty = endy - __piece_xsize;
    }
    
    //read the matrix from the micrograph
    for (int i = 0; i < __piece_xsize; i++)
        for (int j = 0; j < __piece_xsize; j++)
            piece(i, j) = (*__m)(startx + j, starty + i);	    
}

// Get a piece whose top-left corner is at the desired position (if possible)
bool QtWidgetMicrograph::get_corner_piece(
    Matrix2D<double> &piece,
    int _top, int _left, int _skip_y,
    int &_next_skip_x, int &_next_skip_y, int &_next_top, int &_next_left,
    int overlap, bool copyPiece)
{
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    int maxx, maxy;
    __m->size(maxx, maxy);

    if (maxx < _left + __piece_xsize || maxy < _top + __piece_xsize)
        return false;

    _next_skip_x = _next_skip_y = 0;
    bool increase_Y = false;
    if (_left + __piece_xsize != maxx)
    {
        _next_left = _left + __piece_xsize - overlap;
        if (_next_left + __piece_xsize >= maxx)
            _next_left = maxx - __piece_xsize;
    }
    else
    {
        _next_left = 0;
        increase_Y = true;
    }
    if (increase_Y)
    {
        if (_top + __piece_xsize != maxy)
        {
            _next_top = _top + __piece_xsize - overlap;
            if (_next_top + __piece_xsize >= maxy)
                _next_top = maxy - __piece_xsize;
        }
        else
        {
            _next_top = maxy;
        }
    }

    //read the matrix from the micrograph
    if (copyPiece)
    {
        piece.resize(__piece_xsize, __piece_xsize);
        for (int i = 0; i < __piece_xsize; i++)
            for (int j = 0; j < __piece_xsize; j++)
                piece(i, j) = (*__m)(_left + j, _top + i);
    }
	    
    return true;
}

/* Prepare piece ----------------------------------------------------------- */
static pthread_mutex_t preparePieceMutex = PTHREAD_MUTEX_INITIALIZER;
bool QtWidgetMicrograph::prepare_piece(Matrix2D<double> &piece,
    Matrix2D<double> &original_piece)
{
    original_piece = piece;
    #ifdef DEBUG_PREPARE    
        ImageXmipp save;
        save() = piece;
        save.write("PPPpiece0.xmp");
    #endif    

    // Denoise the piece
    Denoising_parameters denoiser;
    denoiser.denoising_type = Denoising_parameters::BAYESIAN;
    denoiser.scale = 3;
    denoiser.output_scale = 1;
//    denoiser.scale = __output_scale + 3;
//    denoiser.output_scale = __output_scale;
    denoiser.produce_side_info();
    denoiser.denoise(piece);
    if (!(piece(0, 0) == piece(0, 0))) return false;
    #ifdef DEBUG_PREPARE    
        save() = piece;
        save.write("PPPpiece1.xmp");
    #endif    
    if (__output_scale==0)
    {
        Matrix2D<double> auxPiece;
        piece.pyramidExpand(auxPiece);
        piece=auxPiece;
        #ifdef DEBUG_PREPARE    
            save() = piece;
            save.write("PPPpiece1.5.xmp");
        #endif
    }

    // Band pass filter
    pthread_mutex_lock( &preparePieceMutex );
    FourierMask Filter;
    Filter.FilterShape = RAISED_COSINE;
    Filter.FilterBand = BANDPASS;
    Filter.w1 = __highpass_cutoff;
    Filter.w2 = 1.0/(__particle_radius/(__reduction*5.0));
    Filter.raised_w = XMIPP_MIN(0.02, __highpass_cutoff);
    Filter.generate_mask(piece);
    Filter.apply_mask_Space(piece);
    STARTINGX(piece) = STARTINGY(piece) = 0;
    #ifdef DEBUG_PREPARE    
        save() = piece;
        save.write("PPPpiece2.xmp");
    #endif    
    pthread_mutex_unlock( &preparePieceMutex );
    
    // Reject 5% of the outliers
    reject_outliers(piece, 5.0);

    // Equalize histogram    
    histogram_equalization(piece, __gray_bins);
    
    #ifdef DEBUG_PREPARE    
        save() = piece;
        save.write("PPPpiece3.xmp");
    #endif    
    original_piece.selfScaleToSize(YSIZE(piece), XSIZE(piece));

    // Reject 5% of the outliers
    reject_outliers(original_piece, 5.0);
    original_piece.statisticsAdjust(0,1);
    original_piece-=original_piece.computeAvg();

    #ifdef DEBUG_PREPARE    
        save() = original_piece;
        save.write("PPPpiece4.xmp");
    #endif
    return true;
}

/* Get neighbours ---------------------------------------------------------- */
//To get the neighbours and their positions in the piece image
void QtWidgetMicrograph::find_neighbour(const Matrix2D<double> &piece,
    std::vector<int> &_idx, int _index,
    int _x, int _y, int _posx, int _posy, Matrix1D<char> &_visited,
    std::vector< Matrix1D<int> > &_nbr)
{
    int piece_xsize = XSIZE(piece);
    int piece_ysize = YSIZE(piece);
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    // If all the particles are visited
    if (_visited.sum() == XSIZE(_visited)) return;

    int current_part = _index + 1;
    int current_part_idx = _idx.at(current_part);
    _nbr.clear();
    _nbr.reserve(XSIZE(_visited));

    int top = CEIL((double)_y / __reduction) - _posy;
    int left = CEIL((double)_x / __reduction) - _posx;
    int bottom = top + FLOOR((double)piece_ysize / __reduction)-1;
    int right = left + FLOOR((double)piece_xsize / __reduction)-1;
    int xmask2 = CEIL(XSIZE(mask) / 2);
    int ymask2 = CEIL(YSIZE(mask) / 2);

    while (current_part < XSIZE(_visited))
    {
        //find the next unvisited particle
        while (_visited(current_part))
        {
            current_part++;
            if (current_part==XSIZE(_visited)) break;
        }
        if (current_part==XSIZE(_visited)) break;
        current_part_idx = _idx.at(current_part);

        //check if it is neighbour or not
        int nx = ROUND((double)__m->coord(current_part_idx).X / __reduction);
        int ny = ROUND((double)__m->coord(current_part_idx).Y / __reduction);
        if ((nx - xmask2 > left) &&
            (ny - ymask2 > top) &&
            (nx + xmask2 < right) &&
            (ny + ymask2 < bottom))
        {
            Matrix1D<int> current_nbr;
            current_nbr.initZeros(3);
            current_nbr(0) = current_part;
            current_nbr(1) = nx - left;
            current_nbr(2) = ny - top;
            _nbr.push_back(current_nbr);
        }
        current_part++;
    }
}

/* Get next scanning position ---------------------------------------------- */
bool QtWidgetMicrograph::get_next_scanning_pos(
    const Matrix2D<double> &piece,
    int &_x, int &_y, int _skip_x, int _skip_y, int overlap)
{
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    if (_x + XSIZE(mask) / 2 > XSIZE(piece) ||
        _y + YSIZE(mask) / 2 > YSIZE(piece)) 
	return false;

    if (_x == 0 && _y == 0)
    {
        _x += XSIZE(mask) - overlap / __reduction;
        _y += YSIZE(mask) - overlap / __reduction;
    }
    else
    {
        int nextx = _x + XSIZE(mask) - overlap / __reduction;
        int nexty = _y + YSIZE(mask) - overlap / __reduction;
	   
	if (nextx + (XSIZE(mask) / 2) > XSIZE(piece))
        {
	    if (nexty + (YSIZE(mask) / 2) > YSIZE(piece))
            {
                _x = nextx;
                _y = nexty;
            }
            else
            {
                _x = _skip_x + (XSIZE(mask) / 2);
                _y = nexty;
            }
        }
        else
        {
            _x = nextx;
        }
    }
    return true;
}

/* Filter particles --------------------------------------------------------*/
//To calculate the euclidean distance between to points
double euclidean_distance(const Particle &p1, const Particle &p2)
{
    return sqrt((double)(p1.x -p2.x)*(p1.x - p2.x) + (double)(p1.y - p2.y)*
                (p1.y - p2.y));
}

int QtWidgetMicrograph::reject_within_distance(
    std::vector<Particle> &_Input, double _min_dist,
    bool _reject_both)
{
    int imax = _Input.size();
    int n = 0;
    for (int i = 0; i < imax; i++)
    {
        if (_Input.at(i).status == 0) continue;
        for (int j = i + 1; j < imax; j++)
        {
            if (_Input.at(j).status == 0) continue;
            double dist = euclidean_distance(_Input.at(i), _Input.at(j));
            if (dist < _min_dist)
            {
                _Input.at(j).status = 0;
                _Input.at(j).cost = -1;
                if (_reject_both)
                {
                    _Input.at(i).status = 0;
                    _Input.at(i).cost = -1;
                }
            }
        }
        if (_Input.at(i).status == 1) n++;
    }
    return n;
}

/* Correct particles ------------------------------------------------------- */
void QtWidgetMicrograph::move_particle(int _idx)
{
    if (__autoselection_done)
    {
        int imax = __auto_candidates.size();
        int idxAuto=-1;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].idx == _idx)
            {
                idxAuto=i;
                break;
            }
        if (idxAuto==-1) return;

        double diffx=__auto_candidates[idxAuto].x - __m->coord(_idx).X;
        double diffy=__auto_candidates[idxAuto].y - __m->coord(_idx).Y;
        if (sqrt(diffx*diffx+diffy*diffy)>2*__particle_radius*0.1)
             __rejected_particles.push_back(__auto_candidates[idxAuto]);
        __auto_candidates[idxAuto].x = __m->coord(_idx).X;
        __auto_candidates[idxAuto].y = __m->coord(_idx).Y;
    }
}

/* Delete particles ---------------------------------------------------------*/
void QtWidgetMicrograph::delete_particle(int _idx)
{
    if (__autoselection_done)
    {
        // Look for this particle in the autolist
        int imax = __auto_candidates.size();
        int idxAuto=-1;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].idx == _idx)
            {
                idxAuto=i;
                break;
            }
        if (idxAuto==-1) return;

        __auto_candidates[idxAuto].status = 0;
        __auto_candidates[idxAuto].cost = -1;
        __rejected_particles.push_back(__auto_candidates[idxAuto]);
    }
}

/* Get the false-positive particles----------------------------------------- */
void QtWidgetMicrograph::getAutoFalsePositives(
    Classification_model &_training_model)
{
    // Add the false positives given by the threshold
    int Nrejected=0;
    int imax = __auto_candidates.size();
    for (int i = 0; i < imax; i++)
        if (__auto_candidates[i].status == 1 &&
            __auto_candidates[i].cost<__minCost)
        {
            _training_model.addParticleTraining(__auto_candidates[i], 2);
            Nrejected++;
        }

    // Add the manually selected false positives
    imax = __rejected_particles.size();
    _training_model.addFalsePositives(imax);
    for (int i = 0; i < imax; i++)
        _training_model.addParticleTraining(__rejected_particles[i], 2);
    Nrejected+=imax;
    std::cout << Nrejected << " false positives are considered\n";
}

/* Get the false-positive particles----------------------------------------- */
void QtWidgetMicrograph::getAutoTruePositives(
    Classification_model &_training_model)
{
    // Add the true positives given by the threshold
    int Naccepted=0;
    int imax = __auto_candidates.size();
    for (int i = 0; i < imax; i++)
        if (__auto_candidates[i].status == 1 &&
            __auto_candidates[i].cost>=__minCost)
        {
            _training_model.addParticleTraining(__auto_candidates[i], 0);
            Naccepted++;
        }
    std::cout << Naccepted << " true positives are considered\n";
}

/* Load models ------------------------------------------------------------- */
void QtWidgetMicrograph::loadModels()
{
    // Get the rootname
    bool ok;
    QString qfn_root = QInputDialog::getText("Loading model",
                       "Model", QLineEdit::Normal,
                       "Model", &ok);
    if (!ok || qfn_root.isEmpty()) return;
    loadModels(qfn_root.ascii());
}

void QtWidgetMicrograph::loadModels(const FileName &fn)
{
    __modelRootName = fn;

    // Load parameters
    std::string dummy;
    std::ifstream fh_params;
    fh_params.open((__modelRootName + ".param").c_str());
    if (!fh_params) return;
    fh_params >> dummy >> __gray_bins
              >> dummy >> __radial_bins
              >> dummy >> __piece_xsize
              >> dummy >> __highpass_cutoff
              >> dummy >> __particle_radius
              >> dummy >> __min_distance_between_particles
              >> dummy >> __scan_overlap
              >> dummy >> __penalization
    ;
    fh_params.close();
    __piece_overlap=2*__particle_radius;
// COSS:    if (__particle_radius>100) __output_scale=1;
// COSS:   else __output_scale=0;
    __reduction=(int)std::pow(2.0, __output_scale);

    // Load the mask
    __mask.type = READ_MASK;
    __mask.fn_mask = __modelRootName + ".mask";
    __mask.generate_2Dmask();
    __mask.get_binary_mask2D().setXmippOrigin();
    __mask_size = XSIZE(__mask.get_binary_mask2D()) * __reduction;
    classifyMask();

    // Load training vectors
    std::ifstream fh_training;
    fh_training.open((__modelRootName + ".training").c_str());
    if (!fh_training)
        REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                     __modelRootName + ".training" + " for input");
    
    fh_training >> __training_loaded_model;
    fh_training.close();

    // If we have already learnt or loaded the model, clear the training
    if (__learn_particles_done || __is_model_loaded)
        __training_model.clear();

    // Particles have not been learnt but loaded from a file
    __learn_particles_done = false;
    __is_model_loaded = true;
    buildSelectionModel();
    std::cout << "The model has been loaded..." << std::endl;
}

/* Save particles ---------------------------------------------------------- */
void QtWidgetMicrograph::saveAutoParticles()
{
    if (__autoselection_done && __auto_label!=-1)
        __m->write_coordinates(__auto_label, __m->micrograph_name() +
                               ".auto.pos");
}

/* Save models ------------------------------------------------------------- */
void QtWidgetMicrograph::saveModels(bool askFilename)
{
    if (__particle_radius!=0) {
        // Get the rootname
        std::string fn_root;
        if (askFilename)
        {
            bool ok;
            QString qfn_root = QInputDialog::getText("Saving model",
                               "Model", QLineEdit::Normal,
                               "Model", &ok);
            if (!ok || qfn_root.isEmpty()) return;
            fn_root = qfn_root.ascii();
        }
        else
            fn_root=__modelRootName;

        // Save the automatically selected particles
        saveAutoParticles();

        // Save the mask
        ImageXmipp save;
        typeCast(__mask.get_binary_mask2D(), save());
        save.write(fn_root + ".mask");

        // Save parameters
        std::ofstream fh_params;
        fh_params.open((fn_root + ".param").c_str());
        if (!fh_params)
            REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                         fn_root + ".param" + " for output");
        fh_params << "gray_bins=                      " << __gray_bins                      << std::endl
                  << "radial_bins=                    " << __radial_bins                    << std::endl
                  << "piece_xsize=                    " << __piece_xsize                    << std::endl
                  << "highpass=                       " << __highpass_cutoff                << std::endl
                  << "particle_radius=                " << __particle_radius                << std::endl
                  << "min_distance_between_particles= " << __min_distance_between_particles << std::endl
                  << "particle_overlap=               " << __scan_overlap                   << std::endl
                  << "penalization=                   " << __penalization                   << std::endl
        ;
        fh_params.close();

        // Save training vectors
        Classification_model aux_model;
        aux_model = __training_model;
        std::ofstream fh_training;
        fh_training.open((fn_root + ".training").c_str());
        if (!fh_training)
            REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                         fn_root + ".training" + " for output");
        fh_training << aux_model << std::endl;
        fh_training.close();
        std::cout << "The model has been saved..." << std::endl;
    }
}

/* Configure auto ---------------------------------------------------------- */
void QtWidgetMicrograph::configure_auto()
{

    QDialog   setPropertiesDialog(this, 0, TRUE);
    setPropertiesDialog.setCaption("Configure AutoSelect");
#ifdef QT3_SUPPORT
    Q3Grid    qgrid(2, &setPropertiesDialog);
#else
    QGrid     qgrid(2, &setPropertiesDialog);
#endif
    qgrid.setMinimumSize(250, 100);

    QLabel    lparticle_radius("Particle radius: ", &qgrid);
    QLineEdit particle_radius(&qgrid);
    particle_radius.setText(integerToString(__particle_radius).c_str());

    QLabel    lpenalization("Penalization: ", &qgrid);
    QLineEdit penalization(&qgrid);
    penalization.setText(removeSpaces(floatToString(__penalization)).c_str());

    QPushButton okButton("Ok", &qgrid);
    QPushButton cancelButton("Cancel", &qgrid);

    connect(&okButton, SIGNAL(clicked(void)),
            &setPropertiesDialog, SLOT(accept(void)));
    connect(&cancelButton, SIGNAL(clicked(void)),
            &setPropertiesDialog, SLOT(reject(void)));

    if (setPropertiesDialog.exec())
    {
        __particle_radius = particle_radius.text().toInt();
        __penalization = penalization.text().toFloat();

        __gray_bins = 8;
        __radial_bins = 16;
        if (4*__particle_radius<512)
            __piece_xsize = 512;
        else
            __piece_xsize = NEXT_POWER_OF_2(6*__particle_radius);
        __highpass_cutoff = 0.02;
        __mask_size = 2*__particle_radius;
        __min_distance_between_particles = __particle_radius/2;
        __scan_overlap = ROUND(__mask_size*0.9);
        __learn_overlap = __particle_radius;
// COSS:        if (__particle_radius>100) __output_scale=1;
// COSS:        else __output_scale=0;
        __reduction=(int)std::pow(2.0, __output_scale);
    }
}

void QtWidgetMicrograph::changeContrast(int _mingray, int _maxgray, float _gamma)
{
    __mingray = _mingray;
    __maxgray = _maxgray;
    __gamma  = _gamma;
    __mImage->changeContrast(_mingray, _maxgray, _gamma);
    __mImageOverview->changeContrast(_mingray, _maxgray, _gamma);
}

void QtWidgetMicrograph::changeMarkType(int _type)
{
    __ellipse_type = _type;
    __mImage->__ellipse_type = __ellipse_type;
    __mImage->repaint(true);
}

void QtWidgetMicrograph::changeCircleRadius(float _circle_radius)
{
    __ellipse_radius = _circle_radius;
    __mImage->__ellipse_radius = __ellipse_radius;
    __mImage->repaint(true);
}

void QtWidgetMicrograph::repaint()
{
    __mImage->repaint(true);
    __mImageOverview->repaint(true);
}

void QtWidgetMicrograph::slotDrawEllipse(int _x, int _y, int _f)
{
    __mImage->drawEllipse(_x, _y, _f, __ellipse_radius, __ellipse_type);
    __mImageOverview->drawEllipse(_x, _y, _f);
}

void QtWidgetMicrograph::slotDrawLastEllipse(int _x, int _y, int _f)
{
    __mImage->drawEllipse(_x, _y, _f, __ellipse_radius, __ellipse_type);
    __mImageOverview->drawEllipse(_x, _y, _f);
    __mImage->drawLastEllipse(_x, _y, _f, __ellipse_radius, __ellipse_type);
}

/* Active family ----------------------------------------------------------- */
void QtWidgetMicrograph::slotActiveFamily(int _f)
{
    __activeFamily = _f;
    emit signalActiveFamily(_f);
}

void QtWidgetMicrograph::slotAddFamily(const char *_familyName)
{
    emit signalAddFamily(_familyName);
}

void QtWidgetMicrograph::slotDeleteMarkOther(int _coord)
{
    __m->coord(_coord).valid = false;
    emit signalRepaint();
}

void QtWidgetMicrograph::slotDeleteAutomatic(int _coord)
{
    __m->coord(_coord).valid = false;
    emit signalRepaint();
}

void QtWidgetMicrograph::slotQuit()
{
    __file_menu->slotQuit();
}

void QtWidgetMicrograph::slotChangeContrast()
{
    AdjustContrastWidget *adjustContrast = new
                                           AdjustContrastWidget(0, 255, 1.0F, this,
                                                                0, "new window", Qt::WDestructiveClose);
    adjustContrast->show();
}

void QtWidgetMicrograph::slotChangeCrop()
{
    CropWidget *crop = new CropWidget(this, 0, "new window", Qt::WDestructiveClose);
    connect(crop, SIGNAL(new_value(std::vector<int>)),
            __mImageOverview, SLOT(slotDrawCropArea(std::vector<int>)));
    crop->show();
}

void QtWidgetMicrograph::slotChangeCircleRadius()
{
    AdjustCircleRadiustWidget *adjustCircleRadius = new
            AdjustCircleRadiustWidget(0, 255, 10, this,
                                      0, "new window", Qt::WDestructiveClose);
    adjustCircleRadius->show();
}

/* ------------------------------------------------------------------------- */
/* Adjust contrast widget                                                    */
/* ------------------------------------------------------------------------- */
/* AdjustContrastWidget ---------------------------------------------------- */
// Constructor
AdjustContrastWidget::AdjustContrastWidget(int min, int max, float gamma,
        QtWidgetMicrograph *_qtwidgetmicrograph,
        QWidget *parent, const char *name, int wflags):
#ifdef QT3_SUPPORT
        QWidget(parent, name, (Qt::WindowFlags) wflags)
#else
        QWidget(parent, name, wflags)
#endif
{
    __qtwidgetmicrograph = _qtwidgetmicrograph;

    // Set this window caption
    setCaption("Adjust Contrast");

    // Create a layout to position the widgets
#ifdef QT3_SUPPORT
    Q3BoxLayout *Layout = new Q3VBoxLayout(this, 10);
#else
    QBoxLayout *Layout = new QVBoxLayout(this, 10);
#endif

    // Create a grid layout to hold most of the widgets
#ifdef QT3_SUPPORT
    Q3GridLayout *grid = new Q3GridLayout(3, 3);
#else
    QGridLayout *grid = new QGridLayout(3, 3);
#endif
    Layout->addLayout(grid, 5);

    // Minimum
    QLabel     *label_min = new QLabel(this, "label");
    label_min->setFont(QFont("times", 12, QFont::Bold));
    label_min->setText("Minimum");
    label_min->setFixedSize(label_min->sizeHint());
    grid->addWidget(label_min, 0, 0, Qt::AlignCenter);

    __scroll_min = new QScrollBar(0, 255, 1, 1, min,
                                  Qt::Horizontal, this, "scroll");
    __scroll_min->setFixedWidth(100);
    __scroll_min->setFixedHeight(15);
    grid->addWidget(__scroll_min, 0, 1, Qt::AlignCenter);
    connect(__scroll_min, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_min = new QLabel(this, "label");
    __label_min->setFont(QFont("courier", 14));
    __label_min->setText(integerToString(min, 3).c_str());
    __label_min->setFixedSize(__label_min->sizeHint());
    grid->addWidget(__label_min, 0, 2, Qt::AlignCenter);

    // Maximum
    QLabel     *label_max = new QLabel(this, "label");
    label_max->setFont(QFont("times", 12, QFont::Bold));
    label_max->setText("Maximum");
    label_max->setFixedSize(label_max->sizeHint());
    grid->addWidget(label_max, 1, 0, Qt::AlignCenter);

    __scroll_max = new QScrollBar(0, 255, 1, 1, max,
                                  Qt::Horizontal, this, "scroll");
    __scroll_max->setFixedWidth(100);
    __scroll_max->setFixedHeight(15);
    grid->addWidget(__scroll_max, 1, 1, Qt::AlignCenter);
    connect(__scroll_max, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_max = new QLabel(this, "label");
    __label_max->setFont(QFont("courier", 14));
    __label_max->setText(integerToString(max, 3).c_str());
    __label_max->setFixedSize(__label_max->sizeHint());
    grid->addWidget(__label_max, 1, 2, Qt::AlignCenter);

    // Gamma
    QLabel     *label_gamma = new QLabel(this, "label");
    label_gamma->setFont(QFont("times", 12, QFont::Bold));
    label_gamma->setText("Gamma");
    label_gamma->setFixedSize(label_gamma->sizeHint());
    grid->addWidget(label_gamma, 2, 0, Qt::AlignCenter);

    __scroll_gamma = new QScrollBar(0, 40, 1, 1, (int)(10*gamma),
                                    Qt::Horizontal, this, "scroll");
    __scroll_gamma->setFixedWidth(100);
    __scroll_gamma->setFixedHeight(15);
    grid->addWidget(__scroll_gamma, 2, 1, Qt::AlignCenter);
    connect(__scroll_gamma, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_gamma = new QLabel(this, "label");
    __label_gamma->setFont(QFont("courier", 14));
    __label_gamma->setText(floatToString(gamma, 3, 2).c_str());
    __label_gamma->setFixedSize(__label_gamma->sizeHint());
    grid->addWidget(__label_gamma, 2, 2, Qt::AlignCenter);

}

// One of the sliders changed ----------------------------------------------
void AdjustContrastWidget::scrollValueChanged(int new_val)
{
    __label_min  ->setText(integerToString(__scroll_min  ->value(), 3).c_str());
    __label_max  ->setText(integerToString(__scroll_max  ->value(), 3).c_str());
    __label_gamma->setText(floatToString((__scroll_gamma->value()) / 10.0, 3, 2).c_str());
    __qtwidgetmicrograph->changeContrast(__scroll_min->value(),
                                         __scroll_max->value(), __scroll_gamma->value() / 10.0);
}

/* ------------------------------------------------------------------------- */
/* Crop widget                                                               */
/* ------------------------------------------------------------------------- */
/* CropWidget -------------------------------------------------------------- */
// Constructor
CropWidget::CropWidget(QtWidgetMicrograph *_qtwidgetmicrograph,
                       QWidget *parent, const char *name, int wflags):
#ifdef QT3_SUPPORT
                       QWidget(parent, name, (Qt::WindowFlags) wflags)
#else
                       QWidget(parent, name, wflags)
#endif
{
    __qtwidgetmicrograph = _qtwidgetmicrograph;
    int Xdim, Ydim;
    __qtwidgetmicrograph->getMicrograph()->size(Xdim, Ydim);

    // Set this window caption
    setCaption("Crop micrograph");

    // Create a layout to position the widgets
#ifdef QT3_SUPPORT
    Q3BoxLayout *Layout = new Q3VBoxLayout(this, 10);
#else
    QBoxLayout *Layout = new QVBoxLayout(this, 10);
#endif

    // Create a grid layout to hold most of the widgets
#ifdef QT3_SUPPORT
    Q3GridLayout *grid = new Q3GridLayout(6, 3);
#else
    QGridLayout *grid = new QGridLayout(6, 3);
#endif
    Layout->addLayout(grid, 5);

    // Layout the four bars
    std::vector<int> min, max, init_value;
    std::vector<char *> prm_name;
    prm_name.push_back("x0");
    min.push_back(0);
    max.push_back(Xdim);
    init_value.push_back(ROUND(0.25*Xdim));
    prm_name.push_back("y0");
    min.push_back(0);
    max.push_back(Ydim);
    init_value.push_back(ROUND(0.25*Ydim));
    prm_name.push_back("xF");
    min.push_back(0);
    max.push_back(Xdim);
    init_value.push_back(ROUND(0.75*Xdim));
    prm_name.push_back("yF");
    min.push_back(0);
    max.push_back(Ydim);
    init_value.push_back(ROUND(0.75*Ydim));
    for (int i = 0; i < min.size(); i++)
    {

        // Add Parameter name
        QLabel     *lab1 = new QLabel(this, "lab1");
        lab1->setFont(QFont("times", 12, QFont::Bold));
        lab1->setText(prm_name[i]);
        lab1->setFixedSize(lab1->sizeHint());
        grid->addWidget(lab1, i, 0, Qt::AlignLeft);

        // Add Scroll Bar
        QScrollBar  *scroll_aux = new QScrollBar(min[i], max[i], 1, 50, (int)init_value[i], Qt::Horizontal, this, "scroll");
        scroll_aux->setFixedWidth(100);
        scroll_aux->setFixedHeight(15);
        grid->addWidget(scroll_aux, i, 1, Qt::AlignCenter);
        __scroll.push_back(scroll_aux);

        // Label for the current value
        QLabel * value_lab_aux;
        value_lab_aux = new QLabel(this, "value_lab");
        value_lab_aux->setFont(QFont("times", 12));
        value_lab_aux->setNum(init_value[i]);
        grid->addWidget(value_lab_aux, i, 2, Qt::AlignLeft);
        __label.push_back(value_lab_aux);

        connect(scroll_aux, SIGNAL(valueChanged(int)), SLOT(scrollValueChanged(int)));
    }

    // Layout the output name
    QLabel     *lab2 = new QLabel(this, "lab2");
    lab2->setFont(QFont("times", 12, QFont::Bold));
    lab2->setText("Output image");
    lab2->setFixedSize(lab2->sizeHint());
    grid->addWidget(lab2, min.size() + 1, 0, Qt::AlignLeft);
    __outputNameLineEdit = new QLineEdit(this, "output name");
    grid->addWidget(__outputNameLineEdit, min.size() + 1, 1, Qt::AlignLeft);

    // Cancel Button
    QPushButton *cancel;
    cancel = new QPushButton(this, "cancel");     // create button 1
    cancel->setFont(QFont("times", 12, QFont::Bold));
    cancel->setText("Cancel");
    cancel->setFixedSize(cancel->sizeHint());
    grid->addWidget(cancel, min.size() + 2, 0, Qt::AlignVCenter);
    connect(cancel, SIGNAL(clicked()), this, SLOT(cancel()));

    // OK button
    QPushButton *do_it;
    do_it = new QPushButton(this, "do_it");     // create button 3
    do_it->setFont(QFont("times", 12, QFont::Bold));
    do_it->setText("Ok");
    do_it->setFixedHeight(do_it->sizeHint().height());
    do_it->setFixedWidth(80);
    grid->addWidget(do_it, min.size() + 2, 2, Qt::AlignVCenter);
    connect(do_it, SIGNAL(clicked()), this, SLOT(accept()));

    __qtwidgetmicrograph->overview()->init_crop_area();
}

// Destructor --------------------------------------------------------------
CropWidget::~CropWidget()
{
    for (int i = 0; i < __label.size(); i++)
    {
        delete __label[i];
        delete __scroll[i];
    }
}

// One of the sliders changed ----------------------------------------------
void CropWidget::scrollValueChanged(int new_val)
{
    std::vector<int> value;
    // Get values
    for (int i = 0; i < __label.size(); i++)
    {
        int v = __scroll[i]->value();
        value.push_back(ROUND((float)v));
    }

    // Check value validity
    value[0] = XMIPP_MIN(value[0], value[2]);
    value[2] = XMIPP_MAX(value[0], value[2]);
    value[1] = XMIPP_MIN(value[1], value[3]);
    value[3] = XMIPP_MAX(value[1], value[3]);

    // Set these values
    for (int i = 0; i < __label.size(); i++)
    {
        __label[i]->setNum(value[i]);
        __scroll[i]->setValue(value[i]);
    }

    emit new_value(value);
}

void CropWidget::accept()
{
    __qtwidgetmicrograph->overview()->finish_crop_area();
    // Get values
    std::vector<int> value;
    for (int i = 0; i < __label.size(); i++)
        value.push_back(__scroll[i]->value());

    // Get output image
    std::string fn_out = __outputNameLineEdit->text().ascii();
    if (fn_out == "")
    {
        QMessageBox::information(this, "Mark",
                                 "The output image is empty\n Cropping is not carried out\n");
        close();
        return;
    }

    // Do the cropping
    int w = value[2] - value[0];
    int h = value[3] - value[1];
    std::string command = (std::string)"xmipp_window_micrograph " +
                     "-i " + __qtwidgetmicrograph->getMicrograph()->micrograph_name() +
                     " -o " + fn_out +
                     " -size " + integerToString(w, 0) + " " + integerToString(h, 0) +
                     " -top_left_corner " + integerToString(value[0], 0) + " " + integerToString(value[1], 0);
    std::cout << "Executing:\n" << command << std::endl;
    system(command.c_str());

    // Close the parameters window
    close();
}

void CropWidget::cancel()
{
    __qtwidgetmicrograph->overview()->finish_crop_area();
    close();
}

/* ------------------------------------------------------------------------- */
/* Adjust circle widget                                                      */
/* ------------------------------------------------------------------------- */
/* AdjustCircleWidget ------------------------------------------------------ */
// Constructor
AdjustCircleRadiustWidget::AdjustCircleRadiustWidget(int min, int max,
        int start_with, QtWidgetMicrograph *_qtwidgetmicrograph,
        QWidget *parent, const char *name, int wflags):
#ifdef QT3_SUPPORT
        QWidget(parent, name, (Qt::WindowFlags) wflags)
#else
        QWidget(parent, name, wflags)
#endif
{
    __qtwidgetmicrograph = _qtwidgetmicrograph;

    // Set this window caption
    setCaption("Change Circle Radius");

    // Create a layout to position the widgets
#ifdef QT3_SUPPORT
    Q3BoxLayout *Layout = new Q3VBoxLayout(this, 3);
#else
    QBoxLayout *Layout = new QVBoxLayout(this, 3);
#endif

    // Create a grid layout to hold most of the widgets
#ifdef QT3_SUPPORT
    Q3GridLayout *grid = new Q3GridLayout(1, 3);
#else
    QGridLayout *grid = new QGridLayout(1, 3);
#endif
    Layout->addLayout(grid);

    // Radius
    QLabel     *label_radius = new QLabel(this, "label");
    label_radius->setFont(QFont("times", 12, QFont::Bold));
    label_radius->setText("Radius");
    label_radius->setFixedSize(label_radius->sizeHint());
    grid->addWidget(label_radius, 0, 0, Qt::AlignCenter);

    __scroll_radius = new QScrollBar(0, 255, 1, 10, start_with,
                                     Qt::Horizontal, this, "scroll");
    __scroll_radius->setFixedWidth(100);
    __scroll_radius->setFixedHeight(15);
    grid->addWidget(__scroll_radius, 0, 1, Qt::AlignCenter);
    connect(__scroll_radius, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_radius = new QLabel(this, "label");
    __label_radius->setFont(QFont("courier", 14));
    __label_radius->setText(integerToString(start_with, 3).c_str());
    __label_radius->setFixedSize(__label_radius->sizeHint());
    grid->addWidget(__label_radius, 0, 2, Qt::AlignCenter);
}

void AdjustCircleRadiustWidget::scrollValueChanged(int new_val)
{
    __label_radius  ->setText(integerToString(__scroll_radius  ->value(), 3).c_str());
    __qtwidgetmicrograph->changeCircleRadius(__scroll_radius->value());
}
