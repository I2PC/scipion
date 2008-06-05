/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
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
#include <data/morphology.h>
#include <data/rotational_spectrum.h>
#include <data/denoise.h>
#include <reconstruction/fourier_filter.h>

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

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const Particle &_p)
{
    _out << _p.x      << " " << _p.y << " "
    << _p.idx    << " "
    << (int)_p.status << " "
    << _p.dist   << " "
    << _p.vec.transpose()
    << std::endl
    ;
    return _out;
}

void Particle::read(std::istream &_in, int _vec_size)
{
    _in >> x >> y
    >> idx
    >> status
    >> dist;
    status -= '0';
    vec.resize(_vec_size);
    _in >> vec;
}

/* Clear ------------------------------------------------------------------- */
void Classification_model::clear()
{
    __avg.clear();
    __sigma.clear();
    __sigma_inv.clear();
    __training_particle.clear();
    __largest_distance = -1;
}

/* Import particles -------------------------------------------------------- */
void Classification_model::import_particles(const Classification_model &_model)
{
    int imax = _model.__training_particle.size();
    for (int i = 0;i < imax;i++)
        if (_model.__training_particle.at(i).status == 1)
            add_particle(_model.__training_particle.at(i));
}

/* Build statistical model ------------------------------------------------- */
void Classification_model::build_model()
{
    int N = __training_particle.size();
    if (N == 0)
    {
        clear();
        return;
    }
    if (N == 1)
    {
        __avg = __training_particle.at(0).vec;
        __sigma.initIdentity(XSIZE(__avg));
        __sigma_inv = __sigma;
        __training_particle.at(0).dist = 0;
        return;
    }

    // Compute the average
    __avg.initZeros(XSIZE(__training_particle.at(0).vec));
    for (int i = 0; i < N; i++) __avg += __training_particle.at(i).vec;
    __avg /= N;

    // Compute the covariance
    __sigma.initZeros(XSIZE(__avg), XSIZE(__avg));
    for (int i = 0;i < N;i++)
    {
        Matrix2D<double> X;
        X.fromVector(__training_particle.at(i).vec - __avg);
        __sigma += (X * X.transpose());
    }
    __sigma /= N - 1;

    __well_posed = well_posed();
    if (__well_posed) __sigma_inv = __sigma.inv();
    else
    {
        // Keep only the diagonal
        std::cout << "The training model is not well posed.\n"
        << "A weak independent model used instead\n";
        __sigma_inv.initZeros(__sigma);
        for (int i = 0; i < XSIZE(__sigma); i++)
            __sigma_inv(i, i) = 1 / __sigma(i, i);
    }

    //compute the distance of each training vector from average
    for (int i = 0;i < N;i++)
    {
        __training_particle.at(i).dist = distance_to_average(
                                             __training_particle.at(i).vec);
    }
}

double Classification_model::distance(const Matrix1D<double> &X,
                                      const Matrix1D<double> &Y)
{
    if (XSIZE(__sigma_inv) == 0) return 1e30;
    Matrix2D<double> dif;
    dif.fromVector(X - Y);
    double dist2;
    dist2 = (dif.transpose() * __sigma_inv * dif)(0, 0);
    return sqrt(dist2);
}

double Classification_model::euclidean_distance_to_average(
    const Matrix1D<double> &my_X)
{
    if (XSIZE(__avg) == 0) return 1e30;
    Matrix1D<double> dif = my_X - __avg;
    return dif.sum2();
}

#define DEBUG
void Classification_model::compute_largest_distance()
{
    __largest_distance = 0.0;
    int N = __training_particle.size();
    Matrix1D<double> dist(N);
#ifdef DEBUG
    std::cout << "Computing largest distance ...\n";
#endif
    for (int i = 0;i < N;i++)
    {
        dist(i) = distance_to_average(__training_particle.at(i).vec);
#ifdef DEBUG
        std::cout << "   Distance of " << i << " = " << dist(i) << std::endl;
#endif
    }
    histogram1D hist;
    compute_hist(dist, hist, 50);
    __largest_distance = hist.percentil(90);

#ifdef DEBUG
    std::cout << "Largest distance = " << __largest_distance << std::endl;
#endif
}
#undef DEBUG

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const Classification_model &_m)
{
    int imax = _m.__training_particle.size();
    int N = 0;
    for (int i = 0; i < imax; i++)
        if (_m.__training_particle.at(i).status != 0) N++;
    _out << "#N.particles= " << N << std::endl;
    if (imax > 0)
        _out << "#Vector_size= " << XSIZE(_m.__training_particle.at(0).vec) << std::endl;
    else
        _out << "#Vector_size= 0\n";
    _out << "#x y index_in_micrograph status dist_to_average vector\n";
    for (int i = 0; i < imax; i++)
        if (_m.__training_particle.at(i).status != 0)
            _out << _m.__training_particle.at(i);
    return _out;
}

std::istream & operator >> (std::istream &_in, Classification_model &_m)
{
    std::string dummy;
    int imax, vec_size;
    _in >> dummy >> imax;
    _in >> dummy >> vec_size;
    for (int i = 0; i < 6; i++) _in >> dummy;
    for (int i = 0; i < imax; i++)
    {
        Particle *P = new Particle;
        P->read(_in, vec_size);
        _m.add_particle(*P);
    }

    return _in;
}

/* Print model ------------------------------------------------------------- */
void Classification_model::print_model(std::ostream &_out)
{
    _out << "Average of the model: " << __avg.transpose() << std::endl;
    _out << "Covariance of the model:\n" << __sigma << std::endl;
    _out << "Largest distance:" << __largest_distance << std::endl;
}

/* Constructor ------------------------------------------------------------- */
QtWidgetMicrograph::QtWidgetMicrograph(QtMainWidgetMark *_mainWidget,
                                       QtFiltersController *_f,
                                       Micrograph *_m) :
        QWidget((QWidget*) _mainWidget)
{
    __filtersController = _f;
    __m              = NULL;
    __activeFamily   = -1;
    __tilted         = FALSE;
    __learn_particles_done = FALSE;
    __autoselection_done = FALSE;
    __auto_label     = -1;
    __use_euclidean_distance_for_errors = true;
    __gray_bins = 8;
    __radial_bins = 24;
    __keep = 0.95;
    __piece_xsize = 512;
    __piece_ysize = 512;
    __output_scale = 1;
    __highpass_cutoff = 0.02;
    __reduction = (int)pow(2.0, __output_scale);
    __particle_radius = 110;
    __min_distance_between_particles = 2 * __particle_radius;
    __Nerror_models = 1;
    __mask_size = 4 * __particle_radius;
    __piece_overlap = ROUND(__mask_size * 2);
    __particle_overlap = (int)(__mask_size * 9.0 / 10.0);
    __numin = 1;
    __numax = 15;
    __use_background = false;

#ifdef QT3_SUPPORT
    __gridLayout     = new Q3VBoxLayout(this);
#else
    __gridLayout     = new QVBoxLayout(this);
#endif
    __menuBar        = new QMenuBar(this);
    __menuBar->setSeparator(QMenuBar::InWindowsStyle);
    __mImage         = new QtImageMicrograph(0);
    __mImage->setWidgetMicrograph(this);
    __mImageOverview = new QtImageOverviewMicrograph(this);
    __mImageOverview->setWidgetMicrograph(this);
    __file_menu      = NULL;
    __ellipse_radius = 5;
    __ellipse_type   = MARK_CIRCLE;

    __mImage->setFiltersController(_f);
    __mImageOverview->setFiltersController(_f);

    connect(__mImageOverview, SIGNAL(signalSetCoords(int, int)),
            __mImage, SLOT(slotSetCoords(int, int)));
    connect(__mImage, SIGNAL(signalSetCoords(int, int)),
            __mImageOverview, SLOT(slotSetCoords(int, int)));
    connect(__mImage, SIGNAL(signalSetWidthHeight(int, int)),
            __mImageOverview, SLOT(slotSetWidthHeight(int, int)));
    connect(__mImage, SIGNAL(signalRepaint(void)),
            __mImageOverview, SLOT(slotRepaint(void)));
    connect(__mImageOverview, SIGNAL(signalRepaint(void)),
            __mImage, SLOT(slotRepaint(void)));
    connect(__mImage, SIGNAL(signalRepaint(void)),
            __mImage, SLOT(slotRepaint(void)));
    connect(__mImageOverview, SIGNAL(signalRepaint(void)),
            __mImageOverview, SLOT(slotRepaint(void)));
    connect(__mImage, SIGNAL(signalAddCoordOther(int, int, int)),
            this, SLOT(slotDrawEllipse(int, int, int)));

    connect(this, SIGNAL(signalActiveFamily(int)),
            __mImage, SLOT(slotActiveFamily(int)));
    connect(this, SIGNAL(signalActiveFamily(int)),
            __mImageOverview, SLOT(slotActiveFamily(int)));

#ifdef QT3_SUPPORT
    Q3Accel *ctrl = new Q3Accel(this);
#else
    QAccel *ctrl = new QAccel(this);
#endif
    ctrl->connectItem(ctrl->insertItem(Qt::Key_G + Qt::CTRL, 200),
                      this, SLOT(slotChangeContrast(void)));
#ifdef QT3_SUPPORT
    Q3Accel *ctrl2 = new Q3Accel(this);
#else
    QAccel *ctrl2 = new QAccel(this);
#endif
    ctrl2->connectItem(ctrl2->insertItem(Qt::Key_R + Qt::CTRL, 200),
                       this, SLOT(slotChangeCircleRadius(void)));

    setMicrograph(_m);

    __mImage->show();
    __gridLayout->setMenuBar(__menuBar);
    __gridLayout->addWidget(__mImageOverview);

    openMenus();
}

QtWidgetMicrograph::~QtWidgetMicrograph()
{
    delete __mImage;
    delete __mImageOverview;
    delete __menuBar;
    delete __gridLayout;
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

/* Open menus -------------------------------------------------------------- */
void QtWidgetMicrograph::openMenus()
{
    __file_menu = new QtFileMenu(this);
    connect(__mImage, SIGNAL(signalAddCoordOther(int, int, int)),
            __file_menu, SLOT(slotCoordChange()));

    QtFilterMenu *filterMenu = new QtFilterMenu(this);

    QtAutoMenu *autoMenu = new QtAutoMenu(this);

    addMenuItem("&File", (QtPopupMenuMark *)(__file_menu));
    addMenuItem("F&ilters", (QtPopupMenuMark *)(filterMenu));
    // Sjors 31oct07: because this only gives core dumps, remove the
    // menu for now (just until it works)
    //addMenuItem("&AutoSelection", (QtPopupMenuMark *)(autoMenu));

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
            this, SLOT(slotRepaint()));

    // *** Add your own menus
}

/* Learn particles --------------------------------------------------------- */
void QtWidgetMicrograph::learnParticles()
{
    std::cerr << "\n------------------Learning Phase---------------------------\n";
    createMask();

    std::vector<int> all_idx;
    int num_part = __m->ParticleNo();
    for (int i = 0; i < num_part; i++)
        if (__m->coord(i).valid && __m->coord(i).label != __auto_label)
            all_idx.push_back(i);
    buildVectors(all_idx, __training_model);
    buildSelectionModel();
    __learn_particles_done = true;

    classify_errors();
}

void QtWidgetMicrograph::buildSelectionModel()
{
    __selection_model = __training_model;
    rebuild_moved_automatic_vectors();
    __selection_model.import_particles(__auto_model);
    __selection_model.import_particles(__training_loaded_model);
    __selection_model.import_particles(__auto_loaded_model);

    std::cerr << "Number of training particles :: "
    << __selection_model.__training_particle.size() << std::endl;

    __selection_model.build_model();
    __selection_model.compute_largest_distance();
}

/* Automatic phase --------------------------------------------------------- */
// This function used for sorting particles
bool sort_criteria(const Particle &p1, const Particle &p2)
{
    return (p1.dist < p2.dist);
}

//#define DEBUG
//#define DEBUG_MORE
//#define DEBUG_EVEN_MORE
void QtWidgetMicrograph::automaticallySelectParticles()
{
    if (XSIZE(__selection_model.__sigma_inv) == 0) return;
    std::cerr << "------------------Automatic Phase---------------------------" << std::endl;

    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    //get the threshold distance
    double threshold = __selection_model.__largest_distance;

    //define a vector to store automatically selected particles
    std::vector< Particle > candidate_vec, all_vec;

    //top,left corner of the piece
    int top = 0, left = 0, next_top = 0, next_left = 0;

    //If the piece available is small then include the scanned part
    //because we need bigger image for denoising but for scanning
    //particles we skip the already scanned part
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    Matrix1D<double> v;
    int N = 1;
    while (get_corner_piece(top, left, skip_y,
                            next_skip_x, next_skip_y, next_top, next_left))
    {
        std::cerr << "Processing piece number " << N << "...\n";
#ifdef DEBUG
        std::cerr << "    (top,left)=" << top << "," << left
        << " skip y,x=" << next_skip_y << "," << next_skip_x
        << " next=" << next_top << "," << next_left << std::endl;
#endif

        // Get a piece and prepare it
        if (!prepare_piece())
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

#ifdef DEBUG_MORE
        std::cerr << "Skip(y,x)=" << skip_y << "," << skip_x << std::endl;
#endif
        while (get_next_scanning_pos(next_posx, next_posy, skip_x, skip_y))
        {
#ifdef DEBUG_MORE
            std::cerr << "Pos(y,x)=" << posy << "," << posx
            << " Micro(y,x)=" << posy*__reduction + top
            << "," << posx*__reduction + left
            << " Next pos(y,x)=" << next_posy << "," << next_posx;
#endif
            bool success = build_vector(posx, posy, v);
            if (success)
            {
                double dist = __selection_model.distance_to_average(v);
#ifdef DEBUG_MORE
                std::cerr << " Success " << dist;
#endif

                // Build the Particle structure
                Particle P;
                P.x = left + posx * __reduction;
                P.y = top + posy * __reduction;
                P.idx = -1;
                P.status = 1;
                P.vec = v;
                P.dist = dist;

                // Insert it in the list of visited positions
                all_vec.push_back(P);

                // If it is likely to belong to the model
                if (dist < threshold)
                {
#ifdef DEBUG_MORE
                    std::cerr << "   Initial";
#endif
                    // Refine the particle position
                    P.x = posx;
                    P.y = posy;
                    refine_center(P);
                    P.x = left + P.x * __reduction;
                    P.y = top + P.y * __reduction;

                    // Insert it in the list of candidates
                    candidate_vec.push_back(P);
#ifdef DEBUG_EVEN_MORE
                    std::cout << "\n   " << P.vec.transpose() << std::endl;
#endif
                }
            }
#ifdef DEBUG_MORE
            std::cerr << std::endl;
#endif

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
        //break; //Uncomment for a single patch
    }
#ifdef DEBUG
    //int iall_scanned=add_family(all_vec,"all_scanned");
    //__m->write_coordinates(iall_scanned, __m->micrograph_name()+
    //   ".all_scanned.pos");
    int iinitial    = add_family(candidate_vec, "initial");
    __m->write_coordinates(iinitial, __m->micrograph_name() +
                           ".initial.pos");
#endif

    // Reject manually selected particles
    reject_previously_selected(__training_model, candidate_vec);
    reject_previously_selected(__auto_model, candidate_vec);

    //sort the particles in order of increasing distance
    sort(candidate_vec.begin(), candidate_vec.end(), sort_criteria);

    // reject the candidates that are pointing to the same particle
    int Nalive = reject_within_distance(candidate_vec, __particle_radius, false);
    std::cout << "Nalive=" << Nalive << std::endl;

    // reject the candidates that are two close to each other
    Nalive = reject_within_distance(candidate_vec, __min_distance_between_particles, true);
    std::cout << "Nalive=" << Nalive << std::endl;

    //insert selected particles in the result
    int imax = candidate_vec.size();
    int Nmax = FLOOR(Nalive * __keep);
    int n = 0;
    for (int i = 0;i < imax && n < Nmax; i++)
        if (candidate_vec.at(i).status == 1)
        {
            // Get the distance of the vector to the training set
            double dist_to_training;
            if (__use_euclidean_distance_for_errors)
                dist_to_training = __selection_model.euclidean_distance_to_average(
                                       candidate_vec.at(i).vec);
            else
                dist_to_training = candidate_vec.at(i).dist;

            bool is_error = false;
            for (int j = 0; j < __Nerror_models && !is_error; j++)
            {
                double dist_to_error_j;
                if (__use_euclidean_distance_for_errors)
                    dist_to_error_j = __error_model.at(j).euclidean_distance_to_average(
                                          candidate_vec.at(i).vec);
                else
                    dist_to_error_j = __error_model.at(j).distance_to_average(
                                          candidate_vec.at(i).vec);
                if (dist_to_error_j < dist_to_training) is_error = true;
            }

            if (!is_error)
            {
                __auto_model.add_particle(candidate_vec.at(i));
                n++;
            }
            else
            {
                candidate_vec.at(i).status = 0;
                std::cout << "Error found:" << candidate_vec.at(i).x << ","
                << candidate_vec.at(i).y << std::endl;
            }
        }

    std::cerr << "Number of automatically selected particles = " << n << std::endl;

    __auto_label = add_family(__auto_model.__training_particle, "auto");

    __autoselection_done = true;
}
#undef DEBUG
#undef DEBUG_MORE
#undef DEBUG_EVEN_MORE

/* Add family -------------------------------------------------------------- */
int QtWidgetMicrograph::add_family(std::vector<Particle> &_list,
                                   const std::string &_family_name)
{
    int ilabel = __m->add_label(_family_name);
    int imax = _list.size();
    for (int i = 0;i < imax;i++)
    {
        int idx = __m->add_coord(_list.at(i).x, _list.at(i).y, ilabel);
        _list.at(i).idx = idx;
    }
    return ilabel;
}

/* Create mask for learning particles -------------------------------------- */
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

/* Classify the mask pixels ------------------------------------------------ */
//#define DEBUG
void QtWidgetMicrograph::classifyMask()
{
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();
    if (XSIZE(mask) == 0) return;
    double max_radius_particle = __particle_radius / __reduction;

    // Determine max_radius
    double max_radius;
    if (__use_background) max_radius =
            sqrt((double)((XSIZE(mask) / 2) * (XSIZE(mask) / 2) +
                          (YSIZE(mask) / 2) * (YSIZE(mask) / 2)));
    else max_radius = max_radius_particle;

    // Initialize some variables
    // 6 is the minimum radius to be informative
    double radial_step = (max_radius - 6) / __radial_bins;

    Matrix2D<int> *classif1 = new Matrix2D<int>;
    classif1->resize(mask);
    classif1->initConstant(-1);
    Matrix1D<int> Nrad(__radial_bins);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
    {
        double radius = sqrt((double)(i * i + j * j));
        double angle = atan2((double)i, (double)j);
        if (angle < 0) angle += 2 * PI;

        // Classif1 is the classification for the radial mass distribution
        if (radius < max_radius)
        {
            int radius_idx;
            if (radius > 6) radius_idx = XMIPP_MIN(__radial_bins - 1, 1 + FLOOR((radius - 6) / radial_step));
            else          radius_idx = 0;
            (*classif1)(i, j) = radius_idx;
            Nrad(radius_idx)++;
        }
    }
    __mask_classification.push_back(classif1);

    // Create the holders for the radius values in classif1
    for (int i = 0; i < __radial_bins; i++)
    {
        Matrix1D<int> *aux = new Matrix1D<int>;
        aux->initZeros(Nrad(i));
        __radial_val.push_back(aux);
    }

#ifdef DEBUG
    ImageXmipp save;
    type_cast(*classif1, save());
    save.write("PPPclassif1.xmp");
#endif
}
#undef DEBUG

/* Build training vectors -------------------------------------------------- */
void QtWidgetMicrograph::buildVectors(std::vector<int> &_idx,
                                      Classification_model &__model)
{
    __model.clear();
    int num_part = _idx.size();
    __model.reserve_examples(num_part);
    Matrix1D<double> v;

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
        get_centered_piece(x, y, posx, posy);

        // Denoise, reduce, reject outliers and equalize histogram
        bool success = prepare_piece();
        if (!success) continue;
        posx = ROUND(posx / __reduction);
        posy = ROUND(posy / __reduction);

        //make vector from this particle
        success = build_vector(posx, posy, v);
        if (success)
        {
            Particle p;
            p.x = x;
            p.y = y;
            p.idx = part_idx;
            p.vec = v;
            p.status = 1;
            p.dist = 0.0;
            __model.add_particle(p);
        }

        //make vector from the neighbours
        std::vector< Matrix1D<int> > nbr;
        nbr.reserve(num_part);
        find_nbr(_idx, part_i, x, y, posx, posy, visited, nbr);
        for (int i = 0;i < nbr.size();i++)
        {
            part_i = nbr.at(i)(0);
            part_idx = _idx.at(part_i);
            posx = nbr.at(i)(1);
            posy = nbr.at(i)(2);
            success = build_vector(posx, posy, v);
            visited(part_i) = 1;
            if (success)
            {
                Particle p;
                p.x = __m->coord(part_idx).X;
                p.y = __m->coord(part_idx).Y;
                p.idx = part_idx;
                p.vec = v;
                p.status = 1;
                p.dist = 0.0;
                __model.add_particle(p);
            }
        }
    }
}

/* Build classification vector --------------------------------------------- */
//#define DEBUG
bool QtWidgetMicrograph::build_vector(int _x, int _y,
                                      Matrix1D<double> &_result)
{
    // First part is the foreground histogram
    // Second part is the background histogram
    // Third part is the radial mass distribution
    // The input image is supposed to be between 0 and bins-1
    _result.initZeros(__radial_bins*(__gray_bins - 1) + (__numax - __numin + 1));
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    if (STARTINGX(mask) + _x < STARTINGX(__piece)) return false;
    if (STARTINGY(mask) + _y < STARTINGY(__piece)) return false;
    if (FINISHINGX(mask) + _x > FINISHINGX(__piece)) return false;
    if (FINISHINGY(mask) + _y > FINISHINGY(__piece)) return false;

#ifdef DEBUG
    bool debug_go = false;
    ImageXmipp save, savefg, savebig;
    if (true)
    {
        save() = __piece;
        save.write("PPP0.xmp");
        std::cout << "Particle is at (y,x)=" << _y << "," << _x << std::endl;
        save().initZeros(YSIZE(mask), XSIZE(mask));
        STARTINGY(save()) = STARTINGY(mask);
        STARTINGX(save()) = STARTINGX(mask);
        savefg() = save();
        savefg().init_constant(-1);
        savebig() = savefg();
        debug_go = true;
    }
#endif

    Matrix1D<int> radial_idx(__radial_bins);
    Matrix2D<double> particle;
    particle.initZeros(YSIZE(mask), XSIZE(mask));
    STARTINGY(particle) = STARTINGY(mask);
    STARTINGX(particle) = STARTINGX(mask);

    FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
    {
        int val = (int)__piece(_y + i, _x + j);
        bool foreground = mask(i, j);

        // Classif 1 -> Histogram of the radial bins
        int idx0 = (*__mask_classification[0])(i, j);
        if (idx0 != -1)
            (*__radial_val[idx0])(radial_idx(idx0)++) = val;

        // Get particle
        if (foreground) particle(i, j) = val;

#ifdef DEBUG
        if (debug_go)
        {
            save(i, j) = val;
            if (foreground) savefg(i, j) = val;
        }
#endif
    }

    // Compute the histogram of the radial bins and store them
    int idx_result = 0;
    for (int i = 0; i < __radial_bins; i++)
    {
        histogram1D hist;
        compute_hist(*__radial_val[i], hist, 0, __gray_bins - 1, __gray_bins);
        for (int j = 0; j < __gray_bins - 1; j++)
            _result(idx_result++) = hist(j);
    }

    // Compute the rotational spectrum
    int dr = 1;
    Rotational_Spectrum spt;
    spt.rl = 3;
    spt.rh = __particle_radius / __reduction;
    spt.dr = dr;
    spt.numin = __numin;
    spt.numax = __numax;
    spt.x0 = (double)XSIZE(particle) / 2;
    spt.y0 = (double)YSIZE(particle) / 2;
    spt.compute_rotational_spectrum(particle, spt.rl, spt.rh, dr, spt.rh - spt.rl);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(spt.rot_spectrum)
    _result(idx_result++) = spt.rot_spectrum(i);

#ifdef DEBUG
    if (debug_go)
    {
        save.write("PPP1.xmp");
        savefg.write("PPP2.xmp");
        std::cout << _result.transpose() << std::endl;
        std::cout << "Distance=" << __selection_model.distance_to_average(_result);
        std::cout << "Press any key\n";
        char c;
        cin >> c;
    }
#endif
    return true;
}
#undef DEBUG

/* Get piece --------------------------------------------------------------- */
// to get the piece containing (x,y) of size xsize,ysize
// return the position of x,y in the piece in posx,posy
void QtWidgetMicrograph::get_centered_piece(int _x, int _y,
        int &_posx, int &_posy)
{
    __piece.resize(__piece_ysize, __piece_xsize);

    int startx = _x - ROUND(__piece_xsize / 2);
    int endx  = _x + ROUND(__piece_xsize / 2);
    int starty = _y - ROUND(__piece_ysize / 2);
    int endy  = _y + ROUND(__piece_ysize / 2);
    int maxx, maxy;
    __m->size(maxx, maxy);
    _posx = ROUND(__piece_xsize / 2);
    _posy = ROUND(__piece_ysize / 2);

    // boundry adjustments
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
        endy = __piece_ysize - 1;
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
        starty = endy - __piece_ysize;
    }

    //read the matrix from the micrograph
    for (int i = 0;i < __piece_ysize;i++)
        for (int j = 0;j < __piece_xsize;j++)
            __piece(i, j) = (*__m)(startx + j, starty + i);
}

// Get a piece whose top-left corner is at the desired position (if possible)
bool QtWidgetMicrograph::get_corner_piece(
    int _top, int _left, int _skip_y,
    int &_next_skip_x, int &_next_skip_y, int &_next_top, int &_next_left)
{
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();
    int maxx, maxy;
    __m->size(maxx, maxy);

    if (maxx < _left + __piece_xsize || maxy < _top + __piece_ysize) return false;

#ifdef NEVER_DEFINED
    _skip_x = _skip_y = 0;
    _next_left = _left + __piece_xsize - __piece_overlap;
    if (_next_left + __piece_xsize >= maxx)
    {
        if (maxx - _next_left < XSIZE(mask)*__reduction)
        {
            _next_left = 0;
            _next_top = _top + __piece_ysize - __piece_overlap;
            _skip_x = _skip_y = 0;
        }
        else
        {
            _skip_x = _next_left + __piece_xsize - maxx;
            _next_left -= _skip_x;
            if (_next_left != _left)_next_top = _top;
            else
            {
                _next_left = 0;
                _next_top = _top + __piece_ysize - __piece_overlap;
                _skip_x = _skip_y = 0;
            }
        }
    }
    else
    {
        _next_top = _top;
    }

    if ((_next_top + __piece_ysize >= maxy) &&
        (maxy - _next_top >= YSIZE(mask)*__reduction))
    {
        _skip_y = _next_top + __piece_ysize - maxy;
        _next_top -= _skip_y;
    }
#endif

    _next_skip_x = _next_skip_y = 0;
    bool increase_Y = false;
    if (_left + __piece_xsize != maxx)
    {
        _next_left = _left + __piece_xsize - __piece_overlap;
        if (_next_left + __piece_xsize >= maxx)
        {
            _next_left = maxx - __piece_xsize;
        }
    }
    else
    {
        _next_left = 0;
        increase_Y = true;
    }

    if (increase_Y)
    {
        if (_top + __piece_ysize != maxy)
        {
            _next_top = _top + __piece_ysize - __piece_overlap;
            if (_next_top + __piece_ysize >= maxy)
                _next_top = maxy - __piece_ysize;
        }
        else
        {
            _next_top = maxy;
        }
    }

    /* COSS: still To fix */
    /*
    if (_next_left+__piece_xsize==maxx)
       _next_skip_x=_left+__piece_xsize-_next_left-XSIZE(mask)/2;
    if (_next_top +__piece_ysize==maxy)
       if (_top==_next_top) _next_skip_y=_skip_y;
       else _next_skip_y=_top +__piece_ysize-_next_top-YSIZE(mask)/2;
    */

    //read the matrix from the micrograph
    __piece.resize(__piece_ysize, __piece_xsize);
    for (int i = 0;i < __piece_ysize;i++)
        for (int j = 0;j < __piece_xsize;j++)
            __piece(i, j) = (*__m)(_left + j, _top + i);
    return true;
}


/* Prepare piece ----------------------------------------------------------- */
bool QtWidgetMicrograph::prepare_piece()
{
    // High pass filter
    FourierMask Filter;
    Filter.FilterShape = RAISED_COSINE;
    Filter.FilterBand = HIGHPASS;
    Filter.w1 = __highpass_cutoff;
    Filter.raised_w = XMIPP_MIN(0.02, __highpass_cutoff);
    __piece.setXmippOrigin();
    Filter.generate_mask(__piece);
    Filter.apply_mask_Space(__piece);
    STARTINGX(__piece) = STARTINGY(__piece) = 0;

    // Denoise the piece
    Denoising_parameters denoiser;
    denoiser.denoising_type = Denoising_parameters::BAYESIAN;
    denoiser.scale = __output_scale + 3;
    denoiser.output_scale = __output_scale;
    denoiser.produce_side_info();
    denoiser.denoise(__piece);
    if (!(__piece(0, 0) == __piece(0, 0))) return false;

    // Reject 5% of the outliers
    reject_outliers(__piece, 5.0);

    // Equalize histogram
    histogram_equalization(__piece, __gray_bins);

    return true;
}

/* Get neighbours ---------------------------------------------------------- */
//To get the neighbours and their positions in the piece image
void QtWidgetMicrograph::find_nbr(std::vector<int> &_idx, int _index, int _x, int _y,
                                  int _posx, int _posy, Matrix1D<char> &_visited,
                                  std::vector< Matrix1D<int> > &_nbr)
{
    int piece_xsize = XSIZE(__piece);
    int piece_ysize = YSIZE(__piece);
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();

    //if all the particles are visited
    if (_visited.sum() == XSIZE(_visited)) return;

    int current_part = _index + 1;
    int current_part_idx = _idx.at(current_part);
    _nbr.clear();
    _nbr.reserve(XSIZE(_visited));

    int top = CEIL((double)_y / __reduction) - _posy;
    int left = CEIL((double)_x / __reduction) - _posx;
    int bottom = top + FLOOR((double)piece_ysize / __reduction);
    int right = left + FLOOR((double)piece_xsize / __reduction);
    int xmask2 = CEIL(XSIZE(mask) / 2);
    int ymask2 = CEIL(YSIZE(mask) / 2);

    while (current_part < XSIZE(_visited))
    {
        //find the next unvisited particle
        while (_visited(current_part)) current_part++;
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
    int &_x, int &_y, int _skip_x, int _skip_y)
{
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();
    if (_x + XSIZE(mask) / 2 >= XSIZE(__piece) ||
        _y + YSIZE(mask) / 2 >= YSIZE(__piece)) return false;

    if (_x == 0 && _y == 0)
    {
        _x = _skip_x + XSIZE(mask) / 2;
        _y = _skip_y + YSIZE(mask) / 2;
    }
    else
    {
        int nextx = _x + XSIZE(mask) - __particle_overlap / __reduction; // COSS: he a�adido __reduction
        int nexty = _y + YSIZE(mask) - __particle_overlap / __reduction; // COSS: he a�adido __reduction
        if (nextx + (XSIZE(mask) / 2) >= XSIZE(__piece))
        {
            if (nexty + (YSIZE(mask) / 2) >= YSIZE(__piece))
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
double dist_euc(const Particle &p1, const Particle &p2)
{
    return sqrt((double)(p1.x -p2.x)*(p1.x - p2.x) + (double)(p1.y - p2.y)*(p1.y - p2.y));
}

int QtWidgetMicrograph::reject_within_distance(
    std::vector<Particle> &_Input, double _min_dist,
    bool _reject_both)
{
    int imax = _Input.size();
    int n = 0;
    for (int i = 0;i < imax;i++)
    {
        if (_Input.at(i).status == 0) continue;
        for (int j = i + 1;j < imax;j++)
        {
            if (_Input.at(j).status == 0) continue;
            double dist = dist_euc(_Input.at(i), _Input.at(j));
            if (dist < _min_dist)
            {
                _Input.at(j).status = 0;
                if (_reject_both) _Input.at(i).status = 0;
            }
        }
        if (_Input.at(i).status == 1) n++;
    }
    return n;
}

/* Reject manually selected ------------------------------------------------ */
void QtWidgetMicrograph::reject_previously_selected(
    const Classification_model &_model,
    std::vector<Particle> &_candidate_vec)
{
    int imax = _candidate_vec.size();
    int jmax = _model.__training_particle.size();
    if (jmax == 0) return;
    for (int i = 0;i < imax;i++)
        for (int j = 0; j < jmax; j++)
            if (dist_euc(_candidate_vec.at(i), _model.__training_particle.at(j))
                < __min_distance_between_particles)
            {
                _candidate_vec.at(i).status = 0;
                break;
            }
}

/* Refine center of a particle --------------------------------------------- */
void QtWidgetMicrograph::refine_center(Particle &my_P)
{
    int previous_direction = -1;
    double dist_at_previous = 0, dist_left, dist_right, dist_top, dist_bottom;
    bool improvement = true;
    Matrix1D<double> current_vec;
    double current_dist = my_P.dist;
    int    current_x = my_P.x, current_y = my_P.y;
    const Matrix2D<int> &mask = __mask.get_binary_mask2D();
    bool success;
    do
    {
        double best_dist = current_dist;
        int    source   = -1;
        // Check at left
        int xpos = current_x - 1;
        if (previous_direction != 0 && xpos >= XSIZE(mask) / 2)
        {
            success = build_vector(xpos, current_y, current_vec);
            if (success)
            {
                dist_left = __selection_model.distance_to_average(current_vec);
                if (dist_left < best_dist)
                {
                    best_dist = dist_left;
                    source = 0;
                }
            }
        }

        // Check at right
        xpos = current_x + 1;
        if (previous_direction != 1 && xpos <= XSIZE(__piece) - XSIZE(mask) / 2)
        {
            success = build_vector(xpos, current_y, current_vec);
            if (success)
            {
                dist_right = __selection_model.distance_to_average(current_vec);
                if (dist_right < best_dist)
                {
                    best_dist = dist_right;
                    source = 1;
                }
            }
        }

        // Check at top
        int ypos = current_y - 1;
        if (previous_direction != 2 && ypos >= YSIZE(mask) / 2)
        {
            success = build_vector(current_x, ypos, current_vec);
            if (success)
            {
                dist_top = __selection_model.distance_to_average(current_vec);
                if (dist_top < best_dist)
                {
                    best_dist = dist_top;
                    source = 2;
                }
            }
        }

        // Check at bottom
        ypos = current_y + 1;
        if (previous_direction != 3 && ypos <= YSIZE(__piece) - YSIZE(mask) / 2)
        {
            success = build_vector(current_x, ypos, current_vec);
            if (success)
            {
                dist_bottom = __selection_model.distance_to_average(current_vec);
                if (dist_bottom < best_dist)
                {
                    best_dist = dist_bottom;
                    source = 3;
                }
            }
        }

        // Check which is best
        if (source == -1)
        {
            improvement = false;
        }
        else
        {
            switch (source)
            {
            case 0:
                current_x--;
                previous_direction = 1;
                break;
            case 1:
                current_x++;
                previous_direction = 0;
                break;
            case 2:
                current_y--;
                previous_direction = 3;
                break;
            case 3:
                current_y++;
                previous_direction = 2;
                break;
            }
            current_dist = best_dist;
        }
    }
    while (improvement);
    my_P.x = current_x;
    my_P.y = current_y;
    my_P.dist = current_dist;
}

/* Correct particles ------------------------------------------------------- */
void QtWidgetMicrograph::move_particle(int _idx)
{
    if (!__autoselection_done) return;
    int imax = __auto_model.__training_particle.size();
    for (int i = 0; i < imax; i++)
        if (__auto_model.__training_particle.at(i).idx == _idx)
        {
            __auto_model.__training_particle.at(i).status = 2;
            __auto_model.__training_particle.at(i).x = __m->coord(_idx).X;
            __auto_model.__training_particle.at(i).y = __m->coord(_idx).Y;
            break;
        }
}

void QtWidgetMicrograph::delete_particle(int _idx)
{
    if (!__autoselection_done) return;
    int imax = __auto_model.__training_particle.size();
    for (int i = 0; i < imax; i++)
        if (__auto_model.__training_particle.at(i).idx == _idx)
        {
            __auto_model.__training_particle.at(i).status = 0;
            __error_index.push_back(i);
            break;
        }
}

/* Load models ------------------------------------------------------------- */
void QtWidgetMicrograph::loadModels()
{
    // Get the rootname
    bool ok;
    QString qfn_root = QInputDialog::getText("Loading model",
                       "Model filename root", QLineEdit::Normal,
                       __m->micrograph_name().c_str(), &ok);
    if (!ok || qfn_root.isEmpty()) return;
    std::string fn_root = qfn_root.ascii();

    // Load parameters
    std::string dummy;
    std::ifstream fh_params;
    fh_params.open((fn_root + ".param").c_str());
    if (!fh_params)
    {
        std::cerr << (std::string)"QtWidgetMicrograph::write: Cannot open file " +
        fn_root + ".param for input" << std::endl;
        return;
    }
    fh_params >> dummy >> __gray_bins
              >> dummy >> __radial_bins
              >> dummy >> __keep
              >> dummy >> __piece_xsize
              >> dummy >> __piece_ysize
              >> dummy >> __particle_radius
              >> dummy >> __min_distance_between_particles
              >> dummy >> __output_scale
              >> dummy >> __reduction
              >> dummy >> __piece_overlap
              >> dummy >> __particle_overlap
              >> dummy >> __numin
              >> dummy >> __numax
              >> dummy >> __Nerror_models
    ;
    fh_params.close();

    // Load the mask
    __mask.type = READ_MASK;
    __mask.fn_mask = fn_root + ".mask";
    __mask.generate_2Dmask();
    __mask.get_binary_mask2D().setXmippOrigin();
    __mask_size = XSIZE(__mask.get_binary_mask2D()) * __reduction;
    classifyMask();

    // Load training vectors
    std::ifstream fh_training;
    fh_training.open((fn_root + ".training").c_str());
    if (!fh_training)
        REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                     fn_root + ".training" + " for input");
    fh_training >> __training_loaded_model;
    fh_training.close();

    // Load auto vectors
    std::ifstream fh_auto;
    fh_auto.open((fn_root + ".auto").c_str());
    if (!fh_auto)
        REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                     fn_root + ".auto" + " for input");
    fh_auto >> __auto_loaded_model;
    fh_auto.close();

    // Build the selection model
    buildSelectionModel();

    // Load error vectors
    for (int i = 0; i < __Nerror_models; i++)
    {
        Classification_model error;
        __error_model.push_back(error);
    }

    std::ifstream fh_error;
    __use_euclidean_distance_for_errors = false;
    for (int i = 0; i < __Nerror_models; i++)
    {
        fh_error.open((fn_root + ".error" + integerToString(i, 1)).c_str());
        if (!fh_error)
            REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                         fn_root + ".error" + integerToString(i, 1) + " for input");
        fh_error >> __error_model.at(i);
        fh_error.close();
        __error_model.at(i).build_model();
        __use_euclidean_distance_for_errors |= !__error_model.at(i).well_posed();
    }

    // Particles have not been learnt but loaded from a file
    __learn_particles_done = false;
}

/* Save models ------------------------------------------------------------- */
void QtWidgetMicrograph::saveModels()
{
    if (__autoselection_done)
        __m->write_coordinates(__auto_label, __m->micrograph_name() +
                               ".auto.pos");

    // Rebuild automatic vectors that have been moved
    rebuild_moved_automatic_vectors();

    // Classify errors
    classify_errors();

    // Write results to a file
    write();
}

/* Rebuild automatic vectors that have been moved -------------------------- */
void QtWidgetMicrograph::rebuild_moved_automatic_vectors()
{
    // Rebuild vectors
    std::vector<int> indexes_to_rebuild, indexes_to_rebuild_in_micrograph;
    int imax = __auto_model.__training_particle.size();
    for (int i = 0; i < imax; i++)
        if (__auto_model.__training_particle.at(i).status == 2)
        {
            indexes_to_rebuild.push_back(i);
            indexes_to_rebuild_in_micrograph.push_back(__auto_model.__training_particle.at(i).idx);
        }
    Classification_model rebuilt_vectors;
    buildVectors(indexes_to_rebuild_in_micrograph, rebuilt_vectors);

    // Update the vectors in auto_vec
    imax = indexes_to_rebuild.size();
    for (int i = 0; i < imax; i++)
    {
        int idx = indexes_to_rebuild.at(i);
        __auto_model.__training_particle.at(idx).status = 1;
        __auto_model.__training_particle.at(idx).vec =
            rebuilt_vectors.__training_particle.at(i).vec;
    }
}

/* Error classification ---------------------------------------------------- */
void QtWidgetMicrograph::classify_errors()
{
    // Initialize error models
    if (__error_model.size() == 0)
    {
        // Create the error models
        for (int i = 0; i < __Nerror_models; i++)
        {
            Classification_model error;
            __error_model.push_back(error);
        }
    }

    // Assign randomly the errors to the models
    int imax = __error_index.size();
    for (int i = 0; i < imax; i++)
        __error_model.at(i % __Nerror_models).
        add_particle(__auto_model.__training_particle.at(__error_index.at(i)));

    // Set all the particles to an active status within the error models
    for (int i = 0; i < __Nerror_models; i++)
    {
        int jmax = __error_model.at(i).__training_particle.size();
        for (int j = 0; j < jmax; j++)
            __error_model.at(i).__training_particle.at(j).status = 1;
    }

    // Iterate until no change
    bool change;
    do
    {
        change = false;

        // Compute the average and covariance of each class
        __use_euclidean_distance_for_errors = false;
        for (int i = 0; i < __Nerror_models; i++)
        {
            __error_model.at(i).build_model();
            __use_euclidean_distance_for_errors |= !__error_model.at(i).well_posed();
        }

        // Reassign each of the vectors in each class
        for (int i = 0; i < __Nerror_models && !change; i++)
        {
            Classification_model &error_class = __error_model.at(i);
            int lmax = error_class.__training_particle.size();
            std::vector<Particle>::iterator ptr = error_class.__training_particle.begin();
            for (int l = 0; l < lmax; l++, ptr++)
            {
                double best_dist = 0;
                Matrix1D<double> &current_vec =
                    error_class.__training_particle.at(l).vec;
                int source = -1;
                for (int j = 0; j < __Nerror_models; j++)
                {
                    double dist;
                    if (__use_euclidean_distance_for_errors)
                        dist = __error_model.at(j).euclidean_distance_to_average(current_vec);
                    else
                        dist = __error_model.at(j).distance_to_average(current_vec);
                    if (dist < best_dist || source == -1)
                    {
                        source = j;
                        best_dist = dist;
                    }
                }
                error_class.__training_particle.at(l).dist = best_dist;
                if (source != i)
                {
                    change = true;
                    __error_model.at(source).add_particle(
                        error_class.__training_particle.at(l));
                    error_class.__training_particle.erase(ptr);
                    break;
                }
            }
        }
    }
    while (change);
    __error_index.clear();
}

/* Write to a file --------------------------------------------------------- */
void QtWidgetMicrograph::write()
{
    // Get the rootname
    bool ok;
    std::string fn_root = (QInputDialog::getText("Saving model",
                                            "Model filename root", QLineEdit::Normal, __m->micrograph_name().c_str(), &ok)).ascii();
    if (!ok) return;

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
    fh_params << "gray_bins=                    " << __gray_bins                     << std::endl
    << "radial_bins=                    " << __radial_bins                   << std::endl
    << "keep=                         " << __keep                          << std::endl
    << "piece_xsize=                  " << __piece_xsize                   << std::endl
    << "piece_ysize=                  " << __piece_ysize                   << std::endl
    << "particle_radius=                " << __particle_radius               << std::endl
    << "min_distance_between_particles= " << __min_distance_between_particles << std::endl
    << "output_scale=                   " << __output_scale                   << std::endl
    << "reduction_factor=               " << __reduction                      << std::endl
    << "piece_overlap=                  " << __piece_overlap                  << std::endl
    << "particle_overlap=               " << __particle_overlap               << std::endl
    << "numin=                          " << __numin                          << std::endl
    << "numax=                          " << __numax                          << std::endl
    << "Nerror_models=                  " << __Nerror_models                  << std::endl
    ;
    fh_params.close();

    // Save training vectors
    Classification_model aux_model;
    aux_model = __training_model;
    aux_model.import_particles(__training_loaded_model);
    std::ofstream fh_training;
    fh_training.open((fn_root + ".training").c_str());
    if (!fh_training)
        REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                     fn_root + ".training" + " for output");
    fh_training << aux_model << std::endl;
    fh_training.close();

    // Save auto vectors
    aux_model = __auto_model;
    aux_model.import_particles(__auto_loaded_model);
    std::ofstream fh_auto;
    fh_auto.open((fn_root + ".auto").c_str());
    if (!fh_auto)
        REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                     fn_root + ".auto" + " for output");
    fh_auto << aux_model << std::endl;
    fh_auto.close();

    // Save error vectors
    std::ofstream fh_error;
    for (int i = 0; i < __Nerror_models; i++)
    {
        fh_error.open((fn_root + ".error" + integerToString(i, 1)).c_str());
        if (!fh_error)
            REPORT_ERROR(1, (std::string)"QtWidgetMicrograph::write: Cannot open file " +
                         fn_root + ".error" + integerToString(i, 1) + " for output");
        fh_error << __error_model.at(i) << std::endl;
        fh_error.close();
    }
}

/* Configure auto ---------------------------------------------------------- */
void QtWidgetMicrograph::configure_auto()
{

    QDialog   setPropertiesDialog(this, 0, TRUE);
    setPropertiesDialog.setCaption("Configure AutoSelect");
#ifdef QT3_SUPPORT
    Q3Grid     qgrid(2, &setPropertiesDialog);
#else
    QGrid     qgrid(2, &setPropertiesDialog);
#endif

    qgrid.setMinimumSize(250, 400);

    QLabel    lpiecexsize("Piece X size: ", &qgrid);
    QLineEdit piece_xsize(&qgrid);
    piece_xsize.setText(integerToString(__piece_xsize).c_str());

    QLabel    lpieceysize("Piece Y size: ", &qgrid);
    QLineEdit piece_ysize(&qgrid);
    piece_ysize.setText(integerToString(__piece_ysize).c_str());

    QLabel    lpiece_overlap("Piece overlap: ", &qgrid);
    QLineEdit piece_overlap(&qgrid);
    piece_overlap.setText(integerToString(__piece_overlap).c_str());

    QLabel    lcutoff("High pass cut-off: ", &qgrid);
    QLineEdit cutoff(&qgrid);
    cutoff.setText(floatToString(__highpass_cutoff).c_str());

    QLabel    loutput_scale("Output scale: ", &qgrid);
    QLineEdit output_scale(&qgrid);
    output_scale.setText(integerToString(__output_scale).c_str());

    QLabel    lmask_size("Mask size: ", &qgrid);
    QLineEdit mask_size(&qgrid);
    mask_size.setText(integerToString(__mask_size).c_str());

    QLabel    lgraybins("Gray bins: ", &qgrid);
    QLineEdit graybins(&qgrid);
    graybins.setText(integerToString(__gray_bins).c_str());

    QLabel    lradialbins("Radial bins: ", &qgrid);
    QLineEdit radialbins(&qgrid);
    radialbins.setText(integerToString(__radial_bins).c_str());

    QLabel    lparticle_radius("Particle radius: ", &qgrid);
    QLineEdit particle_radius(&qgrid);
    particle_radius.setText(integerToString(__particle_radius).c_str());

    QLabel    lmask_overlap("Mask overlap: ", &qgrid);
    QLineEdit mask_overlap(&qgrid);
    mask_overlap.setText(integerToString(__particle_overlap).c_str());

    QLabel    lnumin("Min. Harmonic: ", &qgrid);
    QLineEdit numin(&qgrid);
    numin.setText(integerToString(__numin).c_str());

    QLabel    lnumax("Max. Harmonic: ", &qgrid);
    QLineEdit numax(&qgrid);
    numax.setText(integerToString(__numax).c_str());

    QLabel    lmin_dist("Min. Distance: ", &qgrid);
    QLineEdit min_dist(&qgrid);
    min_dist.setText(integerToString(__min_distance_between_particles).c_str());

    QLabel    lkeep("Keep: ", &qgrid);
    QLineEdit keep(&qgrid);
    keep.setText(floatToString(__keep).c_str());

    QLabel    lNerror_models("#Error models: ", &qgrid);
    QLineEdit Nerror_models(&qgrid);
    Nerror_models.setText(integerToString(__Nerror_models).c_str());

    QPushButton okButton("Ok", &qgrid);
    QPushButton cancelButton("Cancel", &qgrid);

    connect(&okButton, SIGNAL(clicked(void)),
            &setPropertiesDialog, SLOT(accept(void)));
    connect(&cancelButton, SIGNAL(clicked(void)),
            &setPropertiesDialog, SLOT(reject(void)));

    if (setPropertiesDialog.exec())
    {
        __piece_xsize = piece_xsize.text().toInt();
        __piece_ysize = piece_ysize.text().toInt();
        __piece_overlap = piece_overlap.text().toInt();
        __highpass_cutoff = cutoff.text().toFloat();
        __output_scale = output_scale.text().toInt();
        __reduction = (int)pow(2.0, __output_scale);
        __mask_size = mask_size.text().toInt();
        __gray_bins = graybins.text().toInt();
        __radial_bins = radialbins.text().toInt();
        __particle_radius = particle_radius.text().toInt();
        __particle_overlap = mask_overlap.text().toInt();
        __numin = numin.text().toInt();
        __numax = numax.text().toInt();
        __min_distance_between_particles = min_dist.text().toInt();
        __keep = keep.text().toFloat();
        __Nerror_models = Nerror_models.text().toInt();
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
    __mImage->repaint(FALSE);
}

void QtWidgetMicrograph::changeCircleRadius(float _circle_radius)
{
    __ellipse_radius = _circle_radius;
    __mImage->__ellipse_radius = __ellipse_radius;
    __mImage->repaint(FALSE);
}

void QtWidgetMicrograph::repaint(int t)
{
    __mImage->repaint(FALSE);
    __mImageOverview->repaint(FALSE);
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
    repaint();
}

void QtWidgetMicrograph::slotChangeFamilyOther(int _coord, int _f)
{
    __m->coord(_coord).label = _f;
    repaint();
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
