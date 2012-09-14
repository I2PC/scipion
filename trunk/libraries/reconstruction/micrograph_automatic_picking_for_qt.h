/***************************************************************************
*
* Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
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

#ifndef MICROGRAPH_AUTOMATIC_PICKING_QT_H
#define MICROGRAPH_AUTOMATIC_PICKING_QT_H

#include <data/micrograph.h>
#include <data/mask.h>
#include <data/xmipp_program.h>
#include <classification/naive_bayes.h>

/// @defgroup AutomaticPickingQt Automatic Particle Picking for Qt
/// @ingroup ReconsLibrary
//@{
/* Particle ---------------------------------------------------------------- */
class ParticleQt
{
public:
    int x, y;                  // position in micrograph
    int idx;                   // Index of this particle within the micrograph
                               // list of coordinates
    char status;               // rejected=0, selected=1 or moved=2
    Matrix1D<double> vec;      // vector of that particle
    double cost;               // Associated cost

    // Print
    friend std::ostream & operator << (std::ostream &_out, const ParticleQt &_p);

    // Read
    void read(std::istream &_in, int _vec_size);
};

/* Classification model ---------------------------------------------------- */
class Classification_modelQt
{
public:
    // Maximum number of elements per training class
    int __maxTrainingVectors;

    // Example vectors
    std::vector< std::vector< ParticleQt > >        __training_particles;
    int                                           __classNo;
    int                                           __micrographs_number;
    std::vector<int>                              __micrographs_scanned;
    std::vector<int>                              __particles_picked;
    std::vector<int>                              __falsePositives;

public:
    EnsembleNaiveBayes *                          __bayesEnsembleNet;

public:
    /// Constructor
    Classification_modelQt(int _classNo=3, int _maxTrainingVectors=20000);

    // Clear
    void clear();

    /// Is empty
    inline bool isEmpty()
    {
        return __micrographs_number==0;
    }

    /// Initialize
    void init(int _classNo=3);

    /// Different additions
    inline void addMicrographScanned(int micrographScanned)
    {
        __micrographs_scanned.push_back(micrographScanned);
    }

    /// Add training particle
    bool addParticleTraining(const ParticleQt &p, int classIdx);

    /// Add number of particles picked
    inline void addParticlePicked(int particlePicked)
    {
        __particles_picked.push_back(particlePicked);
    }

    /// Add number of false positives
    inline void addFalsePositives(int falsePositives)
    {
        __falsePositives.push_back(falsePositives);
    }

    /// Add micrograph
    inline void addMicrographItem()
    {
        __micrographs_number++;
    }

    /// Import classification model
    void import_model(const Classification_modelQt &_model);

    /// Is a particle?
    inline int isParticle(const Matrix1D<double> &new_features, double &cost,
    					  Matrix1D<double> &aux1, Matrix1D<double> &aux2)
    {
        return __bayesEnsembleNet->doInferenceForClass(0,new_features,cost,aux1,aux2);
    }

    /// Init the naive bayesian network
    void initNaiveBayesEnsemble(const std::vector < MultidimArray<double> >
                                &features, const Matrix1D<double> &probs,
                                int discreteLevels, double penalization,
                                int numberOfClassifiers,
                                double samplingFeatures, double samplingIndividuals,
                                const std::string &newJudgeCombination);

    /// Print
    friend std::ostream & operator << (std::ostream &_out,
                                       const Classification_modelQt &_m);

    /// Read
    friend std::istream & operator >> (std::istream &_in,
                                       Classification_modelQt &_m);
};

/* Automatic particle picking ---------------------------------------------- */
/** Class to perform the automatic particle picking */
class AutoParticlePickingQt
{
public:
    Micrograph                *__m;
    FileName                   __modelRootName;
    FileName                   __outputRoot;
    int                        __numThreads;
    bool                       __learn_particles_done;
    bool                       __autoselection_done;
    Mask                       __mask;
    Classification_modelQt       __training_model;
    Classification_modelQt       __training_loaded_model;
    Classification_modelQt       __selection_model;
    Classification_modelQt       __selection_model2;
    Classification_modelQt       __selection_model3;
    int                        __auto_label;
    int                        __gray_bins;
    int                        __radial_bins;
    double                     __highpass_cutoff;
    double                     __penalization;
    float                      __minCost;
    int                        __piece_xsize;
    int                        __particle_radius;
    int                        __mask_size;
    int                        __min_distance_between_particles;
    int                        __output_scale;
    int                        __reduction; // Of the piece with respect to the micrograph
    int                        __piece_overlap;
    int                        __scan_overlap;
    int                        __learn_overlap;
    int                        __classNo;
    std::vector<ParticleQt>      __auto_candidates;
    std::vector<ParticleQt>      __rejected_particles;
    bool                       __is_model_loaded;
    std::vector < MultidimArray<int> * >    __mask_classification;
    std::vector < MultidimArray<int> * >    __radial_val;
    std::vector < MultidimArray<double> * > __sector;
    std::vector < MultidimArray<double> * > __ring;
    std::vector < MultidimArray<int> * >    __Nsector;
public:
    /// Empty constructor
    AutoParticlePickingQt(Micrograph *_m);

    /// Set the number of threads
    inline void setNumThreads(int _numThreads)
    {
        __numThreads=_numThreads;
    }

    /// Set output directory
    inline void setOutputRoot(const FileName &outputRoot)
    {
        __outputRoot=outputRoot;
    }

    /// Learn particles
    void learnParticles(int _ellipse_radius);

    /// Create mask for learning particles
    void createMask();

    /// Classify mask
    void classifyMask();

    /// Build vectors
    void buildPositiveVectors(std::vector<int> &_idx, Classification_modelQt &_model);

    /** Build vector from non particles */
    void buildNegativeVectors(Classification_modelQt &__model,
                              bool checkForPalsePostives);

    /** Build classification vector.
        x,y are in the coordinate system of the piece (that might be
        a reduced version of a piece in the micrograph)
        (0,0) is the top-left corner
        Returns true if the vector is successfully built */
    bool build_vector(const MultidimArray<int> &piece,
                      const MultidimArray<double> &original_piece,
                      int _x, int _y, Matrix1D<double> &_result);

    /** Get a piece of the micrograph centered at position x,y (if possible)
        the position of (x,y) in the piece is returned in (posx, posy) */
    void get_centered_piece(MultidimArray<double> &piece,
                            int _x, int _y, int &_posx, int &_posy);

    /** Get a piece whose top-left corner is at the desired position (if possible)
        Returns true if the piece could be taken, and false if the whole
        micrograph has been scanned. To scan the full micrograph,
        Top and left should be initialized to 0,0 and this function should
        succesive times. It returns the next top and left coordinates to
        get the next piece along with the skips. The skips indicate which
        part of the piece has been already scanned. This happens towards the
        right and bottom boundaries of the micrograph, where the pieces
        need to be shifted in order to fit with the required size.
        The overlap parameter defines what piece overlap we want. */
    bool get_corner_piece(MultidimArray<double> &piece,
                          int _top, int _left, int _skip_y,
                          int &_next_skip_x, int &_next_skip_y, int &_next_top,
                          int &_next_left, int overlap, bool copyPiece) const;

    /** Filter, denoise, reject outliers and equalize histogram.
        Returns true, if successful. False if unsuccessful (skip this piece)
        Usually, it is unsuccessful if the denoising fails to work because
        some "weird" features of the piece */
    bool prepare_piece(MultidimArray<double> &piece,
    				   MultidimArray<int> &ipiece,
                       MultidimArray<double> &original_piece);

    /** To get the neighbours of the particle at position (x,y) in the micrograph
        (with actual coordinates in the piece posx,posy)
        and their positions in the piece image.
        */
    void find_neighbour(const MultidimArray<double> &piece,
                        std::vector<int> &_idx, int _index,
                        int _x, int _y,
                        int _posx, int _posy, MultidimArray<char> &_visited,
                        std::vector< Matrix1D<int> > &_nbr);

    /** Automatically Select Particles.
       Returns the number of particles selected.
       */
    int automaticallySelectParticles();

    /// check if there are any particles in the actual scanning position
    bool anyParticle(int posx, int posy, int rect_size);

    /** Count the number of scanning positions */
    int count_scanning_pos();

    /** Given a current scanning position, this function returns
        the next scanning position whithin the current piece.
        The skips are given by get_corner_piece.
        It returns true, if the next scanning position can be computed.
        Otherwise, if the piece has been completely scanned, it returns false
        Initialize _x,_y to 0,0 to scan the full piece (even if there are skips).
        The overlap parameter defines what particle overlap we want.
        */
    bool get_next_scanning_pos(
        const MultidimArray<double> &piece,
        int &_x, int &_y, int _skip_x, int _skip_y, int overlap) const;

    /** Run over the list sorted by distances.
     * If two particles are within
        a given distance then either reject both or the one with largest distance
        depending upon _reject_both. This function returns the number of particles
       that are still candidates.
       */
    int reject_within_distance(
        std::vector<ParticleQt> &_Input, double _min_dist, bool _reject_both);

    /// load models with a name
    void loadModels(const FileName &fn);

    /// Save models
    void saveModels(const FileName &fn) const;

    /// Save automatically selected particles
    void saveAutoParticles() const;

    /// Save the feature vectors of the automatically selected particles
    void saveAutoFeatureVectors(const FileName &fn) const;

    /// Load the feature vectors of the automatically selected particles
    void loadAutoFeatureVectors(const FileName &fn);

    /// Get Features of the classification model
    void getFeatures(const Classification_modelQt &_model,
                         std::vector < MultidimArray<double> > &_features);

    /// Get classes probabilities
    void getClassesProbabilities(const Classification_modelQt &_model,
                                     Matrix1D<double> &probabilities);

    /// Get false positives automatically selected
    void getAutoFalsePositives(Classification_modelQt &_training_model);

    /// Get true positives automatically selected
    void getAutoTruePositives(Classification_modelQt &_training_model);

    /// Restrict selection
    void restrictSelection(float _cost);

    /** Move particle.
        The input index is the index of the moved particle in the micrograph list */
    void move_particle(int _idx);

    /** Delete particle.
        The input index is the index of the moved particle in the micrograph list */
    void delete_particle(int _idx);
};

/// AutomaticallySelectThreadParams
struct AutomaticallySelectThreadParamsQt
{
    AutoParticlePickingQt *autoPicking;
    int idThread;
};

/// Automatically Select Particles thread
void * automaticallySelectParticlesThreadQt(void *);

//@}
#endif
