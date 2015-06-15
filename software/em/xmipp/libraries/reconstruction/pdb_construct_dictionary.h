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
#ifndef _PROG_CONSTRUCT_DICTIONARY_HH
#define _PROG_CONSTRUCT_DICTIONARY_HH

#include <data/xmipp_program.h>

/**@defgroup PDBConstructDictionary Construct a low and high resolution dictionary
   @ingroup ReconsLibrary */
//@{
/** Generic class to handle PDB Low and High resolution dictionary*/
class ProgPDBDictionary: public XmippProgram
{
public:
	/** Dictionary rootname */
    FileName fnRoot;

    /** Patch is of size size x size x size */
    int patchSize;

    /** A patch is candidate if its standard deviation is at least this factor of the total standard deviation */
    double stdThreshold;

    double angleThreshold;

    // Regularization
    double lambda;

    // Number of iterations
    int iterations;

    // 3D or 2D mode
    int mode;

public:
    /** Low resolution and high resolution dictionary */
    std::vector< MultidimArray<double> > dictionaryLow, dictionaryHigh;

    /** Signature dictionary */
    std::vector< Matrix1D<double> > dictionarySignature;

    /** Rotation group */
    std::vector< Matrix2D<double> > rotationGroup;

    /** Auxiliary patch */
    MultidimArray<double> auxPatch;

    /** Signature */
    Matrix1D<double> auxSignature;

    /* Some variables to approximate patches */
    Matrix2D<double> Ui, UitUi;
	Matrix1D<double> wi, v1, v2, y, yp, Uity;
	int patchSize_2;

public:
    virtual void defineParams();
    virtual void readParams();
    virtual void show();
    virtual void run()=0;

    /** Extract Patch from volume */
    void extractPatch(const MultidimArray<double> &V, MultidimArray<double> &patch, int k, int i, int j);

    /** Insert Patch into volume */
    void insertPatch(MultidimArray<double> &Vhigh, MultidimArray<double> &weightHigh, const MultidimArray<double> &patchHigh,
    		int k, int i, int j, double R2);

    /** Construct rotation group 2D */
    void constructRotationGroup2D();

    /** Construct rotation group 3D */
    void constructRotationGroup3D();

    /** Construct rotation group */
    void constructRotationGroup();

    /** Orient canonically a patch */
    size_t canonicalOrientation2D(const MultidimArray<double> &patch, MultidimArray<double> &canonicalPatch,
    		Matrix1D<double> &patchSignature);

    /** Orient canonically a patch */
    size_t canonicalOrientation3D(const MultidimArray<double> &patch, MultidimArray<double> &canonicalPatch,
    		Matrix1D<double> &patchSignature);

    /** True if the patch is not already in the low resolution dictionary */
    bool notInDictionary(const MultidimArray<double> &candidatePatch, MultidimArray<double> &canonicalPatch,
    		Matrix1D<double> &canonicalSignature, size_t &canonicalIdx);

    /** Select dictionary patches for a low resolution patch */
    void selectDictionaryPatches(const MultidimArray<double> &lowResolutionPatch, Matrix1D<double> &lowResolutionPatchSignature,
    		std::vector<size_t> &selectedPatchesIdx, std::vector<double> &weight);

    /** Approximate patch.
     * It returns the R2 of the approximation. */
    double approximatePatch(const MultidimArray<double> &lowResolutionPatch,
    		std::vector< size_t > &selectedPatchesIdx, std::vector<double> &weight, Matrix1D<double> &alpha);

    /** Reconstruct patch */
    void reconstructPatch(size_t idxTransf, std::vector< size_t > &selectedPatchesIdx, Matrix1D<double> &alpha,
    		MultidimArray<double> &highResolutionPatch);

    /** Load dictionaries */
    void loadDictionaries();

    /** Save dictionaries */
    void saveDictionaries() const;
};

/** Construct Low and High resolution dictionary*/
class ProgConstructPDBDictionary: public ProgPDBDictionary
{
	/** Metadata with the low resolution volumes */
    FileName fnLow;

    /** Metadata with the high resolution volumes */
    FileName fnHigh;

    double R2threshold;
public:
    void defineParams();
    void readParams();
    void show();
    void run();
};
//@}
#endif
