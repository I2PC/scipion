/***************************************************************************
 *
 * Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#ifndef METADATACONTAINER_H
#define METADATACONTAINER_H

#include <map>
#include "strings.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "funcs.h"

// This enum defines what MetaDataLabels this class can manage, if
// you need a new one add it here and modify affected methods:
//
//  - static MetaDataLabel codifyLabel( std::string strLabel );
// - static std::string decodeLabel( MetaDataLabel inputLabel );
// - void writeValuesToFile( std::ofstream &outfile, MetaDataLabel inputLabel );
// - void addValue( std::string name, std::string value );
//
// Keep this special structure (using MDL_FIRSTLABEL and MDL_LAST_LABEL) so the
// programmer can iterate through it like this:
//
//  for( MetaDataLabel mdl = MDL_FIRST_LABEL ; mdl < MDL_LAST_LABEL ; MetaDataLabel( mdl+1 ) )
//
enum MetaDataLabel
{
    MDL_UNDEFINED = -1,
    MDL_FIRST_LABEL,
    MDL_ANGLEPSI = MDL_FIRST_LABEL, // Psi angle of an image (double)
    MDL_ANGLEROT, // Rotation angle of an image (double)
    MDL_ANGLETILT, // Tilting angle of an image (double)
    MDL_ANGLEPSI2, // Psi angle of an image (double)
    MDL_ANGLEROT2, // Rotation angle of an image (double)
    MDL_ANGLETILT2, // Tilting angle of an image (double)
    MDL_COMMENT, // A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
    MDL_CTFINPUTPARAMS, // Parameters file for the CTF Model (std::string)
    MDL_CTFMODEL, // Name for the CTF Model (std::string)
    MDL_ENABLED, // Is this image enabled? (int [-1 or 1])
    MDL_FLIP, // Flip the image? (bool)
    MDL_IMAGE, // Name of an image (std::string)
    MDL_MAXCC, // Cross-correlation for the image (double)
    MDL_COST, // Cost for the image (double)
    MDL_MICROGRAPH, // Name of a micrograph (std::string)
    MDL_NMA, // Normal mode displacements
    MDL_ORIGINX, // Origin for the image in the X axis (double)
    MDL_ORIGINY, // Origin for the image in the Y axis (double)
    MDL_ORIGINZ, // Origin for the image in the Z axis (double)
    MDL_PERIODOGRAM, // A periodogram's file name (std::string)
    MDL_PMAX, // Maximum value of normalized probability function (now called "Pmax/sumP") (double)
    MDL_REF, // Class to which the image belongs (int)
    MDL_SCALE, // scaling factor for an image or volume (double)
    MDL_BGMEAN, // Mean background value for an image
    MDL_INTSCALE, // Intensity scale for an image
    MDL_MODELFRAC, // Model fraction (alpha_k) for a Maximum Likelihood model
    MDL_MIRRORFRAC, // Mirror fraction for a Maximum Likelihood model
    MDL_LL, // contribution of an image to log-likelihood value
    MDL_WROBUST, // Weight of t-student distribution in robust Maximum likelihood
    MDL_SIGNALCHANGE, // Signal change for an image
    MDL_SERIE, // A collection of micrographs, e.g. a tilt serie (std::string)
    MDL_SHIFTX, // Shift for the image in the X axis (double)
    MDL_SHIFTY, // Shift for the image in the Y axis (double)
    MDL_SHIFTZ, // Shift for the image in the Z axis (double)
    MDL_X, // X component (double)
    MDL_Y, // Y component (double)
    MDL_Z, // Z component (double)
    MDL_WEIGHT, // Weight assigned to the image (double)
    MDL_OBJID, // object id (int)
    //add row only label at the end of the enum
    MDL_SAMPLINGRATE, // sampling rate (double)
    MDL_VOLTAGE, // microscope voltage (double)
    MDL_DEFOCUSU, // microscope defocus U direction (double)
    MDL_DEFOCUSV, // microscope defocus V direction (double)
    MDL_IMGMD, // Name of Metadata file for all images
    MDL_REFMD, // Name of Metadata file for all references
    MDL_ITER, // Current iteration number
    MDL_BLOCK, // Current block number (for incremental EM)
    MDL_SIGMANOISE, // Standard deviation of the noise in ML model
    MDL_SIGMAOFFSET, // Standard deviation of the offsets in ML model
    MDL_SUMWEIGHT, // Sum of all weights in ML model
    MDL_RANDOMSEED, // Seed for random number generator
    /*
     MDL_azimuthal_angle=      72.3493
     MDL_spherical_aberration= 5.6
     MDL_chromatic_aberration= 1.99957
     MDL_energy_loss=          0.0240301
     MDL_lens_stability=       0
     MDL_convergence_cone=     0.000329129
     MDL_longitudinal_displace=8.66588e-05
     MDL_transversal_displace= 4.14845
     MDL_Q0=                   -0.1
     MDL_K=                    2.13333
     MDL_gaussian_K=           2.58626
     MDL_sigmaU=               100000
     MDL_sigmaV=               85359.1
     MDL_cU=                   0.00332111
     MDL_cV=                   0.0132845
     MDL_gaussian_angle=       2.32559e-11
     MDL_sqrt_K=               70.1711
     MDL_sqU=                  20.3219
     MDL_sqV=                  19.1215
     MDL_sqrt_angle=           67.1273
     MDL_base_line=            0.368481
     MDL_gaussian_K2=          0.24108
     MDL_sigmaU2=              6489.64
     MDL_sigmaV2=              5825.19
     MDL_cU2=                  0.06172
     MDL_cV2=                  0.0586838
     MDL_gaussian_angle2=      90
     */
    MDL_LAST_LABEL                       // **** NOTE ****: Do keep this label always at the end
    // it is here for looping purposes
};

inline bool isString(MetaDataLabel lCode)
{
    if (lCode == MDL_COMMENT  || lCode == MDL_IMAGE          || lCode == MDL_MICROGRAPH  ||
        lCode == MDL_CTFMODEL || lCode == MDL_CTFINPUTPARAMS || lCode == MDL_PERIODOGRAM ||
        lCode == MDL_SERIE    || lCode == MDL_IMGMD          || lCode == MDL_REFMD)
        return true;
    else
        return false;
}

inline bool isDouble(MetaDataLabel lCode)
{
    if (lCode == MDL_ANGLEROT     || lCode == MDL_ANGLETILT    || lCode == MDL_ANGLEPSI    ||
        lCode == MDL_ANGLEROT2    || lCode == MDL_ANGLETILT2   || lCode == MDL_ANGLEPSI2   ||
        lCode == MDL_SHIFTX       || lCode == MDL_SHIFTY       || lCode == MDL_SHIFTZ      ||
        lCode == MDL_ORIGINX      || lCode == MDL_ORIGINY      || lCode == MDL_ORIGINZ     ||
        lCode == MDL_X            || lCode == MDL_Y            || lCode == MDL_Z           ||
        lCode == MDL_WEIGHT       || lCode == MDL_MAXCC        || lCode == MDL_PMAX        ||
        lCode == MDL_SCALE        || lCode == MDL_COST         || lCode == MDL_BGMEAN      ||
        lCode == MDL_INTSCALE     || lCode == MDL_SAMPLINGRATE || lCode == MDL_MODELFRAC   ||
        lCode == MDL_MIRRORFRAC   || lCode == MDL_VOLTAGE      || lCode == MDL_DEFOCUSU    ||
        lCode == MDL_DEFOCUSV     || lCode == MDL_LL           || lCode == MDL_WROBUST     ||
        lCode == MDL_SIGNALCHANGE || lCode == MDL_SIGMANOISE   || lCode == MDL_SIGMAOFFSET ||
        lCode == MDL_SUMWEIGHT)
        return true;
    else
        return false;
}

inline bool isVector(MetaDataLabel lCode)
{
    if (lCode == MDL_NMA)
        return true;
    else
        return false;
}

inline bool isBool(MetaDataLabel lCode)
{
    if (lCode == MDL_FLIP)
        return true;
    else
        return false;
}

inline bool isInt(MetaDataLabel lCode)
{
    if (lCode == MDL_REF || lCode == MDL_ENABLED || lCode == MDL_OBJID || lCode
        == MDL_ITER || lCode == MDL_BLOCK ||lCode == MDL_RANDOMSEED)
        return true;
    else
        return false;
}

class MetaDataContainer
{
    /** Container for pairs "name" and value. Note that void * allows to use
     mixed types */
    std::map<MetaDataLabel, void *> values;

    void insertVoidPtr(MetaDataLabel name, void * value);
    void * getVoidPtr(MetaDataLabel name);

public:

    /**Assignment operator
     *
     */
    MetaDataContainer& operator =(MetaDataContainer &MDc);

    /** Constructor */
    MetaDataContainer();
    /** Copy constructor
     *
     */
    MetaDataContainer(MetaDataContainer &MDc);

    /** Destructor */
    ~MetaDataContainer();

    /** Create a new pair name-value of integer type */
    void addValue(const std::string &name, const std::string &value);

    template<class T>
    void addValue(MetaDataLabel name, const T &value)
    {
        void * newValue = (void *) (new T(value));
        insertVoidPtr(name, newValue);
    }

    template<class T>
    void getValue(MetaDataLabel name, T &value)
    {
        std::map<MetaDataLabel, void *>::iterator element;

        element = values.find(name);

        if (element == values.end())
        {
            REPORT_ERROR(1,(std::string) "Label(int) " + decodeLabel(name) + " not found\n" );
        }
        else
        {
            value = *((T *) element->second);
        }
    }
    bool valueExists(MetaDataLabel name);

    //string is not part of the template because - is not defined for string
    bool pairExists(MetaDataLabel name, const std::string &value);

    template<class T>
    bool pairExists(MetaDataLabel name, T value)
    {
        // Traverse all the structure looking for objects
        // that satisfy search criteria
        std::map<MetaDataLabel, void *>::iterator It;

        It = values.find(name);

        if (It != values.end())
        {
            if (ABS( *((T *)(It->second)) - value )
                < XMIPP_EQUAL_ACCURACY)
            {
                return true;
            }
        }

        return false;
    }



    void deleteValue(MetaDataLabel name);

    bool writeValueToStream(std::ostream &outstream, MetaDataLabel inputLabel);
    bool writeValueToFile(std::ofstream &outfile, MetaDataLabel inputLabel);
    bool writeValueToString(std::string &outString, MetaDataLabel inputLabel);

    static MetaDataLabel codifyLabel(std::string strLabel);
    static std::string decodeLabel(MetaDataLabel inputLabel);
    static bool isValidLabel(MetaDataLabel inputLabel);
    static bool isValidLabel(std::string inputLabel);
};

#endif
