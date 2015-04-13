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
#ifndef _PROG_MICROSCOPE_HH
#define _PROG_MICROSCOPE_HH

#include "fourier_filter.h"

#include <data/metadata_extension.h>

/**@defgroup MicroscopeProgram phantom_simulate_microscope (Microscope simulation)
   @ingroup ReconsLibrary */
//@{
/* Microscope Program Parameters ------------------------------------------- */
/** Parameter class for the project program */
class ProgSimulateMicroscope: public XmippMetadataProgram
{
public:
    /// Filename with the CTF
    FileName fn_ctf;
    /// Total noise power
    double   sigma;
    /// Low pass frequency before CTF
    double   low_pass_before_CTF;
    /// Filename with the root squared spectrum for noise after CTF
    bool     after_ctf_noise;
    /// Defocus change (%)
    double   defocus_change;
    /// Target SNR
    double targetSNR;
    /// Estimate SNR
    bool estimateSNR;
    /// CTFpresent
    bool CTFpresent;
    /// CTF
    FourierFilter ctf;
    /// Low pass filter, if it is 0 no lowpass filter is applied
    FourierFilter lowpass;
    /// After CTF noise root squared spectrum
    FourierFilter after_ctf;
    /// Noise power before CTF
    double   sigma_before_CTF;
    /// Noise power after CTF
    double   sigma_after_CTF;
    /// Input image Xdim
    size_t   Xdim;
    /// Input image Ydim
    size_t   Ydim;
    /// Particular reference to mdIn to manipulated
    MetaData * pmdIn;
    /** Downsampling factor */
    double downsampling;
    /* save U defocus in case we randomize it */
    double defocusU;
    /* save V defocus in case we randomize it */
    double defocusV;

public:
    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void readParams();

    /** Usage message.
        This function shows the way of introducing these parameters. */
    void defineParams();

    /** Show parameters. */
    void show();

    /** Estimate sigma for a given SNR */
    void estimateSigma();

    /** Prepare fourier filter */
    void setupFourierFilter(FourierFilter &filter, bool isBackground, double &power);

    /** Update ctfs */
    void updateCtfs();

    /** Produce side information. */
    void preProcess();


    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Apply to a single image. The image is modified.
        If the CTF is randomly selected then a new CTF is generated
        for each image */
    void apply(MultidimArray<double> &I);

};
//@}
#endif
