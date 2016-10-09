/***************************************************************************
 *
 * Authors:    Roberto Marabini roberto@cnb.csic.es           
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
/* This program is strongly based on the umblur code
 * (http://grigoriefflab.janelia.org/unblur)
which has been realeased under the following copyright:

The Janelia Farm Research Campus Software Copyright 1.1

Copyright (c) 2014, Howard Hughes Medical Institute, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:


Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
Neither the name of the Howard Hughes Medical Institute nor the names of
its contributors may be used to endorse or promote products derived from
this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, ANY
IMPLIED WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
REASONABLE ROYALTIES; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
*/
#ifndef _PROG_MOVIE_FILTER_DOSE
#define _PROG_MOVIE_FILTER_DOSE

#include <data/xmipp_program.h>
#include <cmath>
#include <limits>
#include <complex>      // std::complex, std::abs

/**@defgroup MovieAlignmentCorrelation Movie alignment by correlation
   @ingroup ReconsLibrary */
//@{

/** Movie alignment correlation Parameters. */
class ProgMovieFilterDose: public XmippProgram
{
public:

	/**  Set up the critical curve function (summovie */
    double critical_dose_a;
	double critical_dose_b;
	double critical_dose_c;

	double critical_dose_at_dc;

	double voltage_scaling_factor;
	double acceleration_voltage;
	double dose_per_frame;
	double pre_exposure_amount;

    /** input movie */
    FileName user_supplied_input_filename;

    /** output movie */
    FileName user_supplied_output_filename;

	/** frames of interest in movie*/
	int user_supplied_first_frame;
	int user_supplied_last_frame;

    /** Sampling rate */
    double pixel_size;

    /** Max freq. */
    double maxFreq;


	/// Restore noise power after filtering?', 'Renormalise the summed image after filtering
    bool restore_power;

    // Fourier transforms of the input images
	std::vector< MultidimArray<std::complex<double> > * > frameFourierVec;

public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /// Run
    void run();

    /// constructor
    ProgMovieFilterDose(double accelerationVoltage);
    ProgMovieFilterDose(void);
    void init(void);
    void initVoltage(double accelerationVoltage);

    ///Compute the dose filter, which is the signal attenuation
    // factor due to radiation damage (as used by summovie)
    double doseFilter(double dose_at_end_of_frame, double critical_dose);

    ///  Given a spatial frequency, return the critical dose in electrons per square Angstroms
    double criticalDose(double spatial_frequency);


   /// Given the critical dose, return an estimate of the optimal dose (at which the SNR is maximised)
   double optimalDoseGivenCriticalDose(double critical_dose);

   /// Apply a dose filter to the image Fourier transform
   void applyDoseFilterToImage(
		int Ydim, int Xdim,
   		const MultidimArray< std::complex<double> > &FFT1,
		const double dose_start, const double dose_finish
   		);

};
//@}
#endif
