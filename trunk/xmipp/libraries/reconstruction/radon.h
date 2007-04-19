/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
/* This file contains functions related to the Radon Transform */

/* Prototypes -------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* Radon Transform                                                           */
/* ------------------------------------------------------------------------- */
/* The Radon transform is computed in two steps:                             */
/*    1) Projection 2D of the volume in a direction perpendicular to the one */
/*       we want to compute the RT                                           */
/*    2) Projection 1D of the 2D projection along the RT direction           */
/* INPUT                                                                     */
/*    vol:                 volume we want to compute the Radon Transform     */
/*    (x,y,z):             direction along which the Radon Transform will be */
/*                         calculated                                        */
/*    normalize:           if some bits are 1 then the RT will be normalized */
/*                            bit 0, in a statistical way (avg=0, stddev=1)  */
/*                               for the Projection 2D                       */
/*                            bit 1, in a range way (0,1) for the proj. 1D   */
/* OUTPUT                                                                    */
/*    RT:                  Radon Transform. It is supposed to be already     */
/*                         allocated                                         */

#ifndef _RADON_HH
#   define _RADON_HH

#include <data/volume.h>

/**@name Radon Transform */
/** Radon transform of a volume along a direction */
void Radon_Transform(Volume *vol, double rot, double psi,
   matrix1D<double> &RT);

/** Radon transform of a piece of volume along a direction */
void Local_Radon_Transform(Volume *vol, double rot, double tilt,
   int label, Volume *vol_label, matrix1D<double> &RT,
   matrix1D<double> &RT_n);

/** Radon transform of an image. */
void Radon_Transform(const matrix2D<double> &I, double rot_step,
   matrix2D<double> &RT);
#endif
