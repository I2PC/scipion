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
/*****************************************************************************/
/* INTERACTION WITH SPIDER                                                   */
/*****************************************************************************/

#ifndef _RECONS_SPIDER_HH
#define _RECONS_SPIDER_HH

#include <data/metadata.h>
#include <data/funcs.h>

/**@defgroup SpiderRecons Interaction with Spider 
   @ingroup ReconsLibrary */
//@{
/** Call SIRT in Spider.
    The batch file for Spider is called fn_batch+"."+fn_ext. The output
    volume is fn_recons_root+".vol". There is two intermidiate files:
    \\ fn_root+"sel."+fn_ext: that is a selection file for Spider
    \\ fn_root+"ang."+fn_ext: which contains the angles.
    \\Notice that the angles passed to Spider are reversed in order (psi, tilt,
    rot) instead of the usual (rot, tilt, psi).
    All LOG and results file from Spider are not
    removed. The batch file generated is the following:
    @code
      bp rp
      <first file in SF>
      <fn_root>sel.<fn_ext>
      <radius>
      <fn_root>ang.<fn_ext>
      n
      <fn_root>
      (<lambda>,0)
      (<no_it>,1)
      (0.,2.)
      (0.95*1/(1+6*<lambda>))
      en
    @endcode
    The actual process status is written to std::cerr.
*/
void SIRT_Spider(MetaData &SF, double lambda, double no_it, int radius,
                 const FileName &fn_root, const FileName &fn_ext,
                 const FileName &fn_recons_root, const FileName &fn_batch);
//@}
#endif
