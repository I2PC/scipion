/***************************************************************************
 *
 * Authors:      Javier Angel Velazquez Muriel    javi@cnb.uam.es
 *               Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#ifndef _ADJUST_CTF_HH
   #define _ADJUST_CTF_HH

#ifdef _HAVE_VTK
#include "Prog_FourierFilter.hh"

/**@name Adjust parametric CTF */
//@{
/** Adjust CTF parameters. */
class Adjust_CTF_Parameters
{
public:
   /// CTF filename
   FileName             fn_ctf;
   /// Output rootname
   FileName             fn_outroot;
   /// Output CTF parameter file
   FileName             fn_out_CTF_parameters;
   /// CTF amplitude to model
   FourierMask          ctftomodel;
   /// CTF model
   XmippCTF             ctfmodel;
   /// Show convergence values
   bool                 show_optimization;
   /// Do not optimize, only show the adjustment
   bool                 do_not_optimize;
   /// Do not fine search
   bool                 do_not_fine_search;
   /// Allow astigmatic noise
   bool                 astigmatic_noise;
   
   /// Penalty for trespassing CTF in gaussian adjust 
   double               penalty;
   /// Weight for the central (first two CTF lobes) region
   double               central_weight;
   /// Minimum frequency to adjust 
   double               min_freq;
   /// Maximum frequency to adjust 
   double               max_freq;
   /** Gamma correction.
       This factor is applied during the determination of the
       defoci parameters when there is no hint about their
       values. The higher it is the more importance it is given
       to low frequencies. A default value of 3 is recommended.*/
   double               gamma;
   /// Accuracy. Stop criterion for iterations 
   double               accuracy;
   /// Minimum threshold for adjusting
   double               value_th;
   /// Sampling rate
   double               Tm;
   /// Evalaute fitness functions every n pixels, where n comes from evaluation_reduction
   int                  evaluation_reduction;
   /// Set of parameters for the complete adjustment of the cTF
   matrix1D<double>     adjust;
   /// Set of selected parameters to adjust
   matrix1D<double>     steps;

public:
   /// Read parameters from file
   void read(const FileName &fn_param) _THROW;
   
   /// Write to a file
   void write(const FileName &fn, bool rewrite=TRUE) _THROW;

   /// Show parameters
   void show();
   
   /// Usage
   void Usage();

   /// Produce side information
   void produce_side_info();
};

/** Core of the Adjust CTF routine.
    This is the routine which does everything. */
void ROUT_Adjust_CTF(Adjust_CTF_Parameters &prm);
//@}
#endif
#endif
