/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.uam.es)
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

#ifndef _PROG_ASSIGN_CTF
   #define _PROG_ASSIGN_CTF

   #include "Prog_adjust_CTF.hh"
   #include "Prog_SpARMA.hh"
/**@name Assign CTF.
   This program assign different CTFs to the particles in a micrograph */
//@{

/** Assign CTF parameters. */
class Prog_assign_CTF_prm {
public:
   /// Parameters for adjust_CTF program
   Adjust_CTF_Parameters   adjust_CTF_prm;
   /// Parameters for ARMA
   ARMA_parameters         ARMA_prm;
   /// Reversed endian
   bool              	   reversed;
   /// Downsample factor
   int               	   downsampling;
   /// X dimension of micrograph pieces
   int                     N_horizontal;
   /// Y dimension of micrograph pieces
   int                     N_vertical;
   /// X dimension of particle projections
   int			   particle_horizontal; 
   /// Y dimension of particle projections
   int			   particle_vertical;
   /// Selfile with all projections
   FileName                selfile_fn;
   /// Position file
   FileName                picked_fn;
   /// ARMA files root
   FileName                ARMAfile;
   /// CTF root for particle projections
   FileName                CTFfile;
   /// Micrograph filename
   FileName                image_fn;
   /// Do not perform CTF estimation, only interpolate
   bool                    only_interpolate;
   /** the center of the windows in which the CTF is computed
       are the particles (stored at the .pos file) instead of 
       a regular grid. By default this is false.
   */
   bool                    compute_at_particle;
   /** Perform ARMA averaging.
       Like Periodogram averaging.*/
   bool                    ARMA_averaging;
   /** Perform Periodogram averaging */
   bool                    Periodogram_averaging;
public:
   /** Read parameters from file.
       If do_not_read_files is TRUE then all FileNames parameters
       are not read */
   void read(const FileName &fn_prm, bool do_not_read_files=FALSE) _THROW;

   /// Write parameters to file
   void write(const FileName &fn_prm) _THROW;

   /// Process the whole thing
   void process();
};
//@}
#endif
