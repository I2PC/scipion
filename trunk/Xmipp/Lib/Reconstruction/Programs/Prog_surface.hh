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
#ifndef _PROG_SURFACE_HH
#  define _PROG_SURFACE_HH

#include <XmippData/xmippFuncs.hh>
#include <Reconstruction/phantom.hh>

/**@name Surface program */
//@{
/* Surface Program Parameters ---------------------------------------------- */
/** Parameter class for the project program */
class Prog_Surface_Parameters {
public:
   /// Filename with the \Ref{Phantom}.
   FileName fn_phantom;
   /// Probe radius
   double probe_radius;
   /// Enable ztop mode
   bool enable_ztop;
   /// Maximum "ray" height, if ztop==zdim then it is not used
   float ztop;
   /// Enable zbottom mode
   bool enable_zbottom;
   /// Maximum "ray" height, if ztop==zdim then it is not used
   float zbottom;
   /// Z dimension of output volume
   int zdim;
   /// Filename of top surface
   FileName fn_top;
   /// Filename of bottom surface
   FileName fn_bottom;
   /// Filename of output mask
   FileName fn_mask;   
   /**@name Side parameters */
   //@{
   /// Phantom
   Phantom phantom;
   /// Top surface
   ImageXmipp top_surface;
   /// Bottom surface
   ImageXmipp bottom_surface;
   //@}
public:
   /** Read from a command line.
       An exception might be thrown by any of the internal conversions,
       this would mean that there is an error in the command line and you
       might show a usage message.
      
       An exception is thrown if there is no surface specification or if
       there is no phantom and the Z dimension is not given*/
   void read(int argc, char **argv);

   /** Usage message.
       This function shows the way of introducing this parameters. */
   void usage() const;

   /** Produce Side Information.
       Read phantom file and assign ztop and zbottom if they
       are not assigned by user, ie, set it to zdim. */
   void produce_Side_Info();
};

/** Create mask from two surfaces.
    An exception is thrown if the two surfaces are not of the same shape.
    The output volume is resized to the image shape plus the zdim information,
    and the Xmipp origin is set on the Z direction.
    
    This function needs that both surfaces are non empty, if you don't
    want to use one of them set it to NULL, it makes the same effect as 
    covering the whole volume. */
    void create_surface_mask(const Image *top, const Image *bottom,
       int zdim, Volume *V);

/** Run surface.
    This function is who really creates the surface for the phantom. It is
    very simple and lies on \Ref{Phantom::surface} or \Ref{surface_mask}*/
   void ROUT_surface(Prog_Surface_Parameters &prm);
//@}
#endif
