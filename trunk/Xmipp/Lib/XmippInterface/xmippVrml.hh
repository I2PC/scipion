/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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
//Draw spheres in vrml
//lines will be follow soon
#ifndef _XMIPPVRML_HH
#   define _XMIPPVRML_HH

#include <stdlib.h>
#include <fstream>  
#include <iostream>
#include <stdio.h>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippMatrices1D.hh>
/**@name VRML */
//@{
/**  VRML
      This is a class to draw a few vrml objects, designed with the idea
      of visualizing grid points.
    */
class VrmlFile {
  private:
    ofstream fh_out;
    int need_footer; //1=YES, 0=NO
  public:
    /**Open file and write header. Default filename is myvrml.wrl"
      */      
   VrmlFile(FileName vrmlfilename="myvrml" );
   
   /** Draw coordinate system, x axis is red, y axis is green, z axis is blue.
       Parameters are the scale of the coordinate system and the relative size
        between the cone and the cylinder. Defaults for both parameters are 1.
       */ 
   void Axis( double scale=1.,double Cone_scale=1.);
   /** Draw a first sphere.  By default selects color red and places an sphere
       at the center with radious 0.05.
       */ 
   void Sphere( const matrix1D<double> RGB=vector_R3(1.,0.,0.),
                const matrix1D<double> XYZ=vector_R3(0.,0.,0.),
                const double radius=0.05);

   /** Draw a translated copy of the previous sphere. Use this only after
       Sphere has been called once
     */ 
   void Add_sphere(const matrix1D<double> XYZ=vector_R3(0.,0.,0.));

   /** Draw a transparent sphere. Some wrl browser does not handle 
   transparency properly  */
   void Trans_Sphere(const matrix1D<double> XYZ=vector_R3(1.,0.,0.),
                     const matrix1D<double> RGB=vector_R3(0.,0.,0.),
                     const double radius=0.05);
   /** write a footer and Close file
       */
   ~VrmlFile()
   {
   fh_out<< "    ]\n}";         
   fh_out.close(); // all done
   }   
   
   };
//@}

#endif
