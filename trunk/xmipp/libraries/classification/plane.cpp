/***************************************************************************
 *
 * Authors:     Alberto Pascual
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
//-----------------------------------------------------------------------------
// xmippPlanes.cc
//-----------------------------------------------------------------------------

#include "plane.h"

/** getPlane: gets the "plane" (a map representing the effect of the
*   given variable.
*   @param _in: input SOM
*   @param _out: output plane
*   @param _plane: variable
*/

void xmippPlanes::getPlane(const In& _in, Out& _out, unsigned _plane) const
{
  unsigned i, j;


  // Check Maps dimensions

    if ((_out.height() != _in.height()) || (_out.width() != _in.width()))
      throw invalid_argument("xmippPlanes: Invalid Plane dimensions.");

  // Check if _plane is valid

    if (_plane > _in.itemAtPos(SomPos(0, 0)).size())
      throw invalid_argument("xmippPlanes: _plane parameter > Codevector dimension.");

  // Set calibrated tag.

    if (_in.calibrated())
      _out.calibrated(true);

  // find the minimum and maximum values for gray level scaling */

    xmippFeature Max = -MAXFLOAT;
    xmippFeature Min = MAXFLOAT;

    for (i=0;i<_in.width();i++)
      for (j=0;j<_in.height();j++)
        {
          if (_in.itemAtPos(SomPos(i, j))[_plane] > Max)
            Max = _in.itemAtPos(SomPos(i, j))[_plane];
          if (_in.itemAtPos(SomPos(i, j))[_plane] < Min)
            Min = _in.itemAtPos(SomPos(i, j))[_plane];
        }

    xmippFeature bw = Max - Min;




  // Calculate Planes


    for (j=0; j< _in.height(); j++)
		for (i=0; i< _in.width(); i++)
        {
        		
			if (bw != 0.0)
				_out.itemAtPos(SomPos(i, j))[0] = 0.0 + 0.9 * (_in.itemAtPos(SomPos(i, j))[_plane] - Min) / bw;
			else
				_out.itemAtPos(SomPos(i, j))[0] = 0.5;
			
			if (_in.calibrated())
				_out.targetAtPos(SomPos(i, j)) = _in.targetAtPos(SomPos(i, j));
			
        } // for i


} // getPlane


