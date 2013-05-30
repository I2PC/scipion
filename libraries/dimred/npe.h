/***************************************************************************
 *
 * Authors:    Oscar Hueros Fernandez            o.hueros@gmail.com (2013)
 *
 * Universidad San Pablo CEU
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

#ifndef _NPE
#define _NPE

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include "dimred_tools.h"

class NPE: public DimRedAlgorithm
{
public:
	int k;
public:
	/// Set specific parameters
	void setSpecificParameters(int k=12);

	/// Reduce dimensionality
	void reduceDimensionality();
};
//@}
#endif
