/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include "multidim_array_generic.h"

MultidimArrayGeneric::MultidimArrayGeneric(MultidimArrayBase* array, DataType _datatype)
{
    im = array;
    datatype = _datatype;
    destroyData = false;
}

MultidimArrayGeneric::MultidimArrayGeneric(MultidimArrayGeneric &mdim, int select_slice)
{
    init();
    setDatatype(mdim.datatype);

#define ALIAS(type) ((MultidimArray<type>*)(im))->aliasSlice(*((MultidimArray<type>*)mdim.im), select_slice);

    SWITCHDATATYPE(mdim.datatype, ALIAS)
#undef ALIAS
}

MultidimArrayGeneric::~MultidimArrayGeneric()
{
    if (im != NULL  && destroyData)
        delete im;
}
void MultidimArrayGeneric::init()
{
    im = NULL;
    datatype = DT_Unknown;
    destroyData = true;
}

void MultidimArrayGeneric::clear()
{
    if (im != NULL && destroyData)
    {
        im->clear();
        delete im;
        init();
    }
}

void MultidimArrayGeneric::link(MultidimArrayBase* array)
{
    im = array;
    destroyData = false;
}


void MultidimArrayGeneric::setDatatype(DataType imgType)
{
    clear();
    datatype = imgType;
    destroyData = true;

    switch (datatype)
    {
    case DT_Float:
        {
            MultidimArray<float> *imT = new MultidimArray<float>;
            im = imT;
        }
        break;
    case DT_UInt:
        {
            MultidimArray<unsigned int> *imT = new MultidimArray<unsigned int>;
            im = imT;
        }
        break;
    case DT_Int:
        {
            MultidimArray<int> *imT = new MultidimArray<int>;
            im = imT;
        }
        break;
    case DT_UShort:
        {
            MultidimArray<unsigned short> *imT = new MultidimArray<unsigned short>;
            im = imT;
        }
        break;
    case DT_Short:
        {
            MultidimArray<short> *imT = new MultidimArray<short>;
            im = imT;
        }
        break;
    case DT_UHalfByte:
    case DT_UChar:
        {
            MultidimArray<unsigned char> *imT = new MultidimArray<unsigned char>;
            im = imT;
        }
        break;
    case DT_SChar:
        {
            MultidimArray<char> *imT = new MultidimArray<char>;
            im = imT;
        }
        break;
    case DT_Unknown:
        REPORT_ERROR(ERR_IMG_UNKNOWN,"");
        break;
    default:
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Datatype not implemented.");
        break;
    }
}

bool MultidimArrayGeneric::operator==(const MultidimArrayGeneric &mdA) const
{
	if (datatype != mdA.datatype)
	{
		return false;

	}

#define COMPARE(type) return ( ((MultidimArray<type>*)im)->equal(*(MultidimArray<type>*)mdA.im) );
    SWITCHDATATYPE(datatype,COMPARE)
#undef COMPARE
}


bool MultidimArrayGeneric::equal(const MultidimArrayGeneric &op,
           double accuracy) const
{
	if (datatype != op.datatype)
	{
		return false;
	}

#define COMPARE(type) return ((MultidimArray<type>*)im)->equal(*(MultidimArray<type>*)op.im,accuracy);
    SWITCHDATATYPE(datatype,COMPARE)
#undef COMPARE
}



void MultidimArrayGeneric::aliasSlice(MultidimArrayGeneric &mdim, int select_slice)
{
    setDatatype(mdim.datatype);

#define ALIAS(type) ((MultidimArray<type>*)(im))->aliasSlice(*((MultidimArray<type>*)mdim.im), select_slice);

    SWITCHDATATYPE(mdim.datatype, ALIAS)
#undef ALIAS
}
