/***************************************************************************
 *
 * Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#ifndef METADATACONTAINER_H
#define METADATACONTAINER_H

#include <map>
#include "strings.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "funcs.h"
#include "metadata_label.h"

//useful to init values to zero
static double zeroD=0.;
static double    oneD=1.;
static bool  falseb=false;

//FIXME: For now keeping compatibility
typedef MDLabel MetaDataLabel;

inline bool isString(MetaDataLabel lCode)
{
    return MDL::isString(lCode);
}

inline bool isDouble(MetaDataLabel lCode)
{
    return MDL::isDouble(lCode);
}

inline bool isVector(MetaDataLabel lCode)
{
   return MDL::isVector(lCode);
}

inline bool isBool(MetaDataLabel lCode)
{
    return MDL::isBool(lCode);
}

inline bool isInt(MetaDataLabel lCode)
{
    return MDL::isInt(lCode);
}

class MetaDataContainer
{
    /** Container for pairs "name" and value. Note that void * allows to use
     mixed types */
    std::map<MetaDataLabel, void *> values;

    void insertVoidPtr(MetaDataLabel name, void * value);
    void * getVoidPtr(MetaDataLabel name);

public:

    /**Assignment operator
     *
     */
    MetaDataContainer& operator =(const MetaDataContainer &MDc);

    /** Constructor */
    MetaDataContainer();
    /** Copy constructor
     *
     */
    MetaDataContainer(const MetaDataContainer &MDc);

    /** Destructor */
    ~MetaDataContainer();

    /** Create a new pair name-value of integer type */
    void addValue(const std::string &name, const std::string &value);

    template<class T>
    void addValue(MetaDataLabel name, const T &value)
    {
        void * newValue = (void *) (new T(value));
        insertVoidPtr(name, newValue);
    }

    /** clean metadatacontainer
     *
     */
    void clear(void)
    {
    	values.clear();
    }
    template<class T>
    void getValue( const MetaDataLabel name, T &value)
    {
        std::map<MetaDataLabel, void *>::iterator element;

        element = values.find(name);

        if (element == values.end())
        {
            REPORT_ERROR(1,(std::string) "Label(int) " + decodeLabel(name) + " not found\n" );
        }
        else
        {
            value = *((T *) element->second);
        }
    }
    bool valueExists(MetaDataLabel name);

    //string is not part of the template because - is not defined for string
    bool pairExists(MetaDataLabel name, const std::string &value);

    template<class T>
    bool pairExists(MetaDataLabel name, T value)
    {
        // Traverse all the structure looking for objects
        // that satisfy search criteria
        std::map<MetaDataLabel, void *>::iterator It;

        It = values.find(name);

        if (It != values.end())
        {
            if (ABS( *((T *)(It->second)) - value )
                    < XMIPP_EQUAL_ACCURACY)
            {
                return true;
            }
        }

        return false;
    }



    void deleteValue(MetaDataLabel name);

    bool writeValueToStream(std::ostream &outstream, MetaDataLabel inputLabel);
    bool writeValueToFile(std::ofstream &outfile, MetaDataLabel inputLabel);
    bool writeValueToString(std::string &outString, MetaDataLabel inputLabel);

    static MetaDataLabel codifyLabel(std::string strLabel);
    static std::string decodeLabel(MetaDataLabel inputLabel);
    static bool isValidLabel(MetaDataLabel inputLabel);
    static bool isValidLabel(std::string inputLabel);
};

#endif
