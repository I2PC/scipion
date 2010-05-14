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

class MetaDataContainer
{
    /** Container for pairs "name" and value. Note that void * allows to use
     mixed types */
    std::map<MDLabel, void *> values;

    void insertVoidPtr(MDLabel name, void * value);
    void * getVoidPtr(MDLabel name);

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
    void addValue(MDLabel name, const T &value)
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
    void getValue( const MDLabel name, T &value)
    {
        std::map<MDLabel, void *>::iterator element;

        element = values.find(name);

        if (element == values.end())
        {
            REPORT_ERROR(1,(std::string) "Label(int) " + MDL::label2Str(name) + " not found\n" );
        }
        else
        {
            value = *((T *) element->second);
        }
    }
    bool valueExists(MDLabel name);

    //string is not part of the template because - is not defined for string
    bool pairExists(MDLabel name, const std::string &value);

    template<class T>
    bool pairExists(MDLabel name, T value)
    {
        // Traverse all the structure looking for objects
        // that satisfy search criteria
        std::map<MDLabel, void *>::iterator It;

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



    void deleteValue(MDLabel name);

    bool writeValueToStream(std::ostream &outstream, MDLabel inputLabel);
    bool writeValueToFile(std::ofstream &outfile, MDLabel inputLabel);
    bool writeValueToString(std::string &outString, MDLabel inputLabel);
};

#endif
