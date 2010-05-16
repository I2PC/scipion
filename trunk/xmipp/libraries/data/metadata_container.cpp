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
#include "metadata_container.h"

MetaDataContainer::MetaDataContainer()
{
}
;
MetaDataContainer::~MetaDataContainer()
{
}
;

void MetaDataContainer::copy(const MetaDataContainer &MDc)
{

    if (this != &MDc)
    {
        void * aux;
        MDLabel lCode;
        std::map<MDLabel, void *>::const_iterator It;
        for (It = (MDc.values).begin(); It != (MDc.values).end(); It++)
        {
            aux = It->second;
            lCode = It->first;

            if (MDL::isDouble(lCode))
            {
                addValue(lCode, *((double *) aux));
            }
            else if (MDL::isString(lCode))
            {
                addValue(lCode, *((std::string *) aux));
            }
            else if (MDL::isInt(lCode))
            {
                addValue(lCode, *((int *) aux));
            }
            else if (MDL::isBool(lCode))
            {
                addValue(lCode, *((bool *) aux));
            }
            else if (MDL::isVector(lCode))
            {
                addValue(lCode, *((std::vector<double> *) aux));
            }

        }
    }
}

MetaDataContainer& MetaDataContainer::operator =(const MetaDataContainer &MDc)
{
    copy(MDc);
    return *this;
}

MetaDataContainer::MetaDataContainer(const MetaDataContainer &MDc)
{
    copy(MDc);
}

void MetaDataContainer::addValue(const std::string &name,
        const std::string &value)
{
    MDLabel lCode = MDL::str2Label(name);
    std::istringstream i(value);

    // Look for a double value
    if (MDL::isDouble(lCode))
    {
        double doubleValue;

        i >> doubleValue;

        addValue(lCode, doubleValue);
    }
    else if (MDL::isString(lCode))
    {
        addValue(lCode, value);
    }
    else if (MDL::isInt(lCode))
    {
        int intValue;

        i >> intValue;

        addValue(lCode, intValue);
    }
    else if (MDL::isBool(lCode))
    {
        bool boolValue;

        i >> boolValue;

        addValue(lCode, boolValue);
    }
    else if (MDL::isVector(lCode))
    {
        std::vector<double> vectorValue;
        double val;
        std::string ss;
        i >> ss;

        while (i >> val)
            vectorValue.push_back(val);
        addValue(lCode, vectorValue);
    }
}

void MetaDataContainer::insertVoidPtr(MDLabel name, void * value)
{
    values[name] = value;
}
void * MetaDataContainer::getVoidPtr(MDLabel name)
{
    std::map<MDLabel, void *>::iterator element;

    element = values.find(name);

    if (element == values.end())
    {
        REPORT_ERROR(1,(std::string) "Label " + MDL::label2Str(name) + " not found on getVoidPtr()\n" );
    }
    else
    {
        return element->second;
    }
}

bool MetaDataContainer::valueExists(MDLabel name)
{
    if (values.find(name) == values.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}
//A template exists for pairexists different from string
bool MetaDataContainer::pairExists(MDLabel name, const std::string &value)
{
    // Traverse all the structure looking for objects
    // that satisfy search criteria
    std::map<MDLabel, void *>::iterator It;

    It = values.find(name);

    if (It != values.end())
    {
        if (*((std::string *) (It->second)) == value)
        {
            return true;
        }
    }

    return false;
}

void MetaDataContainer::deleteValue(MDLabel name)
{
    values.erase(name);
}

bool MetaDataContainer::writeValueToStream(std::ostream &outstream,
        MDLabel inputLabel)
{
    if (valueExists(inputLabel))
    {
        if (MDL::isDouble(inputLabel))
        {
            double d;
            d = *((double*) (getVoidPtr(inputLabel)));
            if (d!= 0. && ABS(d) < 0.001)
                outstream << std::setw(12) << std::scientific;
            else
                outstream << std::setw(12) << std::fixed;
            outstream << d;
        }
        else if (MDL::isString(inputLabel))
            outstream << *((std::string*) (getVoidPtr(inputLabel)));
        else if (MDL::isInt(inputLabel))
        {
            outstream << std::setw(10) << std::fixed;
            outstream << *((int*) (getVoidPtr(inputLabel)));
        }
        else if (MDL::isBool(inputLabel))
            outstream << *((bool*) (getVoidPtr(inputLabel)));
        else if (MDL::isVector(inputLabel))
        {
            const std::vector<double> &myVector =
                    *((std::vector<double>*) (getVoidPtr(inputLabel)));
            int imax = myVector.size();
            outstream << "** ";
            for (int i = 0; i < imax; i++)
            {
                if (myVector[i] != 0. && ABS(myVector[i]) < 0.001)
                    outstream << std::setw(12) << std::scientific;
                else
                    outstream << std::setw(12) << std::fixed;
                outstream << myVector[i] << " ";
            }
            outstream << "**";
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool MetaDataContainer::writeValueToString(std::string &outString,
        MDLabel inLabel)
{
    std::ostringstream oss;
    outString = writeValueToStream(oss, inLabel) ? oss.str() : std::string("");
}

bool MetaDataContainer::writeValueToFile(std::ofstream &outfile,
        MDLabel inputLabel)
{
    writeValueToStream(outfile, inputLabel);
}

