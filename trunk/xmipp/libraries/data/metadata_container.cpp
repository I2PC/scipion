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
MetaDataContainer& MetaDataContainer::operator =(const MetaDataContainer &MDc)
{

    if (this != &MDc)
    {
        void * aux;
        MetaDataLabel lCode;
        std::map<MetaDataLabel, void *>::const_iterator It;
        for (It = (MDc.values).begin(); It != (MDc.values).end(); It++)
        {
            aux = It->second;
            lCode = It->first;

            if (isDouble(lCode))
            {
                addValue(lCode, *((double *) aux));
            }
            else if (isString(lCode))
            {
                addValue(lCode, *((std::string *) aux));
            }
            else if (isInt(lCode))
            {
                addValue(lCode, *((int *) aux));
            }
            else if (isBool(lCode))
            {
                addValue(lCode, *((bool *) aux));
            }
            else if (isVector(lCode))
            {
                addValue(lCode, *((std::vector<double> *) aux));
            }

        }
    }
    return *this;
}
MetaDataContainer::MetaDataContainer(const MetaDataContainer &MDc)
{
    void * aux;
    MetaDataLabel lCode;
    std::map<MetaDataLabel, void *>::const_iterator It;
    for (It = (MDc.values).begin(); It != (MDc.values).end(); It++)
    {
        aux = It->second;
        lCode = It->first;

        if (isDouble(lCode))
        {
            addValue(lCode, *((double *) aux));
        }
        else if (isString(lCode))
        {
            addValue(lCode, *((std::string *) aux));
        }
        else if (isInt(lCode))
        {
            addValue(lCode, *((int *) aux));
        }
        else if (isBool(lCode))
        {
            addValue(lCode, *((bool *) aux));
        }
        else if (isVector(lCode))
        {
            addValue(lCode, *((std::vector<double> *) aux));
        }
    }
}

void MetaDataContainer::addValue(const std::string &name,
        const std::string &value)
{
    MetaDataLabel lCode = codifyLabel(name);
    std::istringstream i(value);

    // Look for a double value
    if (isDouble(lCode))
    {
        double doubleValue;

        i >> doubleValue;

        addValue(lCode, doubleValue);
    }
    else if (isString(lCode))
    {
        addValue(lCode, value);
    }
    else if (isInt(lCode))
    {
        int intValue;

        i >> intValue;

        addValue(lCode, intValue);
    }
    else if (isBool(lCode))
    {
        bool boolValue;

        i >> boolValue;

        addValue(lCode, boolValue);
    }
    else if (isVector(lCode))
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

void MetaDataContainer::insertVoidPtr(MetaDataLabel name, void * value)
{
    values[name] = value;
}
void * MetaDataContainer::getVoidPtr(MetaDataLabel name)
{
    std::map<MetaDataLabel, void *>::iterator element;

    element = values.find(name);

    if (element == values.end())
    {
        REPORT_ERROR(1,(std::string) "Label " + decodeLabel(name) + " not found on getVoidPtr()\n" );
    }
    else
    {
        return element->second;
    }
}

bool MetaDataContainer::valueExists(MetaDataLabel name)
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
bool MetaDataContainer::pairExists(MetaDataLabel name, const std::string &value)
{
    // Traverse all the structure looking for objects
    // that satisfy search criteria
    std::map<MetaDataLabel, void *>::iterator It;

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

void MetaDataContainer::deleteValue(MetaDataLabel name)
{
    values.erase(name);
}

MetaDataLabel MetaDataContainer::codifyLabel(std::string strLabel)
{
    return MDL::str2Label(strLabel);
}

std::string MetaDataContainer::decodeLabel(MetaDataLabel inputLabel)
{
    return MDL::label2Str(inputLabel);
}

bool MetaDataContainer::writeValueToStream(std::ostream &outstream,
        MetaDataLabel inputLabel)
{
    if (valueExists(inputLabel))
    {
        if (isDouble(inputLabel))
        {
            double d;
            d = *((double*) (getVoidPtr(inputLabel)));
            if (d!= 0. && ABS(d) < 0.001)
                outstream << std::setw(12) << std::scientific;
            else
                outstream << std::setw(12) << std::fixed;
            outstream << d;
        }
        else if (isString(inputLabel))
            outstream << *((std::string*) (getVoidPtr(inputLabel)));
        else if (isInt(inputLabel))
        {
            outstream << std::setw(10) << std::fixed;
            outstream << *((int*) (getVoidPtr(inputLabel)));
        }
        else if (isBool(inputLabel))
            outstream << *((bool*) (getVoidPtr(inputLabel)));
        else if (isVector(inputLabel))
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
        MetaDataLabel inLabel)
{
    std::ostringstream oss;
    outString = writeValueToStream(oss, inLabel) ? oss.str() : std::string("");
}

bool MetaDataContainer::writeValueToFile(std::ofstream &outfile,
        MetaDataLabel inputLabel)
{
    writeValueToStream(outfile, inputLabel);
}

bool MetaDataContainer::isValidLabel(MetaDataLabel inputLabel)
{
    return MDL::isValidLabel(inputLabel);
}

bool MetaDataContainer::isValidLabel(std::string inputLabel)
{
    return MDL::isValidLabel(inputLabel);
}
