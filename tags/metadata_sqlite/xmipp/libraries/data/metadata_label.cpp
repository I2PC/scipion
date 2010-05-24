/***************************************************************************
 *
 * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
#include "metadata_label.h"

//This is needed for static memory allocation
std::map<MDLabel, MDLabelData> MDL::data;
std::map<std::string, MDLabel> MDL::names;
MDLabelStaticInit MDL::initialization; //Just for initialization

void MDL::addLabel(const MDLabel label, const MDLabelType type, std::string name, std::string name2, std::string name3)
{
    data[label] = MDLabelData(type, name);
    names[name] = label;
    if (name2 != "")
        names[name2] = label;
    if (name3 != "")
        names[name3] = label;
}//close function addLabel

MDLabel  MDL::str2Label(const std::string &labelName)
{
    if (names.find(labelName) == names.end())
        return MDL_UNDEFINED;
    return names[labelName];
}//close function str2Label

std::string  MDL::label2Str(const MDLabel label)
{
    if (data.find(label) == data.end())
        return "";
    return data[label].str;
}//close function label2Str

bool MDL::double2Stream(double d, std::ostream &os, bool withFormat)
{
    if (withFormat)
    {
        os << std::setw(12);
        os << ((d != 0. && ABS(d) < 0.001) ? std::scientific : std::fixed);
    }
    os << d;
}

bool MDL::value2Stream(const MDLabel label, const MDValue &value, std::ostream &os, bool withFormat)
{
    switch (labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
        os << value.boolValue;
        break;
    case LABEL_INT:
        if (withFormat)
            os << std::setw(10);
        os << value.intValue;
        break;
    case LABEL_DOUBLE:
        double2Stream(value.doubleValue, os, withFormat);
        break;
    case LABEL_STRING:
        os << value.stringValue;
        break;
    case LABEL_VECTOR:
        os << "** ";
        int size = value.vectorValue.size();
        for (int i = 0; i < size; i++)
        {
            double2Stream(value.vectorValue[i], os, withFormat);
            os << " ";
        }
        os << "**";
        break;
    }
}

std::string MDL::value2Str(const MDLabel label, const MDValue &value, bool withFormat)
{
    std::stringstream ss;
    value2Stream(label, value, ss, withFormat);
    return ss.str();
}

bool MDL::str2Value(const MDLabel label, const std::string &str, MDValue &valueOut)
{
    std::stringstream ss(str);

    switch (labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
        ss >> valueOut.boolValue;
        break;
    case LABEL_INT:
        ss >> valueOut.intValue;
        break;
    case LABEL_DOUBLE:
        ss >> valueOut.doubleValue;
        break;
    case LABEL_STRING:
        valueOut.stringValue = str;
        break;
    case LABEL_VECTOR:
        REPORT_ERROR(-55, "Not yet supported");
        break;
    }
    return true;
}

bool MDL::voidPtr2Value(const MDLabel label, void* ptrValue, MDValue &valueOut)
{
    switch (labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
        valueOut.boolValue = *((bool*)ptrValue);
        break;
    case LABEL_INT:
        valueOut.intValue = *((int*)ptrValue);
        break;
    case LABEL_DOUBLE:
        valueOut.doubleValue = *((double*)ptrValue);
        break;
    case LABEL_STRING:
        valueOut.stringValue = *((std::string*)ptrValue);
        break;
    case LABEL_VECTOR:
        valueOut.vectorValue = *((std::vector<double>*)ptrValue);
        break;
    }
    return true;

}

std::string MDL::label2SqlColumn(const MDLabel label)
{
    std::stringstream ss;
    ss << MDL::label2Str(label) << " ";
    switch (MDL::labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
    case LABEL_INT:
        ss << "INTEGER";
        break;
    case LABEL_DOUBLE:
        ss << "REAL";
        break;
    case LABEL_STRING:
        ss << "TEXT";
        break;
    case LABEL_VECTOR:
        REPORT_ERROR(-55, "Metadata still not suport vectors");
        break;
    }
    return ss.str();
}

bool MDL::isInt(const MDLabel label)
{
    return (data[label].type == LABEL_INT);
}
bool MDL::isBool(const MDLabel label)
{
    return (data[label].type == LABEL_BOOL);
}
bool MDL::isString(const MDLabel label)
{
    return (data[label].type == LABEL_STRING);
}
bool MDL::isDouble(const MDLabel label)
{
    return (data[label].type == LABEL_DOUBLE);
}
bool MDL::isVector(const MDLabel label)
{
    return (data[label].type == LABEL_VECTOR);
}

bool MDL::isValidLabel(const MDLabel label)
{
    return (label > MDL_UNDEFINED && label < MDL_LAST_LABEL);
}

bool MDL::isValidLabel(const std::string &labelName)
{
    return isValidLabel(str2Label(labelName));
}

MDLabelType MDL::labelType(const MDLabel label)
{
    return data[label].type;
}

MDLabelType MDL::labelType(std::string &labelName)
{
    return data[str2Label(labelName)].type;
}
