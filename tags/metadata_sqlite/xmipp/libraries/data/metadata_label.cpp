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



//bool MDL::voidPtr2Value(const MDLabel label, void* ptrValue, MDValue &valueOut)
//{
//    switch (labelType(label))
//    {
//    case LABEL_BOOL: //bools are int in sqlite3
//        valueOut.boolValue = *((bool*)ptrValue);
//        break;
//    case LABEL_INT:
//        valueOut.intValue = *((int*)ptrValue);
//        break;
//    case LABEL_DOUBLE:
//        valueOut.doubleValue = *((double*)ptrValue);
//        break;
//    case LABEL_STRING:
//        valueOut.stringValue = *((std::string*)ptrValue);
//        break;
//    case LABEL_VECTOR:
//        valueOut.vectorValue = *((std::vector<double>*)ptrValue);
//        break;
//    }
//    return true;
//
//}

//bool MDL::value2VoidPtr(const MDLabel label, const MDValue &value, void* valuePtrOut)
//{
//    switch (labelType(label))
//    {
//    case LABEL_BOOL: //bools are int in sqlite3
//        *((bool*)valuePtrOut) = value.boolValue;
//        break;
//    case LABEL_INT:
//        *((int*)valuePtrOut) = value.intValue;
//        break;
//    case LABEL_DOUBLE:
//        *((double*)valuePtrOut) = value.doubleValue;
//        break;
//    case LABEL_STRING:
//        *((std::string*)valuePtrOut) = value.stringValue;
//        break;
//    case LABEL_VECTOR:
//        *((std::vector<double>*)valuePtrOut) = value.vectorValue;
//        break;
//    }
//    return true;
//}

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
