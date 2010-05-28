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
        ss << "TEXT"; //FIXME: This is provisional
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


//----------- Implementation of MDValue -----------------
inline void MDValue::labelTypeCheck(MDLabelType type) const
{
    if (MDL::labelType(this->label) != type)
    {
        std::stringstream ss;
        ss << "Mismatch Label (" << MDL::label2Str(label) << ") and value type(";
        switch (type)
        {
        case LABEL_INT:
            ss << "int";
            break;
        case LABEL_LONG:
            ss << "long int";
            break;
        case LABEL_STRING:
            ss << "string";
            break;
        case LABEL_BOOL:
            ss << "bool";
            break;
        case LABEL_VECTOR:
            ss << "vector";
            break;
        }
        ss << ")";
        REPORT_ERROR(-55, ss.str());
    }
}

//Just a simple constructor with the label
//dont do any type checking as have not value yet
MDValue::MDValue(MDLabel label)
{
    this->label = label;
}
///Constructors for each Label supported type
///these constructor will do the labels type checking
MDValue::MDValue(MDLabel label, const int &intValue)
{
    this->label = label;
    labelTypeCheck(LABEL_INT);
    this->intValue = intValue;
}
MDValue::MDValue(MDLabel label, const double &doubleValue)
{
    this->label = label;
    labelTypeCheck(LABEL_DOUBLE);
    this->doubleValue = doubleValue;
}
MDValue::MDValue(MDLabel label, const bool &boolValue)
{
    this->label = label;
    labelTypeCheck(LABEL_BOOL);
    this->boolValue = boolValue;
}
MDValue::MDValue(MDLabel label, const std::string &stringValue)
{
    this->label = label;
    labelTypeCheck(LABEL_STRING);
    this->stringValue = stringValue;
}
MDValue::MDValue(MDLabel label, const std::vector<double> &vectorValue)
{
    this->label = label;
    labelTypeCheck(LABEL_VECTOR);
    this->vectorValue = vectorValue;
}
MDValue::MDValue(MDLabel label, const long int longintValue)
{
    this->label = label;
    labelTypeCheck(LABEL_LONG);
    this->longintValue = longintValue;

}
MDValue::MDValue(MDLabel label, const float &floatValue)
{
    std::cerr << "Do not use setValue with floats, use double"<< std::endl;
    std::cerr << "Floats are banned from metadata class"<< std::endl;
    exit(1);
}
MDValue::MDValue(MDLabel label, const char &charValue)
{
    std::cerr << "Do not use setValue with char, use string"<< std::endl;
    std::cerr << "chars are banned from metadata class"<< std::endl;
    exit(1);
}

//These getValue also do a compilation type checking
//when expanding templates functions and only
//will allow the supported types
//TODO: think if the type check if needed here
void MDValue::getValue(int &iv) const
{
    labelTypeCheck(LABEL_INT);
    iv = this->intValue;
}
void MDValue::getValue(double &dv) const
{
    labelTypeCheck(LABEL_DOUBLE);
    dv = this->doubleValue;
}
void MDValue::getValue(bool &bv) const
{
    labelTypeCheck(LABEL_BOOL);
    bv = this->boolValue;
}
void MDValue::getValue(std::string &sv) const
{
    labelTypeCheck(LABEL_STRING);
    sv = this->stringValue;
}
void  MDValue::getValue(std::vector<double> &vv) const
{
    labelTypeCheck(LABEL_VECTOR);
    vv = this->vectorValue;
}
void MDValue::getValue(long int &lv) const
{
    labelTypeCheck(LABEL_INT);
    lv = this->longintValue;
}
void MDValue::getValue(float &floatvalue) const
{
    std::cerr << "Do not use setValue with floats, use double"<< std::endl;
    std::cerr << "Floats are banned from metadata class"<< std::endl;
    exit(1);
}
void MDValue::getValue(char  &charvalue) const
{
    std::cerr << "Do not use setValue with char, use string"<< std::endl;
    std::cerr << "chars are banned from metadata class"<< std::endl;
    exit(1);
}

#define DOUBLE2STREAM(d) \
    if (withFormat) {\
            (os) << std::setw(12); \
            (os) << (((d) != 0. && ABS(d) < 0.001) ? std::scientific : std::fixed);\
        } os << d;

#define INT2STREAM(i) \
    if (withFormat) os << std::setw(10); \
    os << i;

void MDValue::toStream(std::ostream &os, bool withFormat) const
{
    switch (MDL::labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
        os << boolValue;
        break;
    case LABEL_INT:
        INT2STREAM(intValue);
        break;
    case LABEL_LONG:
        INT2STREAM(longintValue);
        break;
    case LABEL_DOUBLE:
        DOUBLE2STREAM(doubleValue);
        break;
    case LABEL_STRING:
        os << stringValue;
        break;
    case LABEL_VECTOR:
        os << "** ";
        int size = vectorValue.size();
        for (int i = 0; i < size; i++)
        {
            DOUBLE2STREAM(vectorValue[i]);
            os << " ";
        }
        os << "**";
        break;
    }//close switch
}//close function toStream

std::string MDValue::toString(bool withFormat) const
{
    std::stringstream ss;
    toStream(ss, withFormat);
    return ss.str();
}

bool MDValue::fromStream(std::istream &is)
{
    switch (MDL::labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
        is >> boolValue;
        break;
    case LABEL_INT:
        is >> intValue;
        break;
    case LABEL_LONG:
        is >> longintValue;
        break;
    case LABEL_DOUBLE:
        is >> doubleValue;
        break;
    case LABEL_STRING:
        is >> stringValue;
        break;
    case LABEL_VECTOR:
        is >> stringValue; //Just to remove "**"
        while (is >> doubleValue) //This will stop at ending "**"
            vectorValue.push_back(doubleValue);
        break;
    }
    return true;
}

bool MDValue::fromString(const std::string &str)
{
    std::stringstream ss(str);
    fromStream(ss);
}

