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

#include <stdlib.h>
#include <algorithm>
#include "metadata_label.h"

//This is needed for static memory allocation
//std::map<MDLabel, MDLabelData> MDL::data;
MDLabelData * MDL::data[MDL_LAST_LABEL];
std::map<std::string, MDLabel> MDL::names;
MDRow MDL::emptyHeader;
MDLabelStaticInit MDL::initialization; //Just for initialization


void MDL::addLabel(const MDLabel label, const MDLabelType type, const String &name, int tags)
{
    data[(int)label] = new MDLabelData(type, name, tags);
    names[name] = label;
}//close function addLabel

/**
 * Extra alias can be defined through the environment var XMIPP_EXTRA_ALIASES
 * the syntax is the following
 * XMIPP_EXTRA_ALIASES='anglePsi=otherAnglePsi;shiftX=otherShiftX;shiftY:otherShift'
 * The = sign will add a label, and the alias name will replace the current one (replace=True)
 * if the : sign is used, only a normal alias will added
 */
void MDL::addExtraAliases()
{
  const char * extra_aliases = getenv("XMIPP_EXTRA_ALIASES");

  if (extra_aliases)
  {
      StringVector sv, pair;
      String eq = "=", co = ":";
      tokenize(extra_aliases, sv, ";");
      MDLabel label;
      bool replace;

      for (std::vector<String>::iterator it = sv.begin(); it != sv.end(); ++it)
      {
          if (it->find(eq) != it->npos)
          {
              tokenize(*it, pair, "=");
              replace = true;
          }
          else if (it->find(co) != it->npos)
          {
              tokenize(*it, pair, co);
              replace = false;
          }
          else
              REPORT_ERROR(ERR_ARG_INCORRECT, "Invalid pair separator, use = or :");
          label = MDL::str2Label(pair[0]);
          // Add the label alias
          if (label != MDL_UNDEFINED)
              addLabelAlias(label, pair[1], replace);
          else
              REPORT_ERROR(ERR_ARG_INCORRECT, formatString("Invalid label name: %s found in enviroment var XMIPP_EXTRA_ALIASES", pair[0].c_str()));
      }
  }
}//close function addLabel

void MDL::addLabelAlias(const MDLabel label, const String &alias, bool replace)
{
    names[alias] = label;
    if (replace)
    {
      data[(int)label]->str = alias;
    }
}//close function addLabel

void MDL::str2LabelVector(const String &labelsStr, std::vector<MDLabel> &labels)
{
    labels.clear();
    StringVector parts;
    splitString(labelsStr, " ", parts);
    for (size_t i = 0; i < parts.size(); ++i)
        if (MDL::isValidLabel(parts[i]))
            labels.push_back(MDL::str2Label(parts[i]));
        else
            REPORT_ERROR(ERR_PARAM_INCORRECT, formatString("Unknown label '%s' received.", parts[i].c_str()));
}

MDLabel  MDL::str2Label(const String &labelName)
{
    if (names.find(labelName) == names.end())
        return MDL_UNDEFINED;
    return names[labelName];
}//close function str2Label

String  MDL::label2Str(const MDLabel label)
{
    return  (isValidLabel(label)) ? data[(int)label]->str : "";
}//close function label2Str

String MDL::label2SqlColumn(const MDLabel label)
{
    std::stringstream ss;
    ss << MDL::label2Str(label) << " ";
    switch (MDL::labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
    case LABEL_INT:
    case LABEL_SIZET:
        ss << "INTEGER";
        break;
    case LABEL_DOUBLE:
        ss << "REAL";
        break;
    case LABEL_STRING:
        ss << "TEXT";
        break;
    case LABEL_VECTOR_DOUBLE:
    case LABEL_VECTOR_SIZET:
        ss << "TEXT";
        break;
    case LABEL_NOTYPE:
    	ss << "NO_TYPE";
    	break;
    }
    return ss.str();
}

String MDL::labelType2Str(MDLabelType type)
{
    switch (type)
    {
    case LABEL_STRING:
        return "STRING";
    case LABEL_DOUBLE:
        return "DOUBLE";
    case LABEL_INT:
        return "INT";
    case LABEL_BOOL:
        return "BOOL";
    case LABEL_VECTOR_DOUBLE:
        return "VECTOR(DOUBLE)";
    case LABEL_SIZET:
        return "SIZE_T";
    case LABEL_VECTOR_SIZET:
        return "VECTOR(SIZE_T)";
    case LABEL_NOTYPE:
    	return "NO_TYPE";
    }
    return "UNKNOWN";
}

bool MDL::isInt(const MDLabel label)
{
    return (data[(int)label]->type == LABEL_INT);
}
bool MDL::isLong(const MDLabel label)
{
    return (data[(int)label]->type == LABEL_SIZET);
}
bool MDL::isBool(const MDLabel label)
{
    return (data[(int)label]->type == LABEL_BOOL);
}
bool MDL::isString(const MDLabel label)
{
    return (data[(int)label]->type == LABEL_STRING);
}
bool MDL::isDouble(const MDLabel label)
{
    return (data[(int)label]->type == LABEL_DOUBLE);
}
bool MDL::isVector(const MDLabel label)
{
    return (data[(int)label]->type == LABEL_VECTOR_DOUBLE);
}
bool MDL::isVectorLong(const MDLabel label)
{
    return (data[(int)label]->type == LABEL_VECTOR_SIZET);
}
bool MDL::isValidLabel(const MDLabel label)
{
    return label > MDL_UNDEFINED &&
           label < MDL_LAST_LABEL &&
           data[(int)label] != NULL;
}

bool MDL::isValidLabel(const String &labelName)
{
    return isValidLabel(str2Label(labelName));
}

MDLabelType MDL::labelType(const MDLabel label)
{
    return data[(int)label]->type;
}

MDLabelType MDL::labelType(const String &labelName)
{
    return data[str2Label(labelName)]->type;
}

std::map<String, MDLabel>& MDL::getLabelDict()
{
    return names;
}

bool MDL::hasTag(const MDLabel label, const int tags)
{
    return data[(int)label]->tags & tags;
}

bool MDL::isTextFile(const MDLabel label)
{
    return data[(int)label]->tags & TAGLABEL_TEXTFILE;
}

bool MDL::isMetadata(const MDLabel label)
{
    return data[(int)label]->tags & TAGLABEL_METADATA;
}

bool MDL::isCtfParam(const MDLabel label)
{
    return data[(int)label]->tags & TAGLABEL_CTFPARAM;
}

bool MDL::isImage(const MDLabel label)
{
    return data[(int)label]->tags & TAGLABEL_IMAGE;
}

bool MDL::isStack(const MDLabel label)
{
    return data[(int)label]->tags & TAGLABEL_STACK;
}

bool MDL::isMicrograph(const MDLabel label)
{
    return data[(int)label]->tags & TAGLABEL_MICROGRAPH;
}

bool MDL::isPSD(const MDLabel label)
{
    return data[(int)label]->tags & TAGLABEL_PSD;
}

bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label)
{
    std::vector<MDLabel>::const_iterator location;
    location = std::find(labelsVector.begin(), labelsVector.end(), label);

    return (location != labelsVector.end());
}

void MDObject::copy(const MDObject &obj)
{
    label = obj.label;
    type = obj.type;
    chr = obj.chr;
    if (type == LABEL_STRING)
    {
        delete data.stringValue;
        data.stringValue = new String(*(obj.data.stringValue));
    }
    else if (type == LABEL_VECTOR_DOUBLE)
    {
        delete data.vectorValue;
        data.vectorValue = new std::vector<double>(*(obj.data.vectorValue));
    }
    else if (type == LABEL_VECTOR_SIZET)
    {
        delete data.vectorValueLong;
        data.vectorValueLong = new std::vector<size_t>(*(obj.data.vectorValueLong));
    }
    else
        data = obj.data;
}

MDObject::MDObject(const MDObject & obj)
{
    data.doubleValue = 0;
    copy(obj);
}
MDObject & MDObject::operator = (const MDObject &obj)
{
    data.doubleValue = 0;
    copy(obj);
    return *this;
}

//----------- Implementation of MDValue -----------------
inline void MDObject::labelTypeCheck(MDLabelType checkingType) const
{
    if (this->type != checkingType)
    {
        std::stringstream ss;
        ss << "Mismatch Label (" << MDL::label2Str(label)
        << ") and value type(" << MDL::labelType2Str(checkingType) << ")";
        REPORT_ERROR(ERR_MD_BADLABEL, ss.str());
    }
}

//Just a simple constructor with the label
//dont do any type checking as have not value yet
MDObject::MDObject(MDLabel label)
{
    this->label = label;
    chr = _SPACE;
    if (label != MDL_UNDEFINED)
    {
        type = MDL::labelType(label);
        if (type == LABEL_STRING)
            data.stringValue = new String;
        else if (type == LABEL_VECTOR_DOUBLE)
            data.vectorValue = new std::vector<double>;
        else if (type == LABEL_VECTOR_SIZET)
            data.vectorValueLong = new std::vector<size_t>;
    }
    else
        type = LABEL_NOTYPE;
}
///Constructors for each Label supported type
///these constructor will do the labels type checking
MDObject::MDObject(MDLabel label, const int &intValue)
{
    this->label = label;
    this->type = MDL::labelType(label);
    labelTypeCheck(LABEL_INT);
    this->data.intValue = intValue;
}
MDObject::MDObject(MDLabel label, const double &doubleValue)
{
    this->label = label;
    this->type = MDL::labelType(label);
    labelTypeCheck(LABEL_DOUBLE);
    this->data.doubleValue = doubleValue;
}
MDObject::MDObject(MDLabel label, const bool &boolValue)
{
    this->label = label;
    this->type = MDL::labelType(label);
    labelTypeCheck(LABEL_BOOL);
    this->data.boolValue = boolValue;
}
MDObject::MDObject(MDLabel label, const String &stringValue)
{
    this->label = label;
    this->type = MDL::labelType(label);
    labelTypeCheck(LABEL_STRING);
    this->data.stringValue = new String(stringValue);
}
MDObject::MDObject(MDLabel label, const std::vector<double> &vectorValue)
{
    this->label = label;
    this->type = MDL::labelType(label);
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    this->data.vectorValue = new std::vector<double>(vectorValue);
}
MDObject::MDObject(MDLabel label, const std::vector<size_t> &vectorValueLong)
{
    this->label = label;
    this->type = MDL::labelType(label);
    labelTypeCheck(LABEL_VECTOR_SIZET);
    this->data.vectorValueLong = new std::vector<size_t>(vectorValueLong);
}
MDObject::MDObject(MDLabel label, const size_t &longintValue)
{
    this->label = label;
    this->type = MDL::labelType(label);
    labelTypeCheck(LABEL_SIZET);
    this->data.longintValue = longintValue;

}
MDObject::MDObject(MDLabel label, const float &floatValue)
{
    std::cerr << "Do not use setValue with floats, use double"<< std::endl;
    std::cerr << "Floats are banned from metadata class"<< std::endl;
    exit(1);
}
MDObject::MDObject(MDLabel label, const char* &charValue)
{
    std::cerr << "Do not use setValue with char, use string"<< std::endl;
    std::cerr << "chars are banned from metadata class"<< std::endl;
    exit(1);
}

MDObject::~MDObject()
{
    if (type == LABEL_STRING)
        delete data.stringValue;
    else if (type == LABEL_VECTOR_DOUBLE)
        delete data.vectorValue;
    else if (type == LABEL_VECTOR_SIZET)
        delete data.vectorValueLong;
}

//These getValue also do a compilation type checking
//when expanding templates functions and only
//will allow the supported types
//TODO: think if the type check if needed here
void MDObject::getValue(int &iv) const
{
    labelTypeCheck(LABEL_INT);
    iv = this->data.intValue;
}
void MDObject::getValue(double &dv) const
{
    labelTypeCheck(LABEL_DOUBLE);
    dv = this->data.doubleValue;
}
void MDObject::getValue(bool &bv) const
{
    labelTypeCheck(LABEL_BOOL);
    bv = this->data.boolValue;
}
void MDObject::getValue(String &sv) const
{
    labelTypeCheck(LABEL_STRING);
    sv = *(this->data.stringValue);
}
void  MDObject::getValue(std::vector<double> &vv) const
{
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    vv = *(this->data.vectorValue);
}
void  MDObject::getValue(std::vector<size_t> &vv) const
{
    labelTypeCheck(LABEL_VECTOR_SIZET);
    vv = *(this->data.vectorValueLong);
}
void MDObject::getValue(size_t &lv) const
{
    labelTypeCheck(LABEL_SIZET);
    lv = this->data.longintValue;
}
void MDObject::getValue(float &floatvalue) const
{
    std::cerr << "Do not use setValue with floats, use double"<< std::endl;
    std::cerr << "Floats are banned from metadata class"<< std::endl;
    exit(1);
}
void MDObject::getValue(char*  &charvalue) const
{
    std::cerr << "Do not use setValue with char, use string"<< std::endl;
    std::cerr << "chars are banned from metadata class"<< std::endl;
    exit(1);
}

void MDObject::setValue(const int &iv)
{
    labelTypeCheck(LABEL_INT);
    this->data.intValue = iv;
}
void MDObject::setValue(const double &dv)
{
    labelTypeCheck(LABEL_DOUBLE);
    this->data.doubleValue = dv;
}

void MDObject::setValue(const bool &bv)
{
    labelTypeCheck(LABEL_BOOL);
    this->data.boolValue = bv;
}

void MDObject::setValue(const String &sv)
{
    labelTypeCheck(LABEL_STRING);
    *(this->data.stringValue) = sv;
}
void  MDObject::setValue(const std::vector<double> &vv)
{
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    *(this->data.vectorValue) = vv;
}
void  MDObject::setValue(const std::vector<size_t> &vv)
{
    labelTypeCheck(LABEL_VECTOR_SIZET);
    *(this->data.vectorValueLong) = vv;
}
void MDObject::setValue(const size_t &lv)
{
    labelTypeCheck(LABEL_SIZET);
    this->data.longintValue = lv;
}
void MDObject::setValue(const float &floatvalue)
{
    setValue((double) floatvalue);
}
void MDObject::setValue(const char*  &charvalue)
{
    setValue(String(charvalue));
}

void MDObject::toStream(std::ostream &os, bool withFormat, bool isSql, bool escape) const
{
    if (label == MDL_UNDEFINED) //if undefine label, store as a literal string
        os << data.stringValue;
    else
        switch (MDL::labelType(label))
        {
        case LABEL_BOOL: //bools are int in sqlite3
            os << data.boolValue;
            break;
        case LABEL_INT:
            INT2STREAM(data.intValue);
            break;
        case LABEL_SIZET:
            INT2STREAM(data.longintValue);
            break;
        case LABEL_DOUBLE:
            DOUBLE2STREAM(data.doubleValue);
            break;
        case LABEL_STRING:
            {
                char c = _SPACE;
                if (escape)
                {
                    if (isSql || data.stringValue->find_first_of(_DQUOT) != String::npos)
                        c = _QUOT;
                    else if (data.stringValue->find_first_of(_QUOT) != String::npos)
                        c = _DQUOT;
                    else if (data.stringValue->find_first_of(_SPACE) != String::npos)
                        c = _QUOT;
                    else if (data.stringValue->empty())
                        c = _QUOT;
                }
                if (c == _SPACE)
                    os << *(data.stringValue);
                else
                    os << c << *(data.stringValue) << c;
            }
            break;
        case LABEL_VECTOR_DOUBLE:
            {
                std::vector<double> &vectorDouble = *(data.vectorValue);
                if (escape)
                    os << _QUOT << " ";
                size_t size = vectorDouble.size();
                for (size_t i = 0; i < size; i++)
                {
                    double v = vectorDouble[i];
                    DOUBLE2STREAM(v);
                    os << " ";
                }
                if (escape)
                    os << _QUOT;
            }
            break;
        case LABEL_VECTOR_SIZET:
            {
                std::vector<size_t> &vector = *(data.vectorValueLong);
                if (escape)
                    os << _QUOT << " ";
                size_t size = vector.size();
                for (size_t i = 0; i < size; i++)
                    os << vector[i] << " ";
                if (escape)
                    os << _QUOT;
            }
            break;
        case LABEL_NOTYPE:
        	if (escape) os << _QUOT;
        	os << "No type";
        	if (escape) os << _QUOT;
        	break;
        }//close switch
}//close function toStream

String MDObject::toString(bool withFormat, bool isSql) const
{
    if (type == LABEL_STRING)
    {
        return isSql ? formatString("'%s'", data.stringValue->c_str()) : *data.stringValue;
    }
    std::stringstream ss;
    toStream(ss, withFormat, isSql, isSql);

    return ss.str();
}

//bool MDValue::fromStream(std::istream &is)
std::ostream& operator<< (std::ostream& os, MDObject &value)
{
    value.toStream(os);
    return os;
}

//bool MDValue::fromStream(std::istream &is)
std::istream& operator>> (std::istream& is, MDObject &value)
{
    value.fromStream(is);
    return is;
}

bool MDObject::fromStream(std::istream &is, bool fromString)
{
    if (label == MDL_UNDEFINED) //if undefine label, store as a literal string
    {
        String s;
        is >> s;
    }
    else
    {
        //NOTE: int, bool and long(size_t) are read as double for compatibility with old doc files
        double d;
        size_t value;
        switch (type)
        {
        case LABEL_BOOL: //bools are int in sqlite3
            is >> d;
            data.boolValue = (bool) ((int)d);
            break;
        case LABEL_INT:
            is >> d;
            data.intValue = (int) d;
            break;
        case LABEL_SIZET:
            is >> d;
            data.longintValue = (size_t) d;
            break;
        case LABEL_DOUBLE:
            is >> data.doubleValue;
            break;
        case LABEL_STRING:
            {
                data.stringValue->clear();
                String s;
                is >> s;
                char chr = s[0];
                if (chr == _QUOT || chr == _DQUOT)
                {
                    s = s.substr(1, s.size() - 1); //remove first char '
                    while (s.find_last_of(chr) == String::npos)
                    {
                        data.stringValue->append(s + " ");
                        is >> s;
                    }
                    s = s.substr(0, s.size() - 1); //remove last char '
                }
                data.stringValue->append(s);
            }
            break;
        case LABEL_VECTOR_DOUBLE:
            if (!fromString)
                is.ignore(256, _QUOT);
            //if (data.vectorValue == NULL)
            //  data.vectorValue = new std::vector<double>;
            data.vectorValue->clear();
            while (is >> d) //This will stop at ending "]"
                data.vectorValue->push_back(d);
            if (!fromString)
            {
                is.clear(); //this is for clear the fail state after found ']'
                is.ignore(256, _QUOT); //ignore the ending ']'
            }
            break;
        case LABEL_VECTOR_SIZET:
            if (!fromString)
                is.ignore(256, _QUOT);
            //if (data.vectorValue == NULL)
            //  data.vectorValue = new std::vector<double>;
            data.vectorValueLong->clear();
            while (is >> value) //This will stop at ending "]"
                data.vectorValueLong->push_back(value);
            if (!fromString)
            {
                is.clear(); //this is for clear the fail state after found ']'
                is.ignore(256, _QUOT); //ignore the ending ']'
            }
            break;
        case LABEL_NOTYPE:
        	break;
        }
    }
    return is.good();
}

bool MDObject::fromString(const String& str)
{
    if (type == LABEL_STRING)
        *data.stringValue = str;
    std::stringstream ss(str);
    return fromStream(ss, true);
}

bool MDObject::fromChar(const char * szChar)
{
    std::stringstream ss(szChar);
    return fromStream(ss);
}
//MDObject & MDRow::operator [](MDLabel label)
//{
//    for (iterator it = begin(); it != end(); ++it)
//        if ((*it)->label == label)
//            return *(*it);
//    MDObject * pObj = new MDObject(label);
//    push_back(pObj);
//
//    return *pObj;
//}

void MDRow::clear()
{
    _size = 0;
    //Just initialize all pointers with NULL value
    FOR_ALL_LABELS()
    {
        delete objects[_label];
        objects[_label] = NULL;
    }
}

bool MDRow::empty() const
{
    return _size == 0;
}

void MDRow::resetGeo(bool addLabels)
{
    setValue(MDL_ORIGIN_X,  0., addLabels);
    setValue(MDL_ORIGIN_Y,  0., addLabels);
    setValue(MDL_ORIGIN_Z,  0., addLabels);
    setValue(MDL_SHIFT_X,   0., addLabels);
    setValue(MDL_SHIFT_Y,   0., addLabels);
    setValue(MDL_SHIFT_Z,   0., addLabels);
//    setValue(MDL_ANGLE_ROT, 0., addLabels);
//    setValue(MDL_ANGLE_TILT,0., addLabels);
    setValue(MDL_ANGLE_PSI, 0., addLabels);
    setValue(MDL_WEIGHT,   1., addLabels);
    setValue(MDL_FLIP,     false, addLabels);
    setValue(MDL_SCALE,    1., addLabels);
}

int MDRow::size() const
{
    return _size;
}

bool MDRow::containsLabel(MDLabel label) const
{
    return objects[label] != NULL;
}


void MDRow::addLabel(MDLabel label)
{
    if (objects[label] == NULL)
    {
        objects[label] = new MDObject(label);
        order[_size] = label;
        ++_size;
    }
}

MDObject * MDRow::getObject(MDLabel label) const
{
    return objects[label];
}

/** Get value */
bool MDRow::getValue(MDObject &object) const
{
    int _label = object.label;
    if (objects[_label] == NULL)
        return false;
    object.copy(*(objects[_label]));
    return true;
}

/** Usefull macro for copy values */
/** Set value */
void MDRow::setValue(const MDObject &object)
{
    int _label = object.label;
    if (objects[_label] == NULL)
    {
        objects[_label] = new MDObject(object);
        order[_size] = object.label;
        ++_size;
    }
    else
        objects[_label]->copy(object);
}

void MDRow::setValueFromStr(MDLabel label, const String &value)
{
    MDObject mdValue(label);
    mdValue.fromString(value);
    setValue(mdValue);
}

MDRow::~MDRow()
{
    MDObject ** ptrObjectsLabel=&(objects[0]);
    FOR_ALL_LABELS()
    {
        delete *ptrObjectsLabel;
        ++ptrObjectsLabel;
    }
}

MDRow::MDRow(const MDRow & row)
{
    _size = 0;
    //Just initialize all pointers with NULL value
    memset(objects, 0, MDL_LAST_LABEL * sizeof(size_t));
    copy(row);
}

MDRow::MDRow()
{
    _size = 0;
    //Just initialize all pointers with NULL value
    memset(objects, 0, MDL_LAST_LABEL * sizeof(size_t));
}

MDRow& MDRow::operator = (const MDRow &row)
{
    copy(row);
    return *this;
}

void MDRow::copy(const MDRow &row)
{
    //Copy existing MDObjects from row
    //and delete unexisting ones
    _size = row._size;
    MDObject ** ptrObjectsLabel=&(objects[0]);
    MDObject * const * ptrRowObjectsLabel=&(row.objects[0]);
    FOR_ALL_LABELS()
    {
        if (*ptrRowObjectsLabel == NULL)
        {
            delete *ptrObjectsLabel;
            *ptrObjectsLabel = NULL;
        }
        else
        {
            if (*ptrObjectsLabel == NULL)
                *ptrObjectsLabel = new MDObject(*(*ptrRowObjectsLabel));
            else
                (*ptrObjectsLabel)->copy(*(*ptrRowObjectsLabel));
        }
        ++ptrObjectsLabel;
        ++ptrRowObjectsLabel;
    }
    //copy the order of labels
    memcpy(order, row.order, sizeof(int)*_size);
}

std::ostream& operator << (std::ostream &out, const MDRow &row)
{
    for (int i = 0; i < row._size; ++i)
    {
        row.objects[row.order[i]]->toStream(out);
        out << " ";
    }
    return out;
}
