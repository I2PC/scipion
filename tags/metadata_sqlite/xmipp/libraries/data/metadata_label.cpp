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
StaticInitialization MDL::initialization; //Just for initialization

void MDL::addLabel(MDLabel label, MDLabelType type, std::string name, std::string name2, std::string name3)
{
    data[label] = MDLabelData(type, name);
    names[name] = label;
    if (name2 != "")
        names[name2] = label;
    if (name3 != "")
        names[name3] = label;
}//close function addLable

MDLabel  MDL::str2Label(const std::string &labelName)
{
    if (names.find(labelName) == names.end())
        return MDL_UNDEFINED;
    return names[labelName];
}//close function str2Label

std::string  MDL::label2Str(const MDLabel &label)
{
    if (data.find(label) == data.end())
            return "";
    return data[label].str;
}//close function label2Str

bool MDL::isInt(const MDLabel &label)
{
    return (data[label].type == LABEL_INT);
}
bool MDL::isBool(const MDLabel &label)
{
    return (data[label].type == LABEL_BOOL);
}
bool MDL::isString(const MDLabel &label)
{
    return (data[label].type == LABEL_STRING);
}
bool MDL::isDouble(const MDLabel &label)
{
    return (data[label].type == LABEL_DOUBLE);
}
bool MDL::isVector(const MDLabel &label)
{
    return (data[label].type == LABEL_VECTOR);
}

bool MDL::isValidLabel(const MDLabel &label)
{
    return (label > MDL_UNDEFINED && label < MDL_LAST_LABEL);
}

bool MDL::isValidLabel(const std::string &labelName)
{
    MDLabel label = MDL::str2Label(labelName);
    return MDL::isValidLabel(label);
}
