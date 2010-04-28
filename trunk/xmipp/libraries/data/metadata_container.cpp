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
MetaDataContainer& MetaDataContainer::operator =(MetaDataContainer &MDc)
{

    if (this != &MDc)
    {
        void * aux;
        MetaDataLabel lCode;
        std::map<MetaDataLabel, void *>::iterator It;
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
MetaDataContainer::MetaDataContainer(MetaDataContainer &MDc)
{
    void * aux;
    MetaDataLabel lCode;
    std::map<MetaDataLabel, void *>::iterator It;
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
/**
 void MetaDataContainer::addValue(MetaDataLabel name, int value)
 {
 void * newValue = (void *) (new int(value));
 insertVoidPtr(name, newValue);
 }

 void MetaDataContainer::addValue(MetaDataLabel name, double value)
 {
 void * newValue = (void *) (new double(value));
 insertVoidPtr(name, newValue);
 }

 void MetaDataContainer::addValue(MetaDataLabel name, bool value)
 {
 void * newValue = (void *) (new bool(value));
 insertVoidPtr(name, newValue);
 }

 void MetaDataContainer::addValue(MetaDataLabel name, const std::string &value)
 {
 void * newValue = (void *) (new std::string(value));
 insertVoidPtr(name, newValue);
 }

 void MetaDataContainer::addValue(MetaDataLabel name,
 const std::vector<double> &value)
 {
 void * newValue = (void *) (new std::vector<double>(value));
 insertVoidPtr(name, newValue);
 }
 */
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
        while (i >> val)
            vectorValue.push_back(val);
        addValue(lCode, vectorValue);
    }
}

void MetaDataContainer::insertVoidPtr(MetaDataLabel name, void * value)
{
    values[name] = value;
}
/**
 void MetaDataContainer::getValue(MetaDataLabel name, int &value)
 {
 std::map<MetaDataLabel, void *>::iterator element;

 element = values.find(name);

 if (element == values.end())
 {
 REPORT_ERROR(1,(std::string) "Label(int) " + decodeLabel(name) + " not found\n" );
 }
 else
 {
 value = *((int *) element->second);
 }
 }
 void MetaDataContainer::getValue(MetaDataLabel name, double &value)
 {
 std::map<MetaDataLabel, void *>::iterator element;

 element = values.find(name);

 if (element == values.end())
 {
 REPORT_ERROR(1,(std::string) "Label(double) " + decodeLabel(name) + " not found\n" );
 }
 else
 {
 value = *((double *) element->second);
 }
 }

 void MetaDataContainer::getValue(MetaDataLabel name, std::string &value)
 {
 std::map<MetaDataLabel, void *>::iterator element;

 element = values.find(name);

 if (element == values.end())
 {
 REPORT_ERROR(1,(std::string) "Label(string) " + decodeLabel(name) + " not found\n" );
 }
 else
 {
 value = *((std::string *) element->second);
 }
 }

 void MetaDataContainer::getValue(MetaDataLabel name, bool &value)
 {
 std::map<MetaDataLabel, void *>::iterator element;

 element = values.find(name);

 if (element == values.end())
 {
 REPORT_ERROR(1,(std::string) "Label(bool) " + decodeLabel(name) + " not found\n" );
 }
 else
 {
 value = *((bool *) element->second);
 }
 }

 void MetaDataContainer::getValue(MetaDataLabel name, std::vector<double> &value)
 {
 std::map<MetaDataLabel, void *>::iterator element;

 element = values.find(name);

 if (element == values.end())
 {
 REPORT_ERROR(1,(std::string) "Label(vector) " + decodeLabel(name) + " not found\n" );
 }
 else
 {
 value = *((std::vector<double> *) element->second);
 }
 }
 */
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
    // NOTE: For the "||" cases, the left MetaDataLabel is the XMIPP_3 version's and
    // the right value is the old-styled one

    if (strLabel == "angleRot" || strLabel == "rot")
    {
        return MDL_ANGLEROT;
    }
    else if (strLabel == "angleTilt" || strLabel == "tilt")
    {
        return MDL_ANGLETILT;
    }
    else if (strLabel == "anglePsi" || strLabel == "psi")
    {
        return MDL_ANGLEPSI;
    }
    else if (strLabel == "angleRot2" || strLabel == "rot2")
    {
        return MDL_ANGLEROT2;
    }
    else if (strLabel == "angleTilt2" || strLabel == "tilt2")
    {
        return MDL_ANGLETILT2;
    }
    else if (strLabel == "anglePsi2" || strLabel == "psi2")
    {
        return MDL_ANGLEPSI2;
    }
    else if (strLabel == "image")
    {
        return MDL_IMAGE;
    }
    else if (strLabel == "micrograph")
    {
        return MDL_MICROGRAPH;
    }
    else if (strLabel == "CTFModel")
    {
        return MDL_CTFMODEL;
    }
    else if (strLabel == "scale" || strLabel == "Scale")
    {
        return MDL_SCALE;
    }
    else if (strLabel == "intScale")
    {
        return MDL_INTSCALE;
    }
    else if (strLabel == "modelFraction")
    {
        return MDL_MODELFRAC;
    }
    else if (strLabel == "mirrorFraction")
    {
        return MDL_MIRRORFRAC;
    }

    else if (strLabel == "logLikelihood")
    {
        return MDL_LL;
    }
    else if (strLabel == "signalChange")
    {
        return MDL_SIGNALCHANGE;
    }
    else if (strLabel == "wRobust")
    {
        return MDL_WROBUST;
    }
    else if (strLabel == "bgMean")
    {
        return MDL_BGMEAN;
    }
    else if (strLabel == "shiftX" || strLabel == "Xoff")
    {
        return MDL_SHIFTX;
    }
    else if (strLabel == "shiftY" || strLabel == "Yoff")
    {
        return MDL_SHIFTY;
    }
    else if (strLabel == "shiftZ" || strLabel == "Zoff")
    {
        return MDL_SHIFTZ;
    }
    else if (strLabel == "X")
    {
        return MDL_X;
    }
    else if (strLabel == "Y")
    {
        return MDL_Y;
    }
    else if (strLabel == "Z")
    {
        return MDL_Z;
    }
    else if (strLabel == "enabled")
    {
        return MDL_ENABLED;
    }
    else if (strLabel == "originX")
    {
        return MDL_ORIGINX;
    }
    else if (strLabel == "originY")
    {
        return MDL_ORIGINY;
    }
    else if (strLabel == "originZ")
    {
        return MDL_ORIGINZ;
    }
    else if (strLabel == "weight" || strLabel == "Weight")
    {
        return MDL_WEIGHT;
    }
    else if (strLabel == "objId")
    {
        return MDL_OBJID;
    }
    else if (strLabel == "flip" || strLabel == "Flip")
    {
        return MDL_FLIP;
    }
    else if (strLabel == "ref" || strLabel == "Ref")
    {
        return MDL_REF;
    }
    else if (strLabel == "maxCC")
    {
        return MDL_MAXCC;
    }
    else if (strLabel == "cost")
    {
        return MDL_COST;
    }
    else if (strLabel == "serie")
    {
        return MDL_SERIE;
    }
    else if (strLabel == "pMax" || strLabel == "Pmax" || strLabel == "sumP")
    {
        return MDL_PMAX;
    }
    else if (strLabel == "CTFInputParams")
    {
        return MDL_CTFINPUTPARAMS;
    }
    else if (strLabel == "periodogram")
    {
        return MDL_PERIODOGRAM;
    }
    else if (strLabel == "NMADisplacements")
    {
        return MDL_NMA;
    }
    else if (strLabel == "sampling_rate")
    {
        return MDL_SAMPLINGRATE;
    }
    else if (strLabel == "voltage")
    {
        return MDL_VOLTAGE;
    }
    else if (strLabel == "defocusU")
    {
        return MDL_DEFOCUSU;
    }
    else if (strLabel == "defocusV")
    {
        return MDL_DEFOCUSV;
    }
    else if (strLabel == "iterationNumber")
    {
        return MDL_ITER;
    }
    else if (strLabel == "blockNumber")
    {
        return MDL_BLOCK;
    }
    else if (strLabel == "imageMetaData")
    {
        return MDL_IMGMD;
    }
    else if (strLabel == "referenceMetaData")
    {
        return MDL_REFMD;
    }
    else if (strLabel == "sigmaNoise")
    {
        return MDL_SIGMANOISE;
    }
    else if (strLabel == "sigmaOffset")
    {
        return MDL_SIGMAOFFSET;
    }
    else if (strLabel == "sumWeight")
    {
        return MDL_SUMWEIGHT;
    }
    else if (strLabel == "randomSeed")
    {
        return MDL_RANDOMSEED;
    }

    else
    {
        return MDL_UNDEFINED;
    }
}

std::string MetaDataContainer::decodeLabel(MetaDataLabel inputLabel)
{
    switch (inputLabel)
    {
        case MDL_ANGLEROT:
            return std::string("angleRot");
            break;
        case MDL_ANGLETILT:
            return std::string("angleTilt");
            break;
        case MDL_ANGLEPSI:
            return std::string("anglePsi");
            break;
        case MDL_ANGLEROT2:
            return std::string("angleRot2");
            break;
        case MDL_ANGLETILT2:
            return std::string("angleTilt2");
            break;
        case MDL_ANGLEPSI2:
            return std::string("anglePsi2");
            break;
        case MDL_IMAGE:
            return std::string("image");
            break;
        case MDL_COMMENT:
            return std::string("comment");
            break;
        case MDL_MICROGRAPH:
            return std::string("micrograph");
            break;
        case MDL_CTFMODEL:
            return std::string("CTFModel");
            break;
        case MDL_SCALE:
            return std::string("scale");
            break;
        case MDL_INTSCALE:
            return std::string("intScale");
            break;
        case MDL_MODELFRAC:
            return std::string("modelFraction");
            break;
        case MDL_MIRRORFRAC:
            return std::string("mirrorFraction");
            break;
        case MDL_LL:
            return std::string("logLikelihood");
            break;
        case MDL_SIGNALCHANGE:
            return std::string("signalChange");
            break;
        case MDL_WROBUST:
            return std::string("wRobust");
            break;
        case MDL_BGMEAN:
            return std::string("bgMean");
            break;
        case MDL_SHIFTX:
            return std::string("shiftX");
            break;
        case MDL_SHIFTY:
            return std::string("shiftY");
            break;
        case MDL_SHIFTZ:
            return std::string("shiftZ");
            break;
        case MDL_X:
            return std::string("X");
            break;
        case MDL_Y:
            return std::string("Y");
            break;
        case MDL_Z:
            return std::string("Z");
            break;
        case MDL_ENABLED:
            return std::string("enabled");
            break;
        case MDL_ORIGINX:
            return std::string("originX");
            break;
        case MDL_ORIGINY:
            return std::string("originY");
            break;
        case MDL_ORIGINZ:
            return std::string("originZ");
            break;
        case MDL_WEIGHT:
            return std::string("weight");
            break;
        case MDL_OBJID:
            return std::string("objId");
            break;
        case MDL_FLIP:
            return std::string("flip");
            break;
        case MDL_REF:
            return std::string("ref");
            break;
        case MDL_MAXCC:
            return std::string("maxCC");
            break;
        case MDL_COST:
            return std::string("cost");
            break;
        case MDL_SERIE:
            return std::string("serie");
            break;
        case MDL_PMAX:
            return std::string("pMax");
            break;
        case MDL_CTFINPUTPARAMS:
            return std::string("CTFInputParams");
            break;
        case MDL_PERIODOGRAM:
            return std::string("periodogram");
            break;
        case MDL_NMA:
            return std::string("NMADisplacements");
            break;
        case MDL_SAMPLINGRATE:
            return std::string("sampling_rate");
            break;
        case MDL_VOLTAGE:
            return std::string("voltage");
            break;
        case MDL_DEFOCUSU:
            return std::string("defocusU");
            break;
        case MDL_DEFOCUSV:
            return std::string("defocusV");
            break;
        case MDL_ITER:
            return std::string("iterationNumber");
            break;
        case MDL_BLOCK:
            return std::string("blockNumber");
            break;
        case MDL_IMGMD:
            return std::string("imageMetaData");
            break;
        case MDL_REFMD:
            return std::string("referenceMetaData");
            break;
        case MDL_SIGMANOISE:
            return std::string("sigmaNoise");
            break;
        case MDL_SIGMAOFFSET:
            return std::string("sigmaOffset");
            break;
        case MDL_SUMWEIGHT:
            return std::string("sumWeight");
            break;
        case MDL_RANDOMSEED:
            return std::string("randomSeed");
            break;

        default:
            return std::string("");
            break;
    }
}

bool MetaDataContainer::writeValueToString(std::string &outString,
        MetaDataLabel inLabel)
{
    MetaDataLabel inputLabel = static_cast<MetaDataLabel> (inLabel);

    if (valueExists(inLabel))
    {
        std::ostringstream oss;

        if (isDouble(inputLabel))
        {
            double d;
            d = *((double*) (getVoidPtr(inputLabel)));
            if (d != 0. && ABS(d) < 0.001)
                oss << std::setw(12) << std::scientific;
            else
                oss << std::setw(12) << std::fixed;
            oss << d;
        }
        else if (isString(inputLabel))
            oss << *((std::string*) (getVoidPtr(inputLabel)));
        else if (isInt(inputLabel))
        {
            oss << std::setw(10) << std::fixed;
            oss << *((int*) (getVoidPtr(inputLabel)));
        }
        else if (isBool(inputLabel))
            oss << *((bool*) (getVoidPtr(inputLabel)));
        else if (isVector(inputLabel))
        {
            const std::vector<double> &myVector =
                    *((std::vector<double>*) (getVoidPtr(inputLabel)));
            int imax = myVector.size();
            oss << "** ";
            for (int i = 0; i < imax; i++)
            {
                oss << std::setw(12) << std::scientific;
                if (myVector[i] != 0. && ABS(myVector[i]) < 0.001)
                    oss << std::setw(12) << std::scientific;
                else
                    oss << std::setw(12) << std::fixed;
                oss << myVector[i] << " ";
            }
            oss << "**";
        }

        outString = oss.str();
        return true;
    }
    else
    {
        outString = std::string("");
        return false;
    }
}

bool MetaDataContainer::writeValueToFile(std::ofstream &outfile,
        MetaDataLabel inputLabel)
{
    if (valueExists(inputLabel))
    {
        if (isDouble(inputLabel))
        {
            double d;
            d = *((double*) (getVoidPtr(inputLabel)));
            if (d!= 0. && ABS(d) < 0.001)
                outfile << std::setw(12) << std::scientific;
            else
                outfile << std::setw(12) << std::fixed;
            outfile << d;
        }
        else if (isString(inputLabel))
            outfile << *((std::string*) (getVoidPtr(inputLabel)));
        else if (isInt(inputLabel))
        {
            outfile << std::setw(10) << std::fixed;
            outfile << *((int*) (getVoidPtr(inputLabel)));
        }
        else if (isBool(inputLabel))
            outfile << *((bool*) (getVoidPtr(inputLabel)));
        else if (isVector(inputLabel))
        {
            const std::vector<double> &myVector =
                    *((std::vector<double>*) (getVoidPtr(inputLabel)));
            int imax = myVector.size();
            outfile << "** ";
            for (int i = 0; i < imax; i++)
            {
                if (myVector[i] != 0. && ABS(myVector[i]) < 0.001)
                    outfile << std::setw(12) << std::scientific;
                else
                    outfile << std::setw(12) << std::fixed;
                outfile << myVector[i] << " ";
            }
            outfile << "**";
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool MetaDataContainer::isValidLabel(MetaDataLabel inputLabel)
{
    if (inputLabel > MDL_UNDEFINED && inputLabel < MDL_LAST_LABEL)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool MetaDataContainer::isValidLabel(std::string inputLabel)
{
    MetaDataLabel label = codifyLabel(inputLabel);

    if (label > MDL_UNDEFINED && label < MDL_LAST_LABEL)
    {
        return true;
    }
    else
    {
        return false;
    }
}
