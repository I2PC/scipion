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

MetaDataContainer::MetaDataContainer(){};
MetaDataContainer::~MetaDataContainer(){};
MetaDataContainer& MetaDataContainer::operator = ( MetaDataContainer &MDc)
{

    if (this != &MDc)
    {
    	void * aux;
    	MetaDataLabel lCode;
    	std::map<MetaDataLabel, void *>::iterator It;
		for( It = (MDc.values).begin( ) ; It != (MDc.values).end( ); It ++ )
		{
			aux   = It->second;
		    lCode = It->first;

			if( IS_DOUBLE(lCode) )
			{
				addValue( lCode, * ((double *) aux) );
			}
			else if( IS_STRING(lCode) )
			{
				addValue( lCode, * ((std::string *) aux) );
			}
			else if( IS_INT(lCode))
			{
				addValue( lCode, * ((int *) aux) );
			}
			else if( IS_BOOL(lCode))
			{
				addValue( lCode, * ((bool *) aux) );
			}
			else if( IS_VECTOR(lCode) )
			{
				addValue( lCode, * ((std::vector<double> *) aux) );
			}

		}
    }
    return *this;
}
MetaDataContainer::MetaDataContainer  ( MetaDataContainer &MDc)
{
	void * aux;
	MetaDataLabel lCode;
	std::map<MetaDataLabel, void *>::iterator It;
	for( It = (MDc.values).begin( ) ; It != (MDc.values).end( ); It ++ )
	{
		aux   = It->second;
		lCode = It->first;

		if( IS_DOUBLE(lCode) )
		{
			addValue( lCode, * ((double *) aux) );
		}
		else if( IS_STRING(lCode) )
		{
			addValue( lCode, * ((std::string *) aux) );
		}
		else if( IS_INT(lCode))
		{
			addValue( lCode, * ((int *) aux) );
		}
		else if( IS_BOOL(lCode))
		{
			addValue( lCode, * ((bool *) aux) );
		}
		else if( IS_VECTOR(lCode) )
		{
			addValue( lCode, * ((std::vector<double> *) aux) );
		}

	}
}

void MetaDataContainer::addValue( MetaDataLabel name, int value )
{
	void * newValue = (void *)(new int(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name, double value )
{
    void * newValue = (void *)(new double(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name, bool value )
{
	void * newValue = (void *)(new bool(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name, const std::string &value )
{
	void * newValue = (void *)(new std::string(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( MetaDataLabel name,
    const std::vector<double> &value )
{
	void * newValue = (void *)(new std::vector<double>(value));
	insertVoidPtr( name, newValue );
}

void MetaDataContainer::addValue( const std::string &name,
    const std::string &value )
{	
	MetaDataLabel lCode = codifyLabel( name );
	std::istringstream i( value );
	
	// Look for a double value
	if( IS_DOUBLE(lCode) )
	{
		double doubleValue;
		
		i >> doubleValue;			
				
		addValue( lCode, doubleValue );
	}
	else if( IS_STRING(lCode) )
	{
		addValue( lCode, value );
	}
	else if( IS_INT(lCode))
	{
		int intValue;
		
		i >> intValue;
				
		addValue( lCode, intValue ); 
	}
	else if( IS_BOOL(lCode))
	{
		bool boolValue;
		
		i >> boolValue;
				
		addValue( lCode, boolValue );
	}
	else if( IS_VECTOR(lCode))
	{
		std::vector<double> vectorValue;
                double val;
                while (i >> val)
                    vectorValue.push_back(val);
		addValue( lCode, vectorValue );
	}
}

void MetaDataContainer::insertVoidPtr( MetaDataLabel name, void * value )
{
	values[ name ] = value;
}

void  MetaDataContainer::getValue( MetaDataLabel name, int &value )
{	
	std::map<MetaDataLabel, void *>::iterator element; 

	element = values.find( name );

	if ( element == values.end( ) )
	{
		REPORT_ERROR(1,(std::string) "Label(int) " + decodeLabel(name) + " not found\n" );
	}
	else
	{
		value = *((int *)element->second);
	}
}
void  MetaDataContainer::getValue( MetaDataLabel name, double &value )
{
	std::map<MetaDataLabel, void *>::iterator element;

	element = values.find( name );

	if ( element == values.end( ) )
	{
		REPORT_ERROR(1,(std::string) "Label(double) " + decodeLabel(name) + " not found\n" );
	}
	else
	{
		value = *((double *)element->second);
	}
}

void  MetaDataContainer::getValue( MetaDataLabel name, std::string &value )
{
	std::map<MetaDataLabel, void *>::iterator element;

	element = values.find( name );

	if ( element == values.end( ) )
	{
		REPORT_ERROR(1,(std::string) "Label(string) " + decodeLabel(name) + " not found\n" );
	}
	else
	{
		value = *((std::string *)element->second);
	}
}

void  MetaDataContainer::getValue( MetaDataLabel name, bool &value )
{
	std::map<MetaDataLabel, void *>::iterator element;

	element = values.find( name );

	if ( element == values.end( ) )
	{
		REPORT_ERROR(1,(std::string) "Label(bool) " + decodeLabel(name) + " not found\n" );
	}
	else
	{
		value = *((bool *)element->second);
	}
}

void  MetaDataContainer::getValue( MetaDataLabel name, std::vector<double> &value )
{
	std::map<MetaDataLabel, void *>::iterator element;

	element = values.find( name );

	if ( element == values.end( ) )
	{
		REPORT_ERROR(1,(std::string) "Label(vector) " + decodeLabel(name) + " not found\n" );
	}
	else
	{
		value = *((std::vector<double> *)element->second);
	}
}

void * MetaDataContainer::getVoidPtr( MetaDataLabel name )
{
	std::map<MetaDataLabel, void *>::iterator element;

	element = values.find( name );

	if ( element == values.end( ) )
	{
		REPORT_ERROR(1,(std::string) "Label(bool) " + decodeLabel(name) + " not found\n" );
	}
	else
	{
		return element->second;
	}
}

bool MetaDataContainer::valueExists( MetaDataLabel name )
{
	if( values.find( name ) == values.end( ) )
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool MetaDataContainer::pairExists( MetaDataLabel name, double value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
    double * currValue = (double *)values[ name ];
	
    if( currValue != NULL )
	{
        double val=value- *currValue;
        if( ABS(val)<XMIPP_EQUAL_ACCURACY )
		{
			return true;
		}
	}
	
    return false;
}

bool MetaDataContainer::pairExists( MetaDataLabel name, int value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
    int * currValue = (int *)values[ name ];
	
    if( currValue != NULL )
	{        
		if( value == *currValue )
		{
			return true;
		}
	}
	
	return false;
}

bool MetaDataContainer::pairExists( MetaDataLabel name, bool value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
    bool * currValue = (bool *)values[ name];
    
	if( currValue != 0 )
	{
	    if( value == *currValue )
		{
			return true;
		}
	}

	return false;
}

bool MetaDataContainer::pairExists( MetaDataLabel name,
    const std::string &value )
{
	// Traverse all the structure looking for objects
	// that satisfy search criteria
    std::string * currValue = (std::string *)values[ name ];	
	
    if( currValue != 0 )
	{
	    if( value == *currValue )
		{
			return true;
		}
	}

	return false;
}

void MetaDataContainer::deleteValue( MetaDataLabel name )
{
	//I do not fully understan this but free can no be correct, may delete is OK
	//ROB?
    if( values[ name ] != NULL )
        delete( values[ name ] );
	
    values.erase( name );
}
	
MetaDataLabel MetaDataContainer::codifyLabel( std::string strLabel )
{
	// NOTE: For the "||" cases, the left MetaDataLabel is the XMIPP_3 version's and
	// the right value is the old-styled one

	if( strLabel == "angleRot" || strLabel == "rot" )
	{
		return MDL_ANGLEROT;
	}
	else if( strLabel == "angleTilt" || strLabel == "tilt" )
	{
		return MDL_ANGLETILT;
	}
	else if( strLabel == "anglePsi" || strLabel == "psi" )
	{
		return MDL_ANGLEPSI;
	}
	else if( strLabel == "image" )
	{
		return MDL_IMAGE;
	}
	else if( strLabel == "micrograph" )
	{
		return MDL_MICROGRAPH;
	}
	else if( strLabel == "CTFModel" )
	{
		return MDL_CTFMODEL;
	}
	else if( strLabel == "shiftX" || strLabel == "Xoff" )
	{
		return MDL_SHIFTX;
	}
	else if( strLabel == "shiftY" || strLabel == "Yoff" )
	{
		return MDL_SHIFTY;
	}
	else if( strLabel == "shiftZ" || strLabel == "Zoff" )
	{
		return MDL_SHIFTZ;
	}
	else if( strLabel == "enabled" )
	{
		return MDL_ENABLED;
	}
	else if( strLabel == "originX" )
	{
		return MDL_ORIGINX;
	}
	else if( strLabel == "originY" )
	{
		return MDL_ORIGINY;
	}
	else if( strLabel == "originZ" )
	{
		return MDL_ORIGINZ;
	}
	else if( strLabel == "weight" || strLabel == "Weight"  )
	{
		return MDL_WEIGHT;
	}
	else if( strLabel == "flip" || strLabel == "Flip" )
	{
		return MDL_FLIP;
	}
	else if( strLabel == "ref" || strLabel == "Ref" )
	{
		return MDL_REF;
	}
	else if( strLabel == "maxCC" )
	{
		return MDL_MAXCC;
	}
	else if( strLabel == "serie" )
	{
		return MDL_SERIE;
	}
	else if( strLabel == "pMax" || strLabel == "Pmax" || strLabel == "sumP" )
	{
		return MDL_PMAX;
	}
	else if( strLabel == "CTFInputParams" )
	{
		return MDL_CTFINPUTPARAMS;
	}
	else if( strLabel == "periodogram" )
	{
		return MDL_PERIODOGRAM;
	}
	else if( strLabel == "NMADisplacements" )
	{
		return MDL_NMA;
	}
	else
	{
		return MDL_UNDEFINED;
	}
}

std::string MetaDataContainer::decodeLabel( MetaDataLabel inputLabel )
{
	switch ( inputLabel ) 
	{
		case MDL_ANGLEROT:
			return std::string( "angleRot" );
			break;
		case MDL_ANGLETILT:
			return std::string( "angleTilt" );
			break;
		case MDL_ANGLEPSI:
			return std::string( "anglePsi" );
			break;
		case MDL_IMAGE:
			return std::string( "image" );
			break;
		case MDL_MICROGRAPH:
			return std::string( "micrograph" );
			break;
		case MDL_CTFMODEL:
			return std::string( "CTFModel" );
			break;
		case MDL_SHIFTX:
			return std::string( "shiftX" );
			break;
		case MDL_SHIFTY:
			return std::string( "shiftY" );
			break;
		case MDL_SHIFTZ:
			return std::string( "shiftZ" );
			break;
		case MDL_ENABLED:
			return std::string( "enabled" );
			break;
		case MDL_ORIGINX:
			return std::string( "originX" );
			break;
		case MDL_ORIGINY:
			return std::string( "originY" );
			break;
		case MDL_ORIGINZ:
			return std::string( "originZ" );
			break;
		case MDL_WEIGHT:
			return std::string( "weight" );
			break;
		case MDL_FLIP: 
			return std::string( "flip" );
			break;
		case MDL_REF: 
			return std::string( "ref" );
			break;
		case MDL_MAXCC: 
			return std::string( "maxCC" );
			break;
        case MDL_SERIE:
            return std::string( "serie" );
            break;
        case MDL_PMAX:
            return std::string( "pMax" );
            break;
        case MDL_CTFINPUTPARAMS:
            return std::string( "CTFInputParams" );
            break;
        case MDL_PERIODOGRAM:
            return std::string( "periodogram" );
            break;
        case MDL_NMA:
            return std::string( "NMADisplacements" );
            break;
		default:
			return std::string( "" );
			break;
	}
}

void MetaDataContainer::writeValueToString( std::string &outString,
    MetaDataLabel inputLabel )
{	
    std::ostringstream oss;

    if( IS_DOUBLE(inputLabel) )
    {
        oss << std::setprecision(10) << std::setw(17) << std::scientific;
	    oss << *((double*)(getVoidPtr( inputLabel )));
    }
    else if( IS_STRING(inputLabel) )
	oss << *((std::string*)(getVoidPtr( inputLabel )));
    else if( IS_INT(inputLabel))
    {
        oss << std::setprecision(10) << std::setw(12) << std::fixed;
	    oss << *((int*)(getVoidPtr( inputLabel )));
    }
    else if( IS_BOOL(inputLabel))
	    oss << *((bool*)(getVoidPtr( inputLabel )));
    else if( IS_VECTOR(inputLabel))
    {
        const std::vector<double> &myVector=*(
            (std::vector<double>*)(getVoidPtr( inputLabel )));
        int imax=myVector.size();
        oss << "** ";
        for (int i=0; i<imax; i++)
        {
            oss << std::setprecision(10) << std::setw(17) << std::scientific;
            oss << myVector[i] << " ";
        }
	    oss << "**";
    }
    
    outString = oss.str( );

}

void MetaDataContainer::writeValueToFile( std::ofstream &outfile,
    MetaDataLabel inputLabel )
{
    if( IS_DOUBLE(inputLabel) )
    {
        outfile << std::setprecision(10) << std::setw(17) << std::scientific;
    	outfile << *((double*)(getVoidPtr( inputLabel )));
    }
    else if( IS_STRING(inputLabel) )
    	outfile << *((std::string*)(getVoidPtr( inputLabel )));
    else if( IS_INT(inputLabel))
    {
        outfile << std::setprecision(10) << std::setw(12) << std::fixed;
	    outfile << *((int*)(getVoidPtr( inputLabel )));
    }
    else if( IS_BOOL(inputLabel))
    	outfile << *((bool*)(getVoidPtr( inputLabel )));
    else if( IS_VECTOR(inputLabel))
    {
        const std::vector<double> &myVector=*(
            (std::vector<double>*)(getVoidPtr( inputLabel )));
        int imax=myVector.size();
        outfile << "** ";
        for (int i=0; i<imax; i++)
        {
            outfile << std::setprecision(10) << std::setw(17) << std::scientific;
            outfile << myVector[i] << " ";
        }
    	outfile << "**";
    }
}

bool MetaDataContainer::isValidLabel( MetaDataLabel inputLabel )
{
	if( inputLabel > MDL_UNDEFINED && inputLabel < MDL_LAST_LABEL )
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool MetaDataContainer::isValidLabel( std::string inputLabel )
{
	MetaDataLabel label = codifyLabel( inputLabel );
	
	if( label > MDL_UNDEFINED && label < MDL_LAST_LABEL )
	{
		return true;
	}
	else
	{
		return false;
	}
}
