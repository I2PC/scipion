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

#include "strings.h"

std::string removeChar( std::string str, char character )
{
	std::string temp;
			
	for( unsigned int i = 0 ; i < str.length( ) ; i++ )
	{
		if ( str[ i ] != character ) 
			temp += str[ i ];
	}
	
	return temp;
}

std::string unescape( std::string str )
{
	std::string temp;

	for( unsigned int i = 0 ; i < str.length( ) ; i++ )
	{
		char current_char = str[ i ];
			
		if( current_char != '\n' && current_char != '\t' && 
			current_char != '\v' && current_char != '\b' &&
			current_char != '\r' && current_char != '\f' &&
			current_char != '\a' )
		{
			temp += str[ i ];
		}	
	}
	
	return temp;
}

std::string simplify( std::string str )
{
	std::string temp;
	
	// First, unescape string
	str = unescape( str );
	
	// Remove spaces from the beginning
	int pos = str.find_first_not_of( ' ' );
	str.erase( 0, pos );
	
	// Trim the rest of spaces
	for( unsigned int i = 0 ; i < str.length( ) ; )
	{
		temp += str[ i ];
		
		if ( str[ i ] == ' ' )
		{			
			while( str[ i ] == ' ' )
			{
				i++;
			}
		}
		else
		{
			i++;
		}
	}
	
	// Remove space left at the end of the string 
	// if needed
	if( temp[ temp.size( ) - 1 ] == ' ' )
	{
		temp.resize( temp.size() - 1 );
	}
	
	return temp;
}

