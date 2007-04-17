/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/


/* INCLUDES ---------------------------------------------------------------- */
#include <Reconstruction/Programs/Prog_assign_CTF.hh>
#include <XmippData/xmippArgs.hh>

/****************************************************************************

   NAME:        Usage function
   
   DESCRIPTION:	Show program usage 
				
   INPUT:        
				
   OUTPUT:      Values are shown in standar output 
   
****************************************************************************/			 
void Usage()
{
   cerr << "This program generates and assigns a CTF to a selfile of particles.\n"
		<< "   -i <parameters file>         : Name of the file of parameters with all info for the program.\n"
   ;
}

/****************************************************************************

   NAME:        main
   
   DESCRIPTION:	This program calculates and assings a CTF for every particle
   				picked in a micrograph
			
   INPUT:       
   				
   OUTPUT:      
   
****************************************************************************/			 
int main (int argc,char *argv[])
{
   Prog_assign_CTF_prm     prm;

/*****************************************************************************
		 Check the command line and get the parameters
/*****************************************************************************/
   try
   {
      // Read the file with all the program parameters  
      prm.read(get_param(argc,argv,"-i"));
   }
   catch (Xmipp_error &XE) {cout << XE; Usage(); exit(1);}

   try
   {
      prm.process();
   }
   catch (Xmipp_error &XE) {cout << XE; exit(1);}
}

/* Menus ------------------------------------------------------------------- */
