/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (1999)
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
#include <XmippData/xmippTypes.hh>

/* PROTOTYPES -------------------------------------------------------------- */
void Usage();

#define PLUS  1
#define MINUS -1

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   string          ang1="+rot",ang2="+tilt",ang3="+psi";
   float           rot, tilt, psi;
   int             sign1=PLUS, sign2=PLUS, sign3=PLUS;
   SelFile         selfile;
   DocFile         angles;
   ImageXmipp      proj;

// Check the command line ==================================================
   try {
      int i;
      // Sel option
      if (check_param(argc,argv,"-sel")) {
         selfile.read(get_param(argc,argv,"-sel"));
         angles.read(get_param(argc,argv,"-ang"));

      // Img option
         // An internal Selection and Document file are generated
      } else if (check_param(argc,argv,"-img")) {
         selfile.insert(get_param(argc,argv,"-img"),SelLine::ACTIVE);
         if ((i=position_param(argc,argv,"-ang"))==-1) {
            cout << "Assign angles: No angular information for single image\n";
            Usage();
         }
         if (i+3>=argc) {
            cout << "Assign angles: Not enough parameters behind -ang\n";
            Usage();
         }
         angles.append_data_line();
         angles.set(1,0,AtoF(argv[i+1]));
         angles.set(1,1,AtoF(argv[i+2]));
         angles.set(1,2,AtoF(argv[i+3]));
	 // second euler angles triplet
         if ((i=position_param(argc,argv,"-ang1"))!=-1) {
            if (i+3>=argc) {
               cout << "Assign angles: Not enough parameters behind -ang1\n";
               Usage();
            }
            angles.set(1,3,AtoF(argv[i+1]));
            angles.set(1,4,AtoF(argv[i+2]));
            angles.set(1,5,AtoF(argv[i+3]));
         }
	 // third euler angles triplet
         if ((i=position_param(argc,argv,"-ang2"))!=-1) {
            if (i+3>=argc) {
               cout << "Assign angles: Not enough parameters behind -ang2\n";
               Usage();
            }
            angles.set(1,6,AtoF(argv[i+1]));
            angles.set(1,7,AtoF(argv[i+2]));
            angles.set(1,8,AtoF(argv[i+3]));
         }
      // Error
      } else {
         cout << "Assign angles: No selection or image file\n";
         Usage();
      }
      if (position_param(argc,argv,"-ang1")==-1 &&
          position_param(argc,argv,"-ang2")!=-1) {
         cout << "-ang2 can not be used without -ang1\n";
         Usage();
      }

      // Angle order
      if ((i=position_param(argc,argv,"-order"))!=-1) {
         if (i+3>=argc) {
            cout << "Assign angles: Not enough parameters behind -ang\n";
            Usage();
         }
         ang1=argv[i+1];
         ang2=argv[i+2];
         ang3=argv[i+3];
      }

   } catch (Xmipp_error XE) {cout << XE; Usage();}

   // Check signs of angles
   if      (ang1[0]=='-') sign1=MINUS;
   else if (ang1[0]=='+') sign1=PLUS;
   else    EXIT_ERROR(1,"Unrecognized sign for first angle");

   if      (ang2[0]=='-') sign2=MINUS;
   else if (ang2[0]=='+') sign2=PLUS;
   else    EXIT_ERROR(1,"Unrecognized sign for second angle");

   if      (ang3[0]=='-') sign3=MINUS;
   else if (ang3[0]=='+') sign3=PLUS;
   else    EXIT_ERROR(1,"Unrecognized sign for third angle");

   // Check they are "rot", "tilt", and "psi"
   check_angle_descr(ang1.substr(1));
   check_angle_descr(ang2.substr(1));
   check_angle_descr(ang3.substr(1));
   if (ang1[1]==ang2[1] || ang1[1]==ang3[1] || ang2[1]==ang3[1])
      EXIT_ERROR(1,"Assign angles: There is an angle twice in the angle order");

// Set angles ==============================================================
   try {
      selfile.go_first_ACTIVE();
      angles.go_first_data_line();
      int FirstLine_ColNo=angles.FirstLine_ColNo();
      while (!selfile.eof()) {
         // Check that there are still angles
         if (angles.eof())
            EXIT_ERROR(1,"Assign angles: There are less angles than images");

         // Choose projection and read
         proj.read(selfile.get_current_file());
         //clear flag that stores number of valid euler triplets
	 proj.clear_fFlag_flag();
         // Assign angles
         switch (ang1[1]) {
            case 'r': rot  = angles(0)*sign1; break;
            case 't': tilt = angles(0)*sign1; break;
            case 'p': psi  = angles(0)*sign1; break;
         }
         switch (ang2[1]) {
            case 'r': rot  = angles(1)*sign2; break;
            case 't': tilt = angles(1)*sign2; break;
            case 'p': psi  = angles(1)*sign2; break;
         }
         switch (ang3[1]) {
            case 'r': rot  = angles(2)*sign3; break;
            case 't': tilt = angles(2)*sign3; break;
            case 'p': psi  = angles(2)*sign3; break;
         }
         proj.set_eulerAngles(rot,tilt,psi);
         // Show
         cout << proj.name() << "  : rot = " << rot << " tilt = " << tilt
              << " psi = " << psi << endl;
         // Second Euler angle triplet
         if(FirstLine_ColNo >=6){
	      switch (ang1[1]) {
        	 case 'r': rot  = angles(3)*sign1; break;
        	 case 't': tilt = angles(3)*sign1; break;
        	 case 'p': psi  = angles(3)*sign1; break;
              }
              switch (ang2[1]) {
        	 case 'r': rot  = angles(4)*sign2; break;
        	 case 't': tilt = angles(4)*sign2; break;
        	 case 'p': psi  = angles(4)*sign2; break;
              }
              switch (ang3[1]) {
        	 case 'r': rot  = angles(5)*sign3; break;
        	 case 't': tilt = angles(5)*sign3; break;
        	 case 'p': psi  = angles(5)*sign3; break;
              }
         proj.set_eulerAngles1(rot,tilt,psi);
         cout << "\t\t" << "  : rot1= " << rot << " tilt1= " << tilt
              << " psi1= " << psi << endl;
	 }     
         if(FirstLine_ColNo >=9){
		 switch (ang1[1]) {
        	 case 'r': rot  = angles(6)*sign1; break;
        	 case 't': tilt = angles(6)*sign1; break;
        	 case 'p': psi  = angles(6)*sign1; break;
              }
              switch (ang2[1]) {
        	 case 'r': rot  = angles(7)*sign2; break;
        	 case 't': tilt = angles(7)*sign2; break;
        	 case 'p': psi  = angles(7)*sign2; break;
              }
              switch (ang3[1]) {
        	 case 'r': rot  = angles(8)*sign3; break;
        	 case 't': tilt = angles(8)*sign3; break;
        	 case 'p': psi  = angles(8)*sign3; break;
              }
         proj.set_eulerAngles2(rot,tilt,psi);
         cout << "\t\t" << "  : rot2= " << rot << " tilt2= " << tilt
              << " psi2= " << psi << endl;
	 }

         // Save
         proj.write();

         // Next file
         selfile.NextImg();
         angles.next_data_line();
      }
   } catch (Xmipp_error XE) {cout << XE;}
}

/* Usage ------------------------------------------------------------------- */
void Usage() {
    printf("Purpose:\n");
    printf(" This program allows you to set up to three Euler angles triplets\n"
           " in the header of an image\n"
           "   (then the angles must be explicitly given in the command line)\n"
           "   or a set of images (in this case the Euler angles come from a\n"
           "   Spider document file\n\n");
    printf("Usage:\n");
    printf("   assign_angles <options>\n");
    printf("   Where <options> are:\n");
    printf("      (-sel <sel_file>             : selection file with the set of images\n"
           "      -ang <doc_file> )|           : Spider document file with the angles\n"
           "      (-img <img_file>             : for setting the angles of a single image\n"
           "      -ang <val1> <val2> <val3>)   : Three float values for ang1, ang2 and ang3\n"
           "                                     respectively.\n"
           "      -ang1 <val1'> <val2'> <val3'>)   : Three float values for ang1', ang2' and ang3'\n"
           "                                     respectively.\n"
           "      -ang2 <val1''> <val2''> <val3''>)   : Three float values for ang1'', ang2'' and ang3''\n"
           "                                     respectively.\n"
           "      [-order <ang1> <ang2> <ang3>]: where ang1, ang2 and ang3 are\n"
           "                                     either +psi, -psi, +tilt, -tilt,\n"
           "                                     +rot or -rot. The default\n"
           "                                     order is +rot, +tilt and +psi.\n"
           "                                     A plus sign in the angle leaves\n"
           "                                     the angle value as it is, while\n"
           "                                     a minus sign changes its sign.\n");
    printf("\n"
           "Examples:\n"
           "   assign_angles -img g0tA0001 -ang 30 40 20\n"
           "   assign_angles -img g0tA0001 -ang 30 40 20 -ang2 23 43 33\n"
           "   assign_angles -img g0tA0001 -ang 30 40 20   -order +rot +tilt +psi\n"
           "   assign_angles -sel g0t.sel  -ang angles.doc\n"
           "   assign_angles -sel g0t.sel  -ang angles.doc -order -rot +tilt +psi\n");
    exit(1);
}

/* Menu -------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Assign_angles {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Assign_angles/Help/assign_angles.html";
      help="Set the angles of an projection on its header";
      OPEN MENU Assign_angles;
      COMMAND LINES {
         + SelFile: assign_angles -sel $SELFILE -ang $ANGFILE
              [-order $ANG1 $ANG2 $ANG3] [$ANGL1 $ANGL2 $ANGL3]
         + Single_file: assign_angles -img $FILE_IN -ang $VANG1 $VANG2 $VANG3
              OPT($ANG1) OPT($ANGL1)
      }
      PARAMETER DEFINITIONS {
        $SELFILE {
            label="Input SelFile";
            type=FILE EXISTING;
        }
        $FILE_IN {
            label="Input Image";
            type=FILE EXISTING;
        }
        $ANGFILE {
            label="Document file with angles";
            type=FILE EXISTING;
        }
        OPT($ANG1) {shown=no;}
           $ANG1 {shown=no;type=text;}
           $ANG2 {COPY $ANG1;}
           $ANG3 {COPY $ANG2;}
        OPT($ANGL1) {label="Angle order";}
           $ANGL1 {
              label="First angle";
              type=list {
                  "+rot"  {$ANG1="+rot";}
                  "+tilt" {$ANG1="+tilt";}
                  "+psi"  {$ANG1="+psi";}
                  "-rot"  {$ANG1="-rot";}
                  "-tilt" {$ANG1="-tilt";}
                  "-psi"  {$ANG1="-psi";}
              };
              by default=0;
           }
           $ANGL2 {
              label="Second angle";
              type=list {
                  "+rot"  {$ANG2="+rot";}
                  "+tilt" {$ANG2="+tilt";}
                  "+psi"  {$ANG2="+psi";}
                  "-rot"  {$ANG2="-rot";}
                  "-tilt" {$ANG2="-tilt";}
                  "-psi"  {$ANG2="-psi";}
              };
              by default=1;
           }
           $ANGL3 {
              label="Third angle";
              type=list {
                  "+rot"  {$ANG3="+rot";}
                  "+tilt" {$ANG3="+tilt";}
                  "+psi"  {$ANG3="+psi";}
                  "-rot"  {$ANG3="-rot";}
                  "-tilt" {$ANG3="-tilt";}
                  "-psi"  {$ANG3="-psi";}
              };
              by default=2;
           }

        $VANG1 {label="Value of angle 1"; type=FLOAT;}
        $VANG2 {label="Value of angle 2"; type=FLOAT;}
        $VANG3 {label="Value of angle 3"; type=FLOAT;}
      }
   }
   MENU Assign_angles {
      "I/O Parameters"
      $SELFILE
      $FILE_IN
      "Angular parameters"
      $ANGFILE
      $VANG1
      $VANG2
      $VANG3
      OPT($ANGL1)
   }
*/
