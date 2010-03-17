#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <data/args.h>
#include <data/progs.h>
#include <data/docfile.h>
#include <data/macros.h>
#include <interface/aph.h>

#define PIRAD 3.1416

void Usage(char *argv[]);
void EM2Tilt(double * , double * , double *, double *, double *);
void MRC2Euler(double *, double *, double *);

int main(int argc, char **argv)
{
    /////////////////////// Command Line Parameters
    // I/O Files
    std::string filein;
    std::string fileout;
    // Angular System
    std::string AngSystem;
    // Space
    std::string Space;

    /////////////////////// IO Files
    DocFile EMSeriesFile;
    DocFile AnglesFile1;
    DocFile AnglesFile2;
    DocLine Angles;
    /////////////////////// EM Lat. Variables
    double Ref[4];
    double Tilt[4];
    double DetRef, DetTilt;
    /////////////////////// Angle Variables
    double Ang1[3]; //Stores either Euler ang. or Taxa-Tilt
    double Ang2[3]; //Stores either Euler ang. or Taxa-Tilt
    double aux;
    /////////////////////// Scaling
    double Scale[1];


    //////////////////////// Flags
    bool Euler;
    bool Fourier;

    // Get command line parameters ------------------------------------------
    try
    {
        filein  = getParameter(argc, argv, "-i");
        fileout = getParameter(argc, argv, "-o");
        AngSystem = getParameter(argc, argv, "-ang_syst", "Euler");
        Space = getParameter(argc, argv, "-space", "Real");

    }
    catch (Xmipp_error XE)
    {
        Usage(argv);
    }


    Euler = (bool)(strcmp(AngSystem.c_str(), "Euler") == 0);
    Fourier = (bool)(strcmp(Space.c_str(), "Fourier") == 0);


    //Read EM Lat.Parameters
    EMSeriesFile.read(filein.c_str());
    EMSeriesFile.go_first_data_line();

    AnglesFile1.clear();
    AnglesFile2.clear();

    if (Euler)
    {
        AnglesFile1.append_comment((std::string) "LinNo  NoElements  Im_ID  rot tilt psi Scaling");
        AnglesFile2.append_comment((std::string) "LinNo  NoElements  Im_ID  rot tilt psi Scaling");
    }
    else
    {
        AnglesFile1.append_comment((std::string) "LinNo  NoElements  Im_ID  taxa3D taxa2D tilt Scaling");
        AnglesFile2.append_comment((std::string) "LinNo  NoElements  Im_ID  taxa3D taxa2D tilt Scaling");
    }

    Angles.set_type(DocLine::DATALINE);


    //A Vector on Reference Projection
    Ref[0] = EMSeriesFile(1);
    Ref[1] = EMSeriesFile(2);
    //B Vector on Reference Projection
    Ref[2] = EMSeriesFile(3);
    Ref[3] = EMSeriesFile(4);
    //Orientation from B to A
    DetRef = -(Ref[0] * Ref[3] - Ref[1] * Ref[2]);



    EMSeriesFile.next_data_line();

    //Compute Angles
    while (!EMSeriesFile.eof())
    {

        //A Vector on Tilted Projection
        Tilt[0] = EMSeriesFile(1);
        Tilt[1] = EMSeriesFile(2);
        //B Vector on Tilted Projection
        Tilt[2] = EMSeriesFile(3);
        Tilt[3] = EMSeriesFile(4);
        //Orientation from B to A
        DetTilt = -(Tilt[0] * Tilt[3] - Tilt[1] * Tilt[2]);



        //Taxa-TiltAngle
        if (Fourier)
            EM2Tilt(Ref, Tilt, Scale, Ang1, Ang2);
        else
        {
            //Swap Tilt-Ref for Real Space
            EM2Tilt(Tilt, Ref, Scale, Ang1, Ang2);
            aux = Ang1[0];
            Ang1[0] = Ang1[1];
            Ang1[1] = aux;
            aux = Ang2[0];
            Ang2[0] = Ang2[1];
            Ang2[1] = aux;
            Scale[0] = 1 / Scale[0];
        }



        //Euler if required
        if (Euler)
        {
            //Ensure orientation is consistent
            Ang1[0] = SGN(DetRef) * Ang1[0];
            Ang1[1] = SGN(DetTilt) * Ang1[1];
            MRC2Euler(Ref, Tilt, Ang1);
            //Ensure orientation is consistent
            Ang2[0] = SGN(DetRef) * Ang2[0];
            Ang2[1] = SGN(DetTilt) * Ang2[1];
            MRC2Euler(Ref, Tilt, Ang2);
        }




        //Write Angles to file
        Angles.set(0, EMSeriesFile(0));
        Angles.set(1, Ang1[0]);
        Angles.set(2, Ang1[1]);
        Angles.set(3, Ang1[2]);
        Angles.set(4, Scale[0]);

        AnglesFile1.append_line(Angles);


        //Write Angles to file
        Angles.set(0, EMSeriesFile(0));
        Angles.set(1, Ang2[0]);
        Angles.set(2, Ang2[1]);
        Angles.set(3, Ang2[2]);
        Angles.set(4, Scale[0]);

        AnglesFile2.append_line(Angles);



        //Next Projection
        EMSeriesFile.next_data_line();

    }

    std::string filename = fileout + "1";
    AnglesFile1.write(filename.c_str());
    filename = fileout + "2";
    AnglesFile2.write(filename.c_str());


    return 0;

}

void Usage(char *argv[])
{
    std::cout << "Purpose:\n"
    << "    Computes Angles from EMTilt series\n"
    << "Usage:" << argv[0] << " -i filein -o fileout -ang_system AngSystem" <<  std::endl
    << "\t-i               :  Input EMLatVectors file" << std::endl
    << "\t-o               :  Output EulerAngles file" << std::endl
    << "\t-ang_syst        :  Angular System (Euler or MRCTilt)" << std::endl
    << "\t-space           :  Space (Real or Fourier)" << std::endl
    ;
    exit(1);

}

//////////////// Computation of Projection Parameters ////////////////////////
//
//
// EM2Tilt (based on MRC emtilt)
// Computes Scaling, Taxa & Tilt Angle (in degres) from Lat.Vectors
// on two diferent projection planes (Ref,Tilt).

//  CALCULATE TILT ANGLES FROM TILT AND REF
//  CELL DIMENSIONS.
//
//  THE ALGORITHM IS DESCRIBED IN
//              SHAW AND HILLS. MICROS., 12:279-283 (1981).
//
//  CONVENTION FOR MEASURING TILT AXIS TO A IS THAT THE ANGLE IS
//  FROM TILTAXIS TO A IN THE DIRECTION GIVEN BY
//   1.  A TO B FOR MRC TAXA AND TILT. THIS IS NOT THE CONVENTION USED IN THE PAPER BY SHAW.
//   2.  B TO A FOR EULER ANGLES COMPUTATION. THIS IS THE CONVENTION USED IN THE PAPER BY SHAW.
//  BEING POSITIVE.

void EM2Tilt(double * Ref, double * Tilt, double *Scale, double *Ang1, double * Ang2)
{
    //Lat. Dimensions
    double A, B, GAMMA, AT, BT, GAMMAT;
    A = sqrt(Ref[0] * Ref[0] + Ref[1] * Ref[1]);
    B = sqrt(Ref[2] * Ref[2] + Ref[3] * Ref[3]);
    GAMMA = acos((Ref[0] * Ref[2] + Ref[1] * Ref[3]) / (A * B));

    AT = sqrt(Tilt[0] * Tilt[0] + Tilt[1] * Tilt[1]);
    BT = sqrt(Tilt[2] * Tilt[2] + Tilt[3] * Tilt[3]);
    GAMMAT = acos((Tilt[0] * Tilt[2] + Tilt[1] * Tilt[3]) / (AT * BT));

    double COSG = cos(GAMMA);
    double COSGT = cos(GAMMAT);
    double SinG = sin(GAMMA);

    double C1 = (A * A) / (AT * AT);
    double C2 = (B * B) / (BT * BT);
    double C3 = (A * B) / (AT * BT * COSGT);
    double C4 = (A * B * COSG) / (AT * BT * COSGT);

    //Scaling
    double AX = C2 * (C1 / C3) * (C1 / C3) - C1;
    double BX = C2 + 2.0 * C2 * (C1 / C3) * ((C1 - C4) / C3) - C1;
    double CX = C2 * ((C1 - C4) / C3) * ((C1 - C4) / C3);
    double DISC = BX * BX - 4.0 * AX * CX;
    if (DISC < 0.0)
    {
        std::cout << "Discriminant less than zero" << std::endl;
        return;
    }
    double PSQ1 = (-BX + sqrt(DISC)) / (2.0 * AX);
    double  PSQ2 = (-BX - sqrt(DISC)) / (2.0 * AX);
    double PP1, PP2;

    if (PSQ1 > 0.0)
        PP1 = sqrt(PSQ1);
    else if (PSQ2 > 0.0)
        PP1 = sqrt(PSQ2);
    else
    {
        std::cout << "Two Negative Roots.Something Wrong" << std::endl;
        return;
    }

    PP2 = -PP1;
    double QQ1 = (C1 / C3) * PP1 + ((C1 - C4) / C3) * (1.0 / PP1);
    double QQ2 = -QQ1;

    Scale[0] = sqrt(C1 * (PP1 * PP1 + 1.0));



    //Tilt Axes and Tilt Angle
    // Tilt Axe in 3D

    Ang1[0] = atan2(SinG, QQ1 / PP1 - COSG);
    Ang2[0] = atan2(SinG, QQ2 / PP2 - COSG);

    double SINPHI1 = sin(Ang1[0]);
    double COSPHI1 = cos(Ang1[0]);
    double SINPHI2 = sin(Ang2[0]);
    double COSPHI2 = cos(Ang2[0]);
    // Tilt Angle
    Ang1[2] = atan2(PP1, SINPHI1);
    Ang2[2] = atan2(PP2, SINPHI2);


    // Tilt Axe in 2D projection plane
    // THIS IS THE ANGLE ONE WOULD GET DIRECTLY FROM THE FILM,
    // FOR EXAMPLE BY FINDING THE DIRECTION OF ZERO CHANGE IN THE C.T.F.
    //double TANTHE1 = tan(Ang1[2]);
    //COSPHI1 = COSPHI1/(sqrt(1.0 + SINPHI1*SINPHI1*TANTHE1*TANTHE1));
    //Ang1[1] = acos(COSPHI1);
    //double TANTHE2 = tan(Ang2[2]);
    //COSPHI2 = COSPHI2/(sqrt(1.0 + SINPHI2*SINPHI2*TANTHE2*TANTHE2));
    //Ang2[1] = acos(COSPHI2);

    // Computation using atan. By using acos you skip one possible angle
    double COSTHE1 = cos(Ang1[2]);
    double SINPSI1 = SINPHI1 / COSTHE1;
    Ang1[1] = atan2(SINPSI1, COSPHI1);
    double COSTHE2 = cos(Ang2[2]);
    double SINPSI2 = SINPHI2 / COSTHE2;
    Ang2[1] = atan2(SINPSI2, COSPHI2);

    // Convert to degrees
    Ang1[0] = Ang1[0] * 180 / PIRAD;
    Ang1[1] = Ang1[1] * 180 / PIRAD;
    Ang1[2] = Ang1[2] * 180 / PIRAD;
    Ang2[0] = Ang2[0] * 180 / PIRAD;
    Ang2[1] = Ang2[1] * 180 / PIRAD;
    Ang2[2] = Ang2[2] * 180 / PIRAD;

    return;
}

// MRC2Euler.
// Computes Euler Angles (rot,tilt,psi) from Taxa3D, Taxa2D & Tilt Angle (Ang)
// and Lattice Vetors on Ref (Tilt=0) and Tilt projection planes
// Uses vector Ang to store new angular values
//


void MRC2Euler(double *Ref, double *Tilt, double * Ang)
{
    //Vector a in 3D
    double a3D[2];
    a3D[0] = Ref[0];
    a3D[1] = Ref[1];
    //Projected Vector a on Tilt plane
    double a[2];
    a[0] = Tilt[0];
    a[1] = Tilt[1];
    //Vector angles
    double Ang_a3D, Ang_a;
    double rot, tilt, psi;



    //Angle Tilt.
    tilt = Ang[2];


    //Angle psi.
    //Angle between Tilt Axis (given by Taxa2D) and y-axis on Tilt plane
    //Taxa2D
    Ang_a = atan2(a[1], a[0]) * 180 / PIRAD;
    psi = Ang_a + Ang[1];
    psi = 90 - psi;

    //Angle rot.
    //Angle between y-axis Tilt Axis (given by Taxa3D) on 3D space (on Reference Plane)
    Ang_a3D = atan2(a3D[1], a3D[0]) * 180 / PIRAD;
    rot = Ang_a3D + Ang[0];
    rot = rot - 90;


    //Store new angular values
    Ang[0] = rot;
    Ang[1] = tilt;
    Ang[2] = psi;

    return;

}
