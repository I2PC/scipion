/***************************************************************************
 *
 * Authors:     Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
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
/*****************************************************************************/
/* INTERACTION WITH PDBs                                                     */
/*****************************************************************************/

#ifndef _XMIPP_PDB_HH
#define _XMIPP_PDB_HH

#include <string>
#include "matrix1d.h"
#include "projection.h"
#include "histogram.h"

/**@defgroup PDBinterface PDB
   @ingroup InterfaceLibrary */
//@{
/** Returns the charge of an atom.
    Returns 0 if the atom is not within the short list (H, C, N, O, S, P, Fe)
    of valid atoms. */
int atomCharge(const std::string &atom);

/** Returns the radius of an atom.
    Returns 0 if the atom is not within the short list (H, C, N, O, S, P, Fe)
    of valid atoms.
    
    The radius data is taken from http://www.webelements.com as the empirical
    radius. */
double atomRadius(const std::string &atom);

/** Compute the center of mass and limits of a PDB file.
    The intensity column is used only for the pseudoatoms. It specifies
    from which column we should read the intensity. Valid columns are
    Bfactor or occupancy.
*/
void computePDBgeometry(const std::string &fnPDB,
                        Matrix1D<double> &centerOfMass,
                        Matrix1D<double> &limit0, Matrix1D<double> &limitF,
                        const std::string &intensityColumn);

/** Apply geometry transformation to an input PDB.
    The result is written in the output PDB. Set centerPDB if you
    want to compute the center of mass first and apply the transformation
    after centering the PDB. */
void applyGeometryToPDBFile(const std::string &fn_in, const std::string &fn_out,
                   const Matrix2D<double> &A, bool centerPDB=true,
                   const std::string &intensityColumn="occupancy");

/** Atom class. */
class Atom
{
public:
    /// Type
    char atomType;

    /// Position X
    double x;

    /// Position Y
    double y;

    /// Position Z
    double z;
};

/** Phantom description using atoms. */
class PDBPhantom
{
public:
    /// List of atoms
    std::vector<Atom> atomList;

    /// Add Atom
    void addAtom(const Atom &atom)
    {
        atomList.push_back(atom);
    }

    /// Get Atom at position i
    const Atom& getAtom(int i) const
    {
        return atomList[i];
    }

    /// Get number of atoms
    size_t getNumberOfAtoms() const
    {
        return atomList.size();
    }

    /// Read from PDB file
    void read(const FileName &fnPDB);

    /// Apply a shift to all atoms
    void shift(double x, double y, double z);

    /** Produce side info.
        The side info produces the radial profiles of each atom
    and its projections.
    */
    void produceSideInfo();
};

/** Atom class. */
class RichAtom
{
public:
    /// Type
    char atomType;

    /// Position X
    double x;

    /// Position Y
    double y;

    /// Position Z
    double z;

    /// Name
    String name;

    /// Alternate location
    char altloc;

    /// Residue name
    String resname;

    /// ChainId
    char chainid;

    /// Residue sequence
    int resseq;

    /// Icode
    char icode;

    /// Occupancy
    double occupancy;

    /// Bfactor
    double bfactor;
};

/** Phantom description using atoms. */
class PDBRichPhantom
{
public:
	/// List of remarks
	std::vector<String> remarks;
public:
    /// List of atoms
    std::vector<RichAtom> atomList;

    /// Add Atom
    void addAtom(const RichAtom &atom)
    {
        atomList.push_back(atom);
    }

    /// Get number of atoms
    size_t getNumberOfAtoms() const
    {
        return atomList.size();
    }

    /// Read from PDB file
    void read(const FileName &fnPDB);

    /// Write to PDB file
    void write(const FileName &fnPDB);

};

/** Description of the electron scattering factors.
    The returned descriptor is descriptor(0)=Z (number of electrons of the
    atom), descriptor(1-5)=a1-5, descriptor(6-10)=b1-5.
    The electron scattering factor at a frequency f (Angstroms^-1)
    is computed as f_el(f)=sum_i(ai exp(-bi*x^2)). Use the function
    electronFormFactorFourier or 
    
    See Peng, Ren, Dudarev, Whelan. Robust parameterization of elastic and
    absorptive electron atomic scattering factors. Acta Cryst. A52: 257-276
    (1996). Table 3 and equation 3.*/
void atomDescriptors(const std::string &atom, Matrix1D<double> &descriptors);

/** Compute the electron Form Factor in Fourier space.
    The electron scattering factor at a frequency f (Angstroms^-1)
    is computed as f_el(f)=sum_i(ai exp(-bi*x^2)). */
double electronFormFactorFourier(double f,
                                 const Matrix1D<double> &descriptors);

/** Compute the electron Form Factor in Real space.
    Konwing the electron form factor in Fourier space is easy to make an
    inverse Fourier transform and express it in real space. r is the
    distance to the center of the atom in Angstroms. */
double electronFormFactorRealSpace(double r,
                                   const Matrix1D<double> &descriptors);

/** Atom radial profile.
    Returns the radial profile of a given atom, i.e., the electron scattering
    factor convolved with a suitable low pass filter for sampling the volume
    at a sampling rate M*T. The radial profile is sampled at T Angstroms/pixel.
*/
void atomRadialProfile(int M, double T, const std::string &atom,
                       Matrix1D<double> &profile);

/** Atom projection radial profile.
    Returns the radial profile of the atom described by its profileCoefficients
    (Bspline coefficients). */
void atomProjectionRadialProfile(int M,
                                 const Matrix1D<double> &profileCoefficients,
                                 Matrix1D<double> &projectionProfile);

/** Class for Atom interpolations. */
class AtomInterpolator
{
public:
    // Vector of radial volume profiles
    std::vector< MultidimArray<double> > volumeProfileCoefficients;
    // Vector of radial projection profiles
    std::vector< MultidimArray<double> > projectionProfileCoefficients;
    // Vector of atom radii
    std::vector<double> radii;
    // Downsampling factor
    int M;
    // Fine sampling rate
    double highTs;

    /** Setup.
        HighTs is a fine sampling rate, M is an integer number so that
    the final sampling rate is Ts=M*highTs; */
    void setup(int m, double hights, bool computeProjection=false);

    /// Add atom
    void addAtom(const std::string &atomType, bool computeProjection=false);

    /// Get atom index
    int getAtomIndex(char atom) const
    {
        int idx=-1;
        switch (atom)
        {
        case 'H':
            idx=0;
            break;
        case 'C':
            idx=1;
            break;
        case 'N':
            idx=2;
            break;
        case 'O':
            idx=3;
            break;
        case 'P':
            idx=4;
            break;
        case 'S':
            idx=5;
            break;
        case 'F':
            idx=6;
            break;
        default:
            REPORT_ERROR(ERR_VALUE_INCORRECT,(std::string)
                         "AtomInterpolator::getAtomIndex: Atom "+atom+" unknown");
        }
        return idx;
    }

    /** Radius of an atom in the final sampling rate M*highTs. */
    double atomRadius(char atom) const
    {
        return radii[getAtomIndex(atom)];
    }

    /** Volume value at a distance r of the atom whose first letter
        is the one provided as atom. */
    double volumeAtDistance(char atom, double r) const
    {
        int idx=getAtomIndex(atom);
        if (r>radii[idx])
            return 0;
        else
            return volumeProfileCoefficients[idx].
                   interpolatedElementBSpline1D(r*M,3);
    }

    /** Projection value at a distance r of the atom whose first letter
        is the one provided as atom. */
    double projectionAtDistance(char atom, double r) const
    {
        int idx=getAtomIndex(atom);
        if (r>radii[idx])
            return 0;
        else
            return projectionProfileCoefficients[idx].
                   interpolatedElementBSpline1D(r*M,3);
    }
};

/** Project PDB.
    Project the PDB following a certain projection direction. */
void projectPDB(const PDBPhantom &phantomPDB,
                const AtomInterpolator &interpolator, Projection &proj,
                int Ydim, int Xdim, double rot, double tilt, double psi);

/** Compute distance histogram of a PDB phantom.
 * Consider the distance between each atom and its N nearest neighbours. Then, compute the histogram of these distances
 * with Nbin samples.
 */
void distanceHistogramPDB(const PDBPhantom &phantomPDB, size_t Nnearest, double maxDistance, int Nbins, Histogram1D &hist);
//@}
#endif
