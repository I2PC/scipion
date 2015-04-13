#!/usr/bin/env python
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# m_inout.py  - Library for coordinate input and output
# Part of the ModeHunter package, http://modehunter.biomachina.org
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Request to cite published literature:
#
# When using this software in scholarly work please cite the publications(s)
# posted at http://modehunter.biomachina.org/fref.html
# 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Legal notice:
# 
# This software is copyrighted, (c) 2009, by Joseph N. Stember and Willy Wriggers
# under the following terms:
#
# The authors hereby grant permission to use, copy, modify, and re-distribute this
# software and its documentation for any purpose, provided that existing copyright
# notices are retained in all copies and that this notice is included verbatim in
# any distributions. No written agreement, license, or royalty fee is required for
# any of the authorized uses.
#
# IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR DIRECT,
# INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE
# OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF, EVEN IF THE
# AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
# BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN "AS
# IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import numpy 
import sys
import copy

# global list with chemical element names
element_name= [ "H", "HE", "LI", "BE","B","C","N","O","F", "NE", \
                "NA", "MG", "AL", "SI","P","S", "CL", "AR","K", "CA", \
                "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN" ]

# global list with corresponding atom masses
element_mass= [ 1.008 ,4.003, 6.940, 9.012, 10.810, 12.011, 14.007, 16.000, 18.998, 20.179,\
         22.990, 24.305, 26.982, 28.086, 30.974, 32.060, 35.453, 39.948, 39.102, 40.080, \
         44.956, 47.880, 50.942, 51.996, 54.938, 55.847, 58.933, 58.700, 63.546, 65.380 ]

# global PDB class with data required structure
class PDB(object):
    def __init__(self, recd="ATOM  ", serial=1, type="  ", loc="  ", alt=" ", res="    ", chain=" ", seq=1, icode=" ", x=0.0, y=0.0, z=0.0, occupancy=0.0, beta=0.0, footnote=0, segid="    ", element="  ", charge="  ", weight=0.0):
        self.recd = recd             #  1- 6
        self.serial = serial         #  7-11
        self.type = type             # 13-14
        self.loc = loc               # 15-16
        self.alt = alt               # 17
        self.res = res               # 18-21
        self.chain = chain           # 22
        self.seq = seq               # 23-26
        self.icode = icode           # 27
        self.x = x                   # 31-38
        self.y = y                   # 39-46
        self.z = z                   # 47-54
        self.occupancy = occupancy   # 55-60
        self.beta = beta             # 61-66
        self.footnote = footnote     # 68-70
        self.segid = segid           # 73-76
        self.element = element       # 77-78
        self.charge = charge         # 79-80
        self.weight = weight         # mass of atom assigned by read_pdb


# functions

def m_inout_read_pdb(filename):
    "Situs 2.x style pdb parser"
		
    pdb_list=[]
    f = open(filename, 'r')
    for line in f:
        test_recd=line[0:6]
    	   
        #if test_recd.upper().strip() == "ATOM" or test_recd.upper() == "HETATM":
        if ((test_recd.upper().strip() == "ATOM") or (test_recd.upper().strip() == "HETATM")):
            # create new atom and assign standard fields
            new_atom = PDB()
            new_atom.recd = test_recd
            try:
                new_atom.serial=int(line[5:11])
            except:
                pass
            try:
                new_atom.type=line[12:14]
            except:
                pass
            try:
                new_atom.loc=line[14:16]
            except:
                pass
            try:
                new_atom.alt=line[16]
            except:
                pass
            try:
                new_atom.res=line[17:21]
            except:
                pass
            try:
                new_atom.chain=line[21]
            except:
                pass
            try:
                new_atom.seq=int(line[22:26])
            except:
                pass
            try:
                new_atom.icode=line[26]
            except:
                pass
            try:
                new_atom.x=float(line[30:38])
            except:
                pass
            try:
                new_atom.y=float(line[38:46])
            except:
                pass
            try:
                new_atom.z=float(line[46:54])
            except:
                pass
            try:
                new_atom.occupancy=float(line[54:60])
            except:
                pass
            try:
                new_atom.beta=float(line[60:66])
            except:
                pass
            try:
                new_atom.footnote=int(line[67:70])
            except:
                pass
            try:
                new_atom.segid=line[72:76]
            except:
                pass
            try:
                new_atom.element=line[76:78]
            except:
                pass
            try:
                new_atom.charge=line[78:80]
            except:
                pass
            # compute robust Situs 2.x style mass weight and append atom
            if (new_atom.type == "QV" and new_atom.loc == "OL") or (new_atom.type == "QP" and new_atom.loc == "DB"):
                new_atom.weight=0.0 # codebook vectors have zero mass
            elif (new_atom.type == "DE" and new_atom.loc == "NS"):
                new_atom.weight=new_atom.occupancy # volumetric map density mass assigned from occupancy field
            elif (new_atom.type.lstrip()[0] == "H"):
                new_atom.weight=element_mass[0] # quick hack for hydrogens, may return false positives for Hg, Hf, Ho
            elif (new_atom.type.strip in element_name):
                new_atom.weight = element_mass[element_name.index(new_atom.type.strip)] # standard element
            elif (new_atom.type.strip()[0] in element_name):
                new_atom.weight = element_mass[element_name.index(new_atom.type.strip()[0])] # catches frequent elements, C N O P S
            else:
                #Slavica: Mass of other atoms (eg, 1H,2H,3H, 1HH1, 2HH1, 1HD2, 2HD2 etc) to be incorporated later
                #print "Warning: atom %i: unable to identify atom type %s, assigning carbon mass" % (len(pdb_list)+1,new_atom.type)
                new_atom.weight = element_mass[5]
            pdb_list.append(new_atom)
           
    f.close()
    return pdb_list



def m_inout_write_pdb(pdb_list,filename,remarks):
    "pdb writer, Situs 2.x convention, renumbering atoms"

    f = open(filename, 'w')
    f.write("REMARKS "+ remarks + "\n")

    for i in range(len(pdb_list)):
        f.write(pdb_list[i].recd.strip()[:6].ljust(6))
        # renumber atoms, ignoring .serial
        if ((i+1)//100000 == 0):
            f.write("%5d" % (i+1))
        else:
            f.write("%05d" % (i+1-((i+1)//100000)*100000))
        f.write(" ")
        f.write(pdb_list[i].type.strip()[:2].rjust(2))
        f.write(pdb_list[i].loc.strip()[:2].ljust(2))
        f.write(pdb_list[i].alt.strip()[:1].ljust(1))
        if (len(pdb_list[i].res)<4):
            f.write(pdb_list[i].res.strip().rjust(3))
            f.write(" ")
        else:
            f.write(pdb_list[i].res[:4])
        f.write(pdb_list[i].chain.strip()[:1].ljust(1))
        f.write("%4d" % pdb_list[i].seq)       
        f.write(pdb_list[i].icode.strip()[:1].ljust(1))
        f.write("   ")
        f.write("%8.3f" % pdb_list[i].x)
        f.write("%8.3f" % pdb_list[i].y)
        f.write("%8.3f" % pdb_list[i].z)
        f.write("%6.2f" % pdb_list[i].occupancy)
        f.write("%6.2f" % pdb_list[i].beta)       
        f.write(" ")
        f.write("%3d" % pdb_list[i].footnote)
        f.write("  ")
        f.write(pdb_list[i].segid.strip()[:4].ljust(4))
        f.write(pdb_list[i].element.strip()[:2].rjust(2))
        f.write(pdb_list[i].charge.strip()[:2].rjust(2))
        f.write("\n")
    f.close()

def m_inout_write_pdb_sampled(pdb_name,filename,step,thr_mass):
    "pdb writer, Situs 2.x convention, renumbering atoms"

    f = open(filename, 'w')

    pdb_list = m_inout_read_pdb(pdb_name)

    for i in range(0,len(pdb_list),step):
        if pdb_list[i].beta > thr_mass:
            f.write(pdb_list[i].recd.strip()[:6].ljust(6))
            # renumber atoms, ignoring .serial
            if ((i+1)//100000 == 0):
                f.write("%5d" % (i+1))
            else:
                f.write("%05d" % (i+1-((i+1)//100000)*100000))
            f.write(" ")
            f.write(pdb_list[i].type.strip()[:2].rjust(2))
            f.write(pdb_list[i].loc.strip()[:2].ljust(2))
            f.write(pdb_list[i].alt.strip()[:1].ljust(1))
            if (len(pdb_list[i].res)<4):
                f.write(pdb_list[i].res.strip().rjust(3))
                f.write(" ")
            else:
                f.write(pdb_list[i].res[:4])
            f.write(pdb_list[i].chain.strip()[:1].ljust(1))
            f.write("%4d" % pdb_list[i].seq)       
            f.write(pdb_list[i].icode.strip()[:1].ljust(1))
            f.write("   ")
            f.write("%8.3f" % pdb_list[i].x)
            f.write("%8.3f" % pdb_list[i].y)
            f.write("%8.3f" % pdb_list[i].z)
            f.write("%6.2f" % pdb_list[i].occupancy)
            f.write("%6.2f" % pdb_list[i].beta)       
            f.write(" ")
            f.write("%3d" % pdb_list[i].footnote)
            f.write("  ")
            f.write(pdb_list[i].segid.strip()[:4].ljust(4))
            f.write(pdb_list[i].element.strip()[:2].rjust(2))
            f.write(pdb_list[i].charge.strip()[:2].rjust(2))
            f.write("\n")
    f.close()

def m_inout_carbon_alphas(pdb_list):
    "return 0-based indices of carbon alphas in pdb_list"
    
    cas = []
    for i in range(len(pdb_list)):
        if (pdb_list[i].type == " C" and pdb_list[i].loc == "A "):
            cas.append(i)
    return cas

def m_inout_import_coords(filename): 
    "get x y z coords from PDB file"	
    coords = numpy.array([])
    pdb_list = m_inout_read_pdb(filename)
    for i in numpy.arange(len(pdb_list)):
        coords = numpy.append(coords,[pdb_list[i].x,pdb_list[i].y,pdb_list[i].z])
    return coords

def m_inout_import_bfact(filename): 
    "get Bfactor from PDB file"	
    bfact = numpy.array([])
    pdb_list = m_inout_read_pdb(filename)
    for i in numpy.arange(len(pdb_list)):
        bfact = numpy.append(bfact,[pdb_list[i].beta])
    return bfact
