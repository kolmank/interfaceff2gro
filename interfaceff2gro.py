#!/usr/bin/env python

import MDAnalysis
import sys
import subprocess
import os

def atomtypes_from_ff(l):
    '''Extracts atomtypes from .prm charmm forcefield file. l is argument being a list containing the lines of forcefield file.'''
    f,i,s=float,int,str
    g = []
    for m in range(len(l)):
        if len(l[m])>0 and len(l[m])<15 and s(l[m][:10]) == '!INTERFACE':
            n = i(m)
            n += 1
            while len(l[n])>2:
                 if not l[n][:1] == '!':
                     h = s(l[n][:40])                 
                     h = h.strip()
                     h = h.split()
                     if len(h) == 4 and h[1] == '0.0':
                         g.append([s(h[0]),f(h[1]),f(h[2]),f(h[3])])
                 n +=1
    return g

def atomtypes_from_lammps(l):
    '''Extracts atomtypes from .data lammps file. l is argument being a list containing the lines of .data file.'''
    i,s = int,str
    g = []
    for m in range(len(l)):
        if len(l[m])>0 and l[m][:11] == 'Pair Coeffs':
            n = i(m)
            n += 2
            while len(l[n])>1:
                 g.append([(s(l[n][37:]).strip()).upper()])
                 n +=1
    return g

def bonds_from_ff(l):
    '''Extracts bondtypes from .prm charmm forcefield file. l is argument being a list containing the lines of forcefield file.'''
    f,i,s=float,int,str
    g = []
    for m in range(len(l)):
        if len(l[m])>0 and len(l[m])<15 and s(l[m][:10]) == '!INTERFACE':
            n = i(m)
            n += 1
            while len(l[n])>2:
                 if not l[n][:1] == '!':
                     h = s(l[n][:31])                 
                     h = h.strip()
                     h = h.split()
                     if len(h) == 4 and not h[1] == '0.0':
                         g.append([s(h[0]),s(h[1]),f(h[2]),f(h[3])])
                 n +=1
    return g

def bondtypes_from_lammps(l):
    '''Extracts bondtypes from .data lammps file. l is argument being a list containing the lines of .data file.'''
    i,s = int,str
    g = []
    for m in range(len(l)):
        if len(l[m])>0 and l[m][:11] == 'Bond Coeffs':
            n = i(m)
            n += 2
            while len(l[n])>1:
                 h = s(l[n][50:])
                 h = h.strip()
                 h = h.upper()
                 h = h.split('-')
                 g.append(h)
                 n +=1
    return g

def angles_from_ff(l):
    '''Extracts angletypes from .prm charmm forcefield file. l is argument being a list containing the lines of forcefield file.'''
    f,i,s=float,int,str
    g = []
    for m in range(len(l)):
        if len(l[m])>0 and len(l[m])<15 and s(l[m][:10]) == '!INTERFACE':
            n = i(m)
            n += 1
            while len(l[n])>2:
                 if not l[n][:1] == '!':
                     h = s(l[n][:36])                 
                     h = h.strip()
                     h = h.split()
                     if len(h) == 5 and not h[0] == 'X':
                         g.append([s(h[0]),s(h[1]),s(h[2]),f(h[3]),f(h[4])])
                 n +=1
    return g


def angletypes_from_lammps(l):
    '''Extracts angletypes from .data lammps file. l is argument being a list containing the lines of .data file.'''
    i,s = int,str
    g = []
    for m in range(len(l)):
        if len(l[m])>0 and l[m][:12] == 'Angle Coeffs':
            n = i(m)
            n += 2
            while len(l[n])>1:
                 h = s(l[n][50:])
                 h = h.strip()
                 h = h.upper()
                 h = h.split('-')
                 g.append(h)
                 n +=1
    return g

def create_lammps_data_file():
    '''Creates lampps *.data file from Materials studio .car .mdf and .frc file. To do this a utlity from Lammps package is used. J is the name of *.car file (e.g., something.car). The car and mdf file needs to have the same name. The frc file need to be located in ./msi2lmp/frc_files. '''
    bash_command = os.path.dirname(sys.argv[0]) + '/msi2lmp/src/msi2lmp.exe ' + sys.argv[1][:-4] + ' -class 2 -frc ' + os.path.dirname(sys.argv[0]) + '/msi2lmp/frc_files/pcff_interface.frc -ignore'
#    g = subprocess.Popen(bash_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    g = subprocess.Popen(bash_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(g.stdout.read())

    
def clean_lammps_data_file():
    '''Removes .data file created by the lammps utility software '''
    bash_command = 'rm ' + os.getcwd() + '/' + sys.argv[1][:-4] + '.data'
    g = subprocess.Popen(bash_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    print('Cleaning temporary files.')
    print(g.stdout.read())

# Create a data file from car and mdf files.
try:
    create_lammps_data_file();
except:
    print('First command-line argument is not a proper .car file, or .mdf file is missing/having different name than the .car file, or msi2lmp.exe is not located in ' + os.path.dirname(sys.argv[0]) + '/msi2lmp/src.')
    print('Usage: ./interfaceff2gro.py nameoffile.car.')
    sys.exit()

# Open the created data file and create an universe.
try:
    u = MDAnalysis.Universe(sys.argv[1][:-4] + '.data');
except:
    print('The .data file from .car and .mdf files has not been created.')
    sys.exit()

# Open the forcefield file
try:
    forcefield = open(os.path.dirname(sys.argv[0]) + '/forcefield/charmm27_interface_v1_5.prm','r').readlines();
except:
    print('Frocefield .prm file missing. It should be located in ' + os.path.dirname(sys.argv[0]) + '/forcefield.')
    sys.exit()
    
# Read the created lammps data file
lammps = open(sys.argv[1][:-4] + '.data','r').readlines();

# Extract all informations from .prm forcefield file. Dihedrals and Impropers are zero so are not extracted.
atom_ff = atomtypes_from_ff(forcefield)

# Create dicts having episilion and sigma values related to atomtype
atom_ff_sigma, atom_ff_epsilion = {},{};

for m in atom_ff:
    atom_ff_sigma[m[0]]=m[3]
    atom_ff_epsilion[m[0]]=m[2]

# Extract all informations from .prm forcefield file. Dihedrals and Impropers are zero so are not extracted.
bond_ff = bonds_from_ff(forcefield);

angle_ff = angles_from_ff(forcefield);

# Extract types information from .data file
atom_lammps = atomtypes_from_lammps(lammps);

# Create a dict correlating ordr number from .data file with typename
atom_lammps_dict = {}
for i in range(len(atom_lammps)):
    atom_lammps_dict.update({i+1:atom_lammps[i][0]})
    
# Extract types information from .data file
bond_lammps = bondtypes_from_lammps(lammps);

angle_lammps = angletypes_from_lammps(lammps);


# Dictionary for atom numbers
element2atomNumber= {}
element2atomNumber['H']='1'
element2atomNumber['He']='2'
element2atomNumber['C']='6'
element2atomNumber['N']='7'
element2atomNumber['O']='8'
element2atomNumber['F']='9'
element2atomNumber['Ne']='10'
element2atomNumber['Na']='11'
element2atomNumber['Mg']='12'
element2atomNumber['Al']='13'
element2atomNumber['P']='15'
element2atomNumber['S']='16'
element2atomNumber['Cl']='17'
element2atomNumber['K']='19'
element2atomNumber['Ca']='20'
element2atomNumber['Fe']='26'
element2atomNumber['Zn']='30'
element2atomNumber['Br']='35'
element2atomNumber['I']='53'
element2atomNumber['Cs']='55'
element2atomNumber['Si']='14'
element2atomNumber['Ag']='47'
element2atomNumber['Au']='79'
element2atomNumber['Cu']='29'
element2atomNumber['Ni']='28'
element2atomNumber['Pb']='82'
element2atomNumber['Pd']='46'
element2atomNumber['Pt']='78'

# Dictionary for masses
mass2element= {}
mass2element[1.007970]='H'
mass2element[1.008000]='H'
mass2element[15.999400]='O'
mass2element[40.080000]='Ca'
mass2element[22.989770]='Na'
mass2element[26.981530]='Al'
mass2element[26.982000]='Al'
mass2element[28.086000]='Si'
mass2element[30.973800]='P'
mass2element[39.098300]='K'
mass2element[107.868000]='Ag'
mass2element[196.967000]='Au'
mass2element[63.546000]='Cu'
mass2element[58.710000]='Ni'
mass2element[207.200000]='Pb'
mass2element[106.400000]='Pd'
mass2element[195.090000]='Pt'

# conversion constant between kcal and kJ
kcal2kJ = 4.184

# Formatting of a .itp file.
cformat_atomtypes = '%-s\t%-s\t%8.5f\t%-.6f\t%s\t%.12f\t%.6f \n'
cformat_atomtypes_2 = '%-s\t%-s\t%8.6f\t%-.6f\t%s\t%.12f\t%.6f \n'
cformat_bondtypes = '%-s\t%-s\t%d\t%-.4f\t%-.2f\n'
cformat_angletypes = '%-s\t%-s\t%-s\t%d\t%-.2f\t%-.4f\t%.1f\t%.1f\n'
cformat_molecule = '%-5s%4d\n'
cformat_atoms = '%6d%11s%7d%7s%7s%7d%11.3f%11.3f   ;\n'
cformat_bonds = '%5d%6d%6d \n'
cformat_angles = '%5d%6d%6d%6d \n'
gro_format = '%5d%5s%5s%5d%8.3f%8.3f%8.3f\n'

# Check which atomtypes are present in universe
atomtypes_universe = []
for m in set(zip(u.atoms.types,u.atoms.masses)):
    atomtypes_universe.append(m)
    
# Check which bondtypes are present in universe
bondtypes_universe = bond_lammps

# Check which angletypes are present in universe
angletypes_universe = angle_lammps

# Check atoms in universe
atoms_universe = []
for m in zip(u.atoms.indices,u.atoms.types,u.atoms.masses,u.atoms.charges,u.atoms.masses):
    atoms_universe.append(m)
    
# Check bonds in universe
bonds_universe = []
try:
    for m in u.bonds.indices:
        bonds_universe.append(m)
except:
    print('Bonds between atoms have not been found.')
    
# Check angles in universe
angles_universe = []
try:
    for m in u.angles.indices:
        angles_universe.append(m)
except:
    print('Angles between atoms have not been found.')

# Atom properties for gro file
atoms_universe_gro = []
for m in zip(u.atoms.indices,u.atoms.masses,u.atoms.positions):
    atoms_universe_gro.append(m)
  
# Write a gro file
gro_file = open(sys.argv[1][:-4]+'.gro','w')

gro_file.write('Periodic slab: SURF, t= 0.0\n')
gro_file.write('%6d\n'%(atoms_universe_gro[-1][0]+1))

for m in atoms_universe_gro:
    gro_file.write(gro_format%(1,'SURF',mass2element[m[1]],m[0]+1,m[2][0]/10,m[2][1]/10,m[2][2]/10))

gro_file.write('%10.5f%10.5f%10.5f\n'%(u.dimensions[0]/10,u.dimensions[1]/10,u.dimensions[2]/10))
gro_file.close()

# Writing itp file
itp_file = open(sys.argv[1][:-4]+'.itp','w')
#itp_file = open('test.itp','w')

# writing atomtype section
itp_file.write('[ atomtypes ]\n')
itp_file.write(';name	at.num	mass	charge	ptype	sigma	epsilon\n')
for m in atomtypes_universe:
    if m[1] > 10:
        itp_file.write(cformat_atomtypes%(atom_lammps_dict[int(m[0])],element2atomNumber[mass2element[m[1]]],m[1],0,'A',2.0**(-1.0/6.0)*0.2*atom_ff_sigma[atom_lammps_dict[int(m[0])]],-kcal2kJ*atom_ff_epsilion[atom_lammps_dict[int(m[0])]]))
    else:
        itp_file.write(cformat_atomtypes_2%(atom_lammps_dict[int(m[0])],element2atomNumber[mass2element[m[1]]],m[1],0,'A',2.0**(-1.0/6.0)*0.2*atom_ff_sigma[atom_lammps_dict[int(m[0])]],-kcal2kJ*atom_ff_epsilion[atom_lammps_dict[int(m[0])]]))

# writing bondtype section
if len(bondtypes_universe) > 0:
    itp_file.write('\n')
    itp_file.write('[ bondtypes ]\n')
    itp_file.write('; i	j	func	b0	kb\n')
    for m in bondtypes_universe:
        for n in bond_ff:
            if (m[0] == n[0] and m[1] == n[1]) or (m[1] == n[0] and m[0] == n[1]):
                itp_file.write(cformat_bondtypes%(n[0],n[1],1,n[3]/10,n[2]*200*kcal2kJ))

# writing angletype section
if len(angletypes_universe) > 0:
    itp_file.write('\n')
    itp_file.write('[ angletypes ]\n')
    itp_file.write('; i	j	k	func	th0	cth\n')
    for m in angletypes_universe:
        for n in angle_ff:
            if (m[0] == n[0] and m[1] == n[1]and m[2] == n[2]) or (m[2] == n[0] and m[0] == n[2] and m[1] == n[1]):
                itp_file.write(cformat_angletypes%(n[0],n[1],n[2],5,n[4],n[3]*2*kcal2kJ,0,0))

# writing molecule section
itp_file.write('\n')
itp_file.write('[ moleculetype ]\n')
itp_file.write('; molname    nrexcl\n')
itp_file.write(cformat_molecule%('SURF',3.0))

# writing atom section
itp_file.write('\n')
itp_file.write('[ atoms ]\n')
itp_file.write(';   nr       type  resnr residue  atom   cgnr     charge       mass\n')
for m in atoms_universe:
    itp_file.write(cformat_atoms%(m[0]+1,atom_lammps_dict[int(m[1])],1,'SURF',mass2element[m[2]],m[0]+1,m[3],m[4]))

# writing bond section
if len(bonds_universe) > 0:
    itp_file.write('\n')
    itp_file.write('[ bonds ]\n')
    itp_file.write(';  ai    aj funct            c0            c1            c2            c3\n')
    for m in bonds_universe:
        itp_file.write(cformat_bonds%(m[0]+1,m[1]+1,1))
  
# writing angle section
if len(angles_universe) > 0:
    itp_file.write('\n')
    itp_file.write('[ angles ]\n')
    itp_file.write(';  ai    aj    ak funct            c0            c1            c2            c3\n')
    for m in angles_universe:
        itp_file.write(cformat_angles%(m[0]+1,m[1]+1,m[2]+1,5))

# Wriritng position restraints
itp_file.write('\n')
itp_file.write('; Include Position restraint file\n')
itp_file.write('#ifdef POSRES\n')
itp_file.write('#include "posre_SURF.itp"\n')
itp_file.write('#endif\n')
itp_file.write('\n')

# Closing file
itp_file.close()

# Writing top file
top_file = open('system.top','w')
top_file.write('#include "charmm27.ff/forcefield.itp"\n')
top_file.write('#include "%s.itp"\n'%sys.argv[1][:-4])
top_file.write('#include "charmm27.ff/ions.itp"\n')
top_file.write('#include "charmm27.ff/tip4p.itp"\n')
top_file.write('\n')
top_file.write('[ system ]\n')
top_file.write('Surface\n')
top_file.write('\n')
top_file.write('[ molecules ]\n')
top_file.write('SURF 1 ;1 periodic slab\n')
top_file.write('\n')
top_file.close()

# Removing the created data file
clean_lammps_data_file()

# Printing output to terminal
print('Three files have been created: %s.itp, %s.gro and system.top.'%(sys.argv[1][:-4],sys.argv[1][:-4]))
