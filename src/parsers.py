#imports
import re
import os
import pandas
import numpy
from io import StringIO
from ase import Atoms
from ase.data import atomic_masses_iupac2016
from ase.visualize import view
from ase.io import write
from rdkit import Chem
from ase.geometry.analysis import Analysis
from itertools import combinations
from openbabel import pybel
from data import Angstrom, eV

class ORCA_parser:
    '''
    ORCA_parcer object.
    The ORCA_Parcer object can extract information from an orca output file
    '''
    def __init__(self, filename):
        '''
        filename: str full path to the ORCA output file
        '''
        filename = open(str(filename), 'r')
        filename = filename.read()
        self.filename = filename

    def info(self):
        '''
        Function returns dict with information about functional, basis set and system charge
        '''
        pattern = r'!(.+?)\n'
        info = re.findall(pattern, self.filename, re.IGNORECASE)[0]
        functionals = ['\W*(b3lyp)\W', '\W*(revpbe)\W', '\W*(revpbe0)\W', '\W*(revpbe38)\W', '\W*(pbe)\W', '\W*(pbe0)\W', '\W*(blyp)\W', '\W*(m06)\W', '\W*(m06l)\W', '\W*(m062x)\W', '\W*(DLPNO-CCSD)\W']
        if any((functional := re.search(item, info, re.IGNORECASE)) for item in functionals):
            functional = functional.group(1)
        basis_sets = ['\W*(ma-def2-tzvpp)\W', '\W*(def2-tzvpp)\W', '\W*(cc-pvtz)\W']
        if any((basis := re.search(item, info, re.IGNORECASE)) for item in basis_sets):
            basis=basis.group(1)
        solvation=['cpcm\((.+?)\)']
        if any((solvation := re.search(item, info, re.IGNORECASE)) for item in solvation):
            solvation = solvation.group(1)
        else:
            solvation = 'gas'
        charge = re.search('Total Charge\s+\w+\s+\.+\s+(.+?)\n', self.filename, re.IGNORECASE).group(1)
        multiplicity = re.search('Multiplicity\s+\w+\s+\.+\s+(.+?)\n', self.filename, re.IGNORECASE).group(1)
        dipole = re.search('Magnitude \W(Debye)\W\s+\:\s+(.+?)\n', self.filename, re.IGNORECASE).group(2)
        temperature = re.findall('Temperature\s+\.+\s+(.+?)\w+\n', self.filename, re.IGNORECASE)
        info = {'functional': functional.upper(),
              'basis': basis.upper(),
              'charge': int(charge),
              'multiplicity': int(multiplicity),
              'solvation': solvation,
              'temperature': temperature[len(temperature)-1],
              'dipole': float('{:.2f}'.format(float(dipole)))}
        return info

    def geometry(self):
        '''
        Function returns molecule geometry
        '''
        pattern = r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n([\s\S]*?)\n-+\nCARTESIAN COORDINATES \(A\.U\.\)"
        geometry = re.findall(pattern,self.filename,re.DOTALL)#.group(1)
        table_assign = ['Element','X','Y','Z']
        geometry = pandas.read_csv(StringIO(geometry[len(geometry)-1]), names = table_assign, sep = r'\s+',engine = 'python')
        positions = [(geometry['X'][i], geometry['Y'][i], geometry['Z'][i]) for i in range(len(geometry['Element']))]
        geometry = Atoms([str(i) for i in geometry['Element']], positions=positions)
        return geometry

    def molar_mass(self):
        '''
        Returns molar mass of the molecule in g/mol
        '''
        molar_mass = [atomic_masses_iupac2016[i] for i in self.geometry().get_atomic_numbers()]
        return float(numpy.sum(molar_mass))

    def mol_vol(self):
        '''
        Returns molecular volume from CPCM model in m^3
        '''
        if self.info()['solvation'] != 'gas':
            pattern="GEPOL Volume\s+\.+\s+\d+.\d+"
            vol=re.findall(pattern, self.filename)
            vol=vol[len(vol)-1]
            pattern="\d+.\d+"
            vol=re.search(pattern, vol).group(0)
        else:
            vol = 'NaN'
        return float(vol)*Angstrom**3

    def visualise(self):
        view(self.geometry())

    def smiles(self):
        '''
        Returns SMILES of the molecule
        '''
        geo = self.geometry()

        #Number of atoms calculation
        symbols = geo.get_chemical_symbols()
        number_of_atoms = len(symbols)

        #Number of bonds calculation
        unique_symbols = []
        for i in symbols:
            if i not in unique_symbols:
                unique_symbols.append(i)
        Bonds = []
        for i in enumerate(unique_symbols):
            Bonds.append((i[1], i[1]))
        for i in combinations(unique_symbols, 2):
            Bonds.append(i)
        geo_analysis = Analysis(geo)
        at_combinations = []
        for i in Bonds:
            Bonds = geo_analysis.get_bonds(i[0], i[1], unique=True)
            if all(Bonds) == True:
                #BondValues = geo_analysis.get_values(Bonds)
                Bonds_list = [item for sublist in Bonds for item in sublist]
                at_combinations.append(Bonds_list)
        Bonds = sorted([item for sublist in at_combinations for item in sublist])
        number_of_bonds = len(Bonds)
        at_number = []
        for i in enumerate(Bonds):
            for j in enumerate(i[1]):
                if j[1] not in at_number:
                    at_number.append(j[1])
        at_number = sorted(at_number)
        #Bonds type
        bond_type = []
        for i in Bonds:
            if float('{:.2f}'.format(geo.get_distance(i[0], i[1]))) < 1.1:
                bond_type.append(1)
            elif (float('{:.2f}'.format(geo.get_distance(i[0], i[1]))) < 1.17) and (float('{:.2f}'.format(geo.get_distance(i[0], i[1]))) > 1.1):
                bond_type.append(3)
            elif (float('{:.2f}'.format(geo.get_distance(i[0], i[1]))) < 1.41) and (float('{:.2f}'.format(geo.get_distance(i[0], i[1]))) > 1.35):
                bond_type.append('ar')
            else:
                bond_type.append(1)
        #File generation
        f = open("gen_structure.mol2", "w")
        Molecule_section_line = """@<TRIPOS>MOLECULE\n*****\n{0} {1} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n"""
        f.write(Molecule_section_line.format(number_of_atoms, number_of_bonds))
        Atom_section_line = """{0} {1}\t\t\t{2}\t{3}\t{4} {5}\t{6}  {7}\n"""
        for i in range(len(symbols)):
            x, y, z = geo.positions[i]
            f.write(Atom_section_line.format(i+1, symbols[i], '{:.4f}'.format(x), '{:.4f}'.format(y), '{:.4f}'.format(z), symbols[i], 1, 'MOL'))
        f.write('@<TRIPOS>BOND\n')
        Bond_section_line = """{0}\t{1}\t{2}\t{3}\n"""
        for i in range(number_of_bonds):
            f.write(Bond_section_line.format(i + 1, Bonds[i][0] + 1, Bonds[i][1] + 1, bond_type[i]))
        f.close()
        mol = next(pybel.readfile('mol2', 'gen_structure.mol2'))
        smiles = mol.write(format="smi").split()[0].strip()
        os.remove('gen_structure.mol2')
        return smiles

    def thermodynamics(self):
        '''
        Returns dict with thermodynamics data of the molecule
        '''
        pattern = r'!(.+?)\n'
        info = re.findall(pattern, self.filename, re.IGNORECASE)[0]
        if re.search('FREQ', info, re.IGNORECASE) is not None:
            pattern = "Total thermal energy\s+(.+?)\s+Eh"
            zpe = re.search(pattern,self.filename,re.DOTALL).group(1)
            pattern = "Total Enthalpy\s+\.+\s+(.+?)\s+Eh"
            enthalpy = re.search(pattern,self.filename,re.DOTALL).group(1)
            pattern = "Final entropy term\s+\.+\s+(.+?)\s+Eh"
            entropy = re.search(pattern,self.filename,re.DOTALL).group(1)
            pattern = "Final Gibbs free energy\s+\.+\s+(.+?)\s+Eh"
            gibbs = re.search(pattern,self.filename,re.DOTALL).group(1)
            thermo_data = {'ZPE, eV':float('{:.2f}'.format(float(zpe)*eV)),
                         'H, eV':float('{:.2f}'.format(float(enthalpy)*eV)),
                         'TS, eV':float('{:.2f}'.format(float(entropy)*eV)),
                         'G, eV':float('{:.2f}'.format(float(gibbs)*eV))}
        elif re.search('FREQ', info, re.IGNORECASE) is None:
            pattern = 'FINAL SINGLE POINT ENERGY\s+(.+?)\s'
            zpe = re.findall(pattern, self.filename, re.DOTALL)
            zpe=zpe[len(zpe)-1]
            thermo_data = {'ZPE, eV':float('{:.2f}'.format(float(zpe)*eV)),
                         'H, eV': 'NaN',
                         'TS, eV': 'NaN',
                         'G, eV': 'NaN'}
        return thermo_data
