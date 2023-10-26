#imports
import os
import sys
sys.path.append('src')

import pandas
import numpy
from parsers import ORCA_parser
from data import electron_charge, vacuum_dielectric_permitivity, solvent_refractive_index, solvent_dielectric_permitivity, J_to_eV, gas_constant, k_Boltzman, solvent_viscosity

def Diffusion(solvent, temperature, mol_vol):
    diff=k_Boltzman*temperature/(6*numpy.pi*solvent_viscosity[solvent]*(3*mol_vol/(4*numpy.pi))**(1/3))
    return diff*10000

def Solvent_reorganization(mol_vol, solvent):
    mol_rad = (3*mol_vol/(4*numpy.pi))**(1/3)
    solvent_reorg_eng = (electron_charge**2)/(8*numpy.pi*vacuum_dielectric_permitivity)*(1/(solvent_refractive_index[solvent]**2)-1/solvent_dielectric_permitivity[solvent])*(1/mol_rad-1/(2*mol_rad))
    solvent_reorg_eng = solvent_reorg_eng*J_to_eV
    return float('{:.2f}'.format(solvent_reorg_eng))

def Collision_factor(temperature, molar_mass):
    Z = numpy.sqrt(gas_constant*temperature/2/numpy.pi/(molar_mass/1000))*100
    return Z


def Find_files(filename, search_path):
   '''This function finds all files with the name in directory'''
   result = []
   for root, dir, files in os.walk(search_path):
      if filename in files:
         result.append(os.path.join(root, filename))
   return result

def SMILES_modifier(mult, charge):
    if mult == 2 and charge == 0:
        mod = '[*]'
        return mod
    elif mult == 2 and charge < 0:
        charge = abs(int(charge))
        mod = '[*{0}-]'.format(charge)
        return mod
    elif mult == 1 and charge < 0:
        charge = abs(int(charge))
        mod = '[{0}-]'.format(charge)
        return mod
    elif charge > 0:
        charge = abs(int(charge))
        mod = '[{0}+]'.format(charge)
        return mod
    else:
        return ''

def Data_collection(files):
    columns = ['SMILES', 'Molar mass', 'Functional', 'Basis', 'Multiplicity', 'Charge', 'Temperature, K', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'Molecular volume', 'path']
    database = pandas.DataFrame(columns=columns)
    for i in enumerate(files):
        file_data = ORCA_parser(i[1])
        if file_data.thermodynamics()['G, eV'] != 'NaN':
            par = file_data.thermodynamics()['G, eV'] - file_data.thermodynamics()['ZPE, eV']
        else:
            par = 'NaN'
        if file_data.info()['temperature'] != None:
            temp = file_data.info()['temperature']
        else:
            temp = 'NaN'
        SMILES = file_data.smiles() + SMILES_modifier(file_data.info()['multiplicity'], file_data.info()['charge'])
        dataframe = {'SMILES': SMILES,
                     'Molar mass': '{:.3f}'.format(file_data.molar_mass()),
                     'Functional': file_data.info()['functional'],
                     'Basis': file_data.info()['basis'],
                     'Charge': file_data.info()['charge'],
                     'Multiplicity': file_data.info()['multiplicity'],
                     'Temperature, K': temp,
                     'Solvent': file_data.info()['solvation'],
                     'ZPE, eV': file_data.thermodynamics()['ZPE, eV'],
                     'H, eV': file_data.thermodynamics()['H, eV'],
                     'TS, eV': file_data.thermodynamics()['TS, eV'],
                     'G, eV': file_data.thermodynamics()['G, eV'],
                     'G-ZPE, eV': par,
                     'Dipole, D': file_data.info()['dipole'],
                     'Molecular volume': file_data.mol_vol(),
                     'path': i[1]}
        database.loc[len(database)] = dataframe
    output_path = input('Output path:\t') + '/' + 'Database.csv'
    return output_path, database.to_csv(output_path)


#Main code
#input path to the calculations results by user
print('Enter path to the folder with your calculations')
input_path = input('Input path:\t')
files = Find_files("out.txt", input_path)

Database = Data_collection(files)
