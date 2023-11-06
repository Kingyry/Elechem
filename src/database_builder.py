#imports
import os
import sys
sys.path.append('src')

import pandas
import numpy
from parsers import ORCA_parser
from data import k_Boltzman, solvent_viscosity

def Diffusion(solvent, temperature, mol_vol):
    diff = k_Boltzman*temperature/(6*numpy.pi*solvent_viscosity[solvent]*(3*mol_vol/(4*numpy.pi))**(1/3))
    return diff*10000


def Find_files(filename, search_path):
   '''This function finds all files with the name in directory'''
   result = []
   for root, dir, files in os.walk(search_path):
      if filename in files:
         result.append(os.path.join(root, filename))
   return result

def Data_collection(files):
    columns = ['SMILES', 'Molar mass', 'Functional', 'Basis', 'Multiplicity', 'Charge', 'Temperature, K', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'Molecular volume', 'Diffusion, cm2/s', 'path']
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
        ele = ['[Cl]', '[Br]', '[I]', '[Fe]', '[O]', '[CH2]']
        for k, j in enumerate(ele):
            if j in file_data.smiles() and file_data.info()['charge'] < 0:
                smiles = file_data.smiles().replace(j, '[' + j.strip("[]") + str(file_data.info()['charge']) + ']')
                break
            elif j in file_data.smiles() and file_data.info()['charge'] > 0:
                smiles = file_data.smiles().replace(j, '[' + j.strip("[]") + '+' + str(file_data.info()['charge']) + ']')
                break
            else:
                smiles = file_data.smiles()
        SMILES = smiles
        #Neweton-Stokes equation for diffusion
        if file_data.info()['solvation'] != 'gas':
            D = Diffusion(file_data.info()['solvation'], temp, file_data.mol_vol())
        else:
            D = 'NaN'
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
                     'Diffusion, cm2/s': D,
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
