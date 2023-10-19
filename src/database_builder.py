#imports
import pandas
import os
import numpy
from parsers import ORCA_parser
from itertools import permutations
from itertools import combinations
from data import electron_charge, vacuum_dielectric_permitivity, solvent_refractive_index, solvent_dielectric_permitivity, J_to_eV, gas_constant, faraday, k_Boltzman, solvent_viscosity

def Unique_comb(a, b):
    unique_combinations = []
    permut1 = permutations(a, len(b))
    permut2 = permutations(b, len(a))
    for comb in permut1:
        zipped = zip(comb, b)
        if zipped not in unique_combinations:
            unique_combinations.append(list(zipped))
    for comb in permut2:
        zipped = zip(a, comb)
        if zipped not in unique_combinations:
            unique_combinations.append(list(zipped))
    return unique_combinations

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

def Drop_duplicated(database):
    duplicated = database.duplicated(keep = False)
    if any(duplicated) == True:
        index = [i[0] for i in enumerate(duplicated) if i[1] == True]
        duplicated = database.iloc[index].drop_duplicates().min()
        database = database.drop_duplicates(keep = False)
        database = pandas.concat([database, duplicated.to_frame().T], ignore_index = True)
    else:
        database = pandas.concat([database], ignore_index = True)
    return database

def Elements_from_SMILES(SMILES):
    return SMILES.split('.')

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

def Geo_comp(geo1, geo2):
    '''This function compares two geometries'''
    geo1 = ORCA_parser(geo1).geometry()
    geo2 = ORCA_parser(geo2).geometry()
    if len(geo1.positions) != len(geo2.positions):
        return False
    elif len(geo1.positions) == len(geo2.positions):
        res = []
        for i in enumerate(geo1.positions):
            tr = []
            for j in enumerate(i[1]):
                if geo1.positions[i[0]][j[0]] == geo2.positions[i[0]][j[0]]:
                    tr.append(True)
                else:
                    tr.append(False)
            res.append(tr)
    return all([all(j[1]) for j in enumerate(res)])

def Find_files(filename, search_path):
   '''This function finds all files with the name in directory'''
   result = []
   for root, dir, files in os.walk(search_path):
      if filename in files:
         result.append(os.path.join(root, filename))
   return result

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

def SMILES_analysis(dataframe):
    '''This function separates data by smiles in the right group'''
    index_molecules = []
    index_adiabatic = []
    index_radicals = []
    index_BDE = []
    for i in enumerate(dataframe['SMILES']):
        if (dataframe['Charge'][i[0]] == 0 and dataframe['Multiplicity'][i[0]] == 1 and ('.' not in i[1])):
            index_molecules.append(i[0])
        elif (dataframe['Charge'][i[0]] != 0 and dataframe['Multiplicity'][i[0]] >= 1):
            index_adiabatic.append(i[0])
        elif (dataframe['Charge'][i[0]] == 0 and dataframe['Multiplicity'][i[0]] == 2):
            index_radicals.append(i[0])
    for i in enumerate(dataframe['SMILES']):
        for j in dataframe['path'].iloc[index_adiabatic]:
            geo1 = dataframe['path'][i[0]]
            geo2 = j
            if (Geo_comp(geo1, geo2) == True and (dataframe['Charge'][i[0]] == 0 and dataframe['Multiplicity'][i[0]] == 1) and (i[0] not in index_BDE)):
                index_BDE.append(i[0])
    index_molecules = [i for i in index_molecules if i not in index_BDE]
    return [index_molecules, index_adiabatic, index_radicals, index_BDE]

def Electrochem_reactions(database):
    x = SMILES_analysis(database)
    database = database.drop(x[3]).reset_index(drop=True)
    all_comb = combinations(database.index, 2)
    all_comb = [i for i in all_comb]
    for i in enumerate(database.index):
        all_comb.append((i[0], i[0]))
    all_comb = [(i[0], i[1]) for i in all_comb if float('{:.3f}'.format(database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])) in database['Molar mass'].values]

    SMILES = pandas.Series([(database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + ' + ' +  database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb])

    molar_mass = [(database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] +  database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    functional = [(database['Functional'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Functional'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Functional'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    basis = [(database['Basis'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Basis'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Basis'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    spin_1 = numpy.array([(database['Multiplicity'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] - 1)/2 for i in all_comb])

    spin_2 = numpy.array([(database['Multiplicity'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]] - 1)/2 for i in all_comb])

    spin = numpy.array([(spin_1[i[0]] + spin_2[i[0]]) if (spin_1[i[0]] + spin_2[i[0]]) != 1 else (spin_1[i[0]] - spin_2[i[0]]) for i in enumerate(spin_1)])

    mult = (2*(spin)+ 1)

    mult_1 = (2*(spin_1)+ 1)

    mult_2 = (2*(spin_2)+ 1)

    charge = [(database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    charge_1 = [(database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb]

    charge_2 = [(database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    temp = [(database['Temperature, K'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Temperature, K'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Temperature, K'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    solvent = [(database['Solvent'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Solvent'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Solvent'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    ZPE = [(database['ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    H = [(database['H, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['H, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    TS = [(database['TS, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['TS, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    G = [(database['G, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['G, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    G_ZPE = [(database['G-ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['G-ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    dip = ['NaN' for i in enumerate(G_ZPE)]
    mol_vol = ['NaN' for i in enumerate(G_ZPE)]
    path = ['NaN' for i in enumerate(G_ZPE)]

    columns = ['SMILES', 'Molar mass', 'Functional', 'Basis', 'Multiplicity 1', 'Charge 1', 'Multiplicity 2', 'Charge 2', 'Multiplicity', 'Charge', 'Temperature, K', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'Molecular volume', 'path']

    df = pandas.DataFrame(numpy.array([SMILES, molar_mass, functional, basis, mult_1, charge_1, mult_2, charge_2, mult, charge, temp, solvent, ZPE, H, TS, G, G_ZPE, dip, mol_vol, path]).T, columns = columns)
    df = df.loc[(df['Charge'] == -1) & (df['Multiplicity'] == 2)]
    database = pandas.concat([database, df], ignore_index = True)

    columns = ['Reagent SMILES', 'Reagent multiplicity', 'Reagent charge', 'Reagent energy, eV', 'Product SMILES','Product multiplicity', 'Product charge', 'Product energy, eV', 'Reaction', 'Type', 'Mechanism', 'Free energy, eV']
    dataframe = pandas.DataFrame(columns = columns)
    electron = (0.5, -1) # (electron spin, charge)
    for i in enumerate(database.index):
        mult = database['Multiplicity'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])]
        charge = database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])]
        #spin change due to the electron transfer
        old_spin = (mult - 1)/2
        new_spin = numpy.array([(i[1] - electron[0]) if i[1] == 0.5 else (i[1] + electron[0]) for i in enumerate(old_spin)])
        new_mult = 2*new_spin + 1
        new_charge = numpy.array([(i[1] + electron[1]) for i in enumerate(charge)])

        reagent_smiles = database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]]) & (database['Molar mass'] == database['Molar mass'].iloc[i[0]]) & (database['Multiplicity'] == database['Multiplicity'].iloc[i[0]]) & (database['Charge'] == database['Charge'].iloc[i[0]])]

        product_smiles = database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]]) & (database['Molar mass'] == database['Molar mass'].iloc[i[0]]) & (database['Multiplicity'] == new_mult[i[0]]) & (database['Charge'] == new_charge[i[0]])]


        if len(reagent_smiles.index) != 0 and len(product_smiles.index) != 0:
            uni_comb = Unique_comb(reagent_smiles.index, product_smiles.index)

            for i in enumerate(uni_comb):
                if '+' in database['SMILES'].iloc[i[1][0][1]]:
                    react_type = 'dissoc.'
                    mech = 'concerted'
                elif '.' in database['SMILES'].iloc[i[1][0][1]]:
                    react_type = 'dissoc.'
                    mech = 'concerted'
                else:
                    react_type = 'adiab.'
                    mech = 'stepwise'
                data = {'Reagent SMILES': database['SMILES'].iloc[i[1][0][0]],
                        'Reagent multiplicity': database['Multiplicity'].iloc[i[1][0][0]],
                        'Reagent charge': database['Charge'].iloc[i[1][0][0]],
                        'Reagent energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][0]])),
                        'Product SMILES': database['SMILES'].iloc[i[1][0][1]],
                        'Product multiplicity': database['Multiplicity'].iloc[i[1][0][1]],
                        'Product charge': database['Charge'].iloc[i[1][0][1]],
                        'Product energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][1]])),
                        'Reaction': 'electrochem.',
                        'Type': react_type,
                        'Mechanism': mech,
                        'Free energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][1]] - database['G, eV'].iloc[i[1][0][0]]))}
                dataframe.loc[len(dataframe)] = data

    dataframe = dataframe.drop_duplicates(subset = ['Reagent SMILES','Product SMILES', 'Free energy, eV'])

    return dataframe


def Thermochem_reactions(database):
    x = SMILES_analysis(database)
    database = database.drop(x[3]).reset_index(drop=True)
    all_comb = combinations(database.index, 2)
    all_comb = [i for i in all_comb]
    for i in enumerate(database.index):
        all_comb.append((i[0], i[0]))
    sup_molar_mass = [float('{:.3f}'.format(database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])) for i in all_comb]

    all_comb = [(i[0], i[1]) for i in all_comb if float('{:.3f}'.format(database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])) in database['Molar mass'].values or float('{:.3f}'.format(database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])) in sup_molar_mass]

    SMILES = pandas.Series([(database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + ' + ' +  database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb])

    molar_mass = [(database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] +  database['Molar mass'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    functional = [(database['Functional'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Functional'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Functional'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    basis = [(database['Basis'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Basis'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Basis'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    spin_1 = numpy.array([(database['Multiplicity'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] - 1)/2 for i in all_comb])

    spin_2 = numpy.array([(database['Multiplicity'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]] - 1)/2 for i in all_comb])

    spin = numpy.array([(spin_1[i[0]] + spin_2[i[0]]) if (spin_1[i[0]] + spin_2[i[0]]) != 1 else (spin_1[i[0]] - spin_2[i[0]]) for i in enumerate(spin_1)])

    mult = (2*(spin)+ 1)

    mult_1 = (2*(spin_1)+ 1)

    mult_2 = (2*(spin_2)+ 1)

    charge = [(database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    charge_1 = [(database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb]

    charge_2 = [(database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    temp = [(database['Temperature, K'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Temperature, K'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Temperature, K'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    solvent = [(database['Solvent'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]]) for i in all_comb if (database['Solvent'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] == database['Solvent'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]])]

    ZPE = [(database['ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    H = [(database['H, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['H, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    TS = [(database['TS, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['TS, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    G = [(database['G, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['G, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    G_ZPE = [(database['G-ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[0]] + database['G-ZPE, eV'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])].iloc[i[1]]) for i in all_comb]

    dip = ['NaN' for i in enumerate(G_ZPE)]
    mol_vol = ['NaN' for i in enumerate(G_ZPE)]
    path = ['NaN' for i in enumerate(G_ZPE)]

    columns = ['SMILES', 'Molar mass', 'Functional', 'Basis', 'Multiplicity 1', 'Charge 1', 'Multiplicity 2', 'Charge 2', 'Multiplicity', 'Charge', 'Temperature, K', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'Molecular volume', 'path']

    df = pandas.DataFrame(numpy.array([SMILES, molar_mass, functional, basis, mult_1, charge_1, mult_2, charge_2, mult, charge, temp, solvent, ZPE, H, TS, G, G_ZPE, dip, mol_vol, path]).T, columns = columns)
    #df = df.loc[(df['Charge'] != -1) & (df['Multiplicity'] != 2)]
    database = pandas.concat([database, df], ignore_index = True)

    columns = ['Reagent SMILES', 'Reagent multiplicity', 'Reagent charge', 'Reagent energy, eV', 'Product SMILES','Product multiplicity', 'Product charge', 'Product energy, eV', 'Reaction', 'Type', 'Mechanism', 'Free energy, eV']
    dataframe = pandas.DataFrame(columns = columns)

    for i in enumerate(database.index):
        mult = database['Multiplicity'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])]
        charge = database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])]


        reagent_smiles = database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]]) & (database['Molar mass'] == database['Molar mass'].iloc[i[0]]) & (database['Multiplicity'] == database['Multiplicity'].iloc[i[0]]) & (database['Charge'] == database['Charge'].iloc[i[0]])]

        product_smiles = database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]]) & (database['Molar mass'] == database['Molar mass'].iloc[i[0]]) & (database['Charge'] == database['Charge'].iloc[i[0]])]


        if len(reagent_smiles.index) != 0 and len(product_smiles.index) != 0:
            uni_comb = Unique_comb(reagent_smiles.index, product_smiles.index)

            for i in enumerate(uni_comb):
                if '+' in database['SMILES'].iloc[i[1][0][1]] or '+' in database['SMILES'].iloc[i[1][0][0]]:
                    react_type = 'dissoc.'
                    mech = 'concerted'
                elif '.' in database['SMILES'].iloc[i[1][0][1]] or '.' in database['SMILES'].iloc[i[1][0][0]]:
                    react_type = 'dissoc.'
                    mech = 'concerted'
                else:
                    react_type = 'adiab.'
                    mech = 'stepwise'
                data = {'Reagent SMILES': database['SMILES'].iloc[i[1][0][0]],
                        'Reagent multiplicity': database['Multiplicity'].iloc[i[1][0][0]],
                        'Reagent charge': database['Charge'].iloc[i[1][0][0]],
                        'Reagent energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][0]])),
                        'Product SMILES': database['SMILES'].iloc[i[1][0][1]],
                        'Product multiplicity': database['Multiplicity'].iloc[i[1][0][1]],
                        'Product charge': database['Charge'].iloc[i[1][0][1]],
                        'Product energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][1]])),
                        'Reaction': 'thermochem.',
                        'Type': react_type,
                        'Mechanism': mech,
                        'Free energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][1]] - database['G, eV'].iloc[i[1][0][0]]))}
                dataframe.loc[len(dataframe)] = data

    dataframe = dataframe.drop_duplicates(subset = ['Free energy, eV'], keep = 'first')
    dataframe = dataframe[dataframe['Free energy, eV'] < 0].reset_index()
    dataframe = dataframe.drop(columns = ['index'])
    return dataframe



#Main code
#input path to the calculations results by user
print('Enter path to the folder with your calculations')
input_path = input('Input path:\t')
files = Find_files("out.txt", input_path)

Database = Data_collection(files)
Newdatabase = pandas.read_csv(Database[0])
#Newdatabase = pandas.read_csv('Data/Electrocarboxylation/Database.csv')
Electrochem = Electrochem_reactions(Newdatabase)
#Reactions = Thermochem_reactions(Newdatabase)
#Database = pandas.concat([Electrochem, Reactions], ignore_index = True)
#Reaction = Reaction_constructor(Newdatabase, SMILES_analysis(Newdatabase))
