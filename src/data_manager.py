#imports
import sys
sys.path.append('src')

import pandas
import numpy

from parsers import ORCA_parser
from itertools import permutations
from itertools import combinations
from data import electron_charge, vacuum_dielectric_permitivity, solvent_refractive_index, solvent_dielectric_permitivity, J_to_eV, gas_constant, faraday, Avogadro, k_Boltzman, Planck

def Choose_reactions(path):
    database = pandas.read_csv(path)
    print('Input index number (1 2 3.. n..)')
    Reaction_indexes = [int(i) for i in input("Enter the reaction index numbers:\t").split()]
    database = database.iloc[Reaction_indexes].reset_index(drop = True)
    database = database.drop(columns = ['Unnamed: 0'])
    return database.to_csv(path)

def Drop_reactions(path):
    database = pandas.read_csv(path)
    print('Input index number (1 2 3.. n..)')
    Reaction_indexes = [int(i) for i in input("Enter the reaction index numbers:\t").split()]
    database = database.drop(index = Reaction_indexes).reset_index(drop = True)
    database = database.drop(columns = ['Unnamed: 0'])
    return database.to_csv(path)

def Reduction_potential(path1):
    path2 = input('Input path to your reference electrode database:\t')
    database1 = pandas.read_csv(path1)
    database2 = pandas.read_csv(path2)
    red_pot = []
    for i in enumerate(database1.index):
        if database1['Reaction'].iloc[i[0]] == 'electrochem.':
            charge = database1['Reagent charge'].iloc[i[0]] - database1['Product charge'].iloc[i[0]]
            Reduction_potential = (-database1['Free energy, eV'].iloc[i[0]]/abs(charge)) + database2['Free energy, eV'].loc[(database2['Functional'] == database1['Functional'].iloc[i[0]]) & (database2['Basis'] == database1['Basis'].iloc[i[0]]) & (database2['Solvent'] == database1['Solvent'].iloc[i[0]])].values[0]
            red_pot.append(Reduction_potential)
        else:
            red_pot.append('NaN')
    database1['Reduction potential, V'] = red_pot
    database1 = database1.drop(columns = ['Unnamed: 0'])
    return database1.to_csv(path)

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
            if (Geo_comp(geo1, geo2) == True and (dataframe['Charge'][i[0]] == 0) and (i[0] not in index_BDE)):# and dataframe['Multiplicity'][i[0]] == 1
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

    all_comb = [(i[0], i[1]) for i in all_comb if (database['Functional'].iloc[i[0]] == database['Functional'].iloc[i[1]] and database['Basis'].iloc[i[0]] == database['Basis'].iloc[i[1]] and database['Solvent'].iloc[i[0]] == database['Solvent'].iloc[i[1]])]
    all_comb = [(i[0], i[1]) for i in all_comb if float('{:.3f}'.format(database['Molar mass'].iloc[i[0]] + database['Molar mass'].iloc[i[1]])) in database['Molar mass'].values]
    SMILES = pandas.Series([(database['SMILES'].iloc[i[0]] + '.' +  database['SMILES'].iloc[i[1]]) for i in all_comb])
    molar_mass = [(database['Molar mass'].iloc[i[0]] +  database['Molar mass'].iloc[i[1]]) for i in all_comb]
    functional = [(database['Functional'].iloc[i[0]]) for i in all_comb]
    basis = [(database['Basis'].iloc[i[0]]) for i in all_comb]
    spin_1 = numpy.array([(database['Multiplicity'].iloc[i[0]] - 1)/2 for i in all_comb])
    spin_2 = numpy.array([(database['Multiplicity'].iloc[i[1]] - 1)/2 for i in all_comb])
    spin = numpy.array([(spin_1[i[0]] + spin_2[i[0]]) if (spin_1[i[0]] + spin_2[i[0]]) != 1 else (spin_1[i[0]] - spin_2[i[0]]) for i in enumerate(spin_1)])
    mult = (2*(spin)+ 1)
    charge = [(database['Charge'].iloc[i[0]] + database['Charge'].iloc[i[1]]) for i in all_comb]
    temp = [(database['Temperature, K'].iloc[i[0]]) for i in all_comb]
    solvent = [(database['Solvent'].iloc[i[0]]) for i in all_comb]
    ZPE = [(database['ZPE, eV'].iloc[i[0]] + database['ZPE, eV'].iloc[i[1]]) for i in all_comb]
    H = [(database['H, eV'].iloc[i[0]] + database['H, eV'].iloc[i[1]]) for i in all_comb]
    TS = [(database['TS, eV'].iloc[i[0]] + database['TS, eV'].iloc[i[1]]) for i in all_comb]
    G = [(database['G, eV'].iloc[i[0]] + database['G, eV'].iloc[i[1]]) for i in all_comb]
    G_ZPE = [(database['G-ZPE, eV'].iloc[i[0]] + database['G-ZPE, eV'].iloc[i[1]]) for i in all_comb]

    dip = ['NaN' for i in enumerate(G_ZPE)]
    mol_vol = ['NaN' for i in enumerate(G_ZPE)]
    path = ['NaN' for i in enumerate(G_ZPE)]

    columns = ['SMILES', 'Molar mass', 'Functional', 'Basis', 'Multiplicity', 'Charge', 'Temperature, K', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'Molecular volume', 'path']

    df = pandas.DataFrame(numpy.array([SMILES, molar_mass, functional, basis, mult, charge, temp, solvent, ZPE, H, TS, G, G_ZPE, dip, mol_vol, path]).T, columns = columns)
    if df.empty == False:
        df = df.loc[(df['Charge'] == -1) & (df['Multiplicity'] == 2)]
        database = pandas.concat([database, df], ignore_index = True)

    columns = ['Reagent SMILES', 'Reagent multiplicity', 'Reagent charge', 'Reagent energy, eV', 'Product SMILES','Product multiplicity', 'Product charge', 'Product energy, eV', 'Functional', 'Basis', 'Reaction SMILES', 'Solvent', 'Temperature, K', 'Reaction', 'Type', 'Mechanism', 'Free energy, eV']
    dataframe = pandas.DataFrame(columns = columns)
    for i in enumerate(database.index):
        mult = database['Multiplicity'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])]
        charge = database['Charge'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]])]

        reagent_smiles = database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]]) & (database['Molar mass'] == database['Molar mass'].iloc[i[0]]) & (database['Multiplicity'] == database['Multiplicity'].iloc[i[0]]) & (database['Charge'] == database['Charge'].iloc[i[0]])]

        product_smiles = database['SMILES'].loc[(database['Functional'] == database['Functional'].iloc[i[0]]) & (database['Basis'] == database['Basis'].iloc[i[0]]) & (database['Solvent'] == database['Solvent'].iloc[i[0]]) & (database['Molar mass'] == database['Molar mass'].iloc[i[0]]) & (database['Multiplicity'] != database['Multiplicity'].iloc[i[0]]) & (database['Charge'] != database['Charge'].iloc[i[0]])]

        if len(reagent_smiles.index) != 0 and len(product_smiles.index) != 0:
            uni_comb = Unique_comb(reagent_smiles.index, product_smiles.index)

            for i in enumerate(uni_comb):
                if '.' in database['SMILES'].iloc[i[1][0][1]] and database['path'].iloc[i[1][0][1]] != 'NaN':
                    react_type = 'adiab.'
                    mech = 'concerted'
                elif '.' in database['SMILES'].iloc[i[1][0][1]] and database['path'].iloc[i[1][0][1]] == 'NaN':
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
                        'Functional': database['Functional'].iloc[i[1][0][0]],
                        'Basis': database['Basis'].iloc[i[1][0][0]],
                        'Reaction SMILES': database['SMILES'].iloc[i[1][0][0]] + '>>' + database['SMILES'].iloc[i[1][0][1]],
                        'Solvent': database['Solvent'].iloc[i[1][0][0]],
                        'Temperature, K': database['Temperature, K'].iloc[i[1][0][0]],
                        'Product energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][1]])),
                        'Reaction': 'electrochem.',
                        'Type': react_type,
                        'Mechanism': mech,
                        'Free energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][1]] - database['G, eV'].iloc[i[1][0][0]]))}
                dataframe.loc[len(dataframe)] = data

    dataframe = dataframe.drop_duplicates(subset = ['Reagent SMILES','Product SMILES', 'Free energy, eV'])
    dataframe = dataframe[dataframe['Free energy, eV'] < 0].reset_index(drop = True)
    return dataframe


def Thermochem_reactions(database):
    x = SMILES_analysis(database)
    database = database.drop(x[3]).reset_index(drop=True)
    all_comb = combinations(database.index, 2)
    all_comb = [i for i in all_comb]
    for i in enumerate(database.index):
        all_comb.append((i[0], i[0]))

    all_comb = [(i[0], i[1]) for i in all_comb if (database['Functional'].iloc[i[0]] == database['Functional'].iloc[i[1]] and database['Basis'].iloc[i[0]] == database['Basis'].iloc[i[1]] and database['Solvent'].iloc[i[0]] == database['Solvent'].iloc[i[1]])]

    all_comb = [(i[0], i[1]) for i in all_comb if float('{:.3f}'.format(database['Molar mass'].iloc[i[0]] + database['Molar mass'].iloc[i[1]])) in database['Molar mass'].values]

    SMILES = pandas.Series([(database['SMILES'].iloc[i[0]] + '.' +  database['SMILES'].iloc[i[1]]) for i in all_comb])

    molar_mass = [(database['Molar mass'].iloc[i[0]] + database['Molar mass'].iloc[i[1]]) for i in all_comb]

    functional = [(database['Functional'].iloc[i[0]]) for i in all_comb]

    basis = [(database['Basis'].iloc[i[0]]) for i in all_comb]

    spin_1 = numpy.array([(database['Multiplicity'].iloc[i[0]] - 1)/2 for i in all_comb])

    spin_2 = numpy.array([(database['Multiplicity'].iloc[i[1]] - 1)/2 for i in all_comb])

    spin = numpy.array([(spin_1[i[0]] + spin_2[i[0]]) if (spin_1[i[0]] + spin_2[i[0]]) != 1 else (spin_1[i[0]] - spin_2[i[0]]) for i in enumerate(spin_1)])

    mult = (2*(spin)+ 1)


    charge = [(database['Charge'].iloc[i[0]] + database['Charge'].iloc[i[1]]) for i in all_comb]


    temp = [(database['Temperature, K'].iloc[i[0]]) for i in all_comb]

    solvent = [(database['Solvent'].iloc[i[0]]) for i in all_comb]

    ZPE = [(database['ZPE, eV'].iloc[i[0]] + database['ZPE, eV'].iloc[i[1]]) for i in all_comb]

    H = [(database['H, eV'].iloc[i[0]] + database['H, eV'].iloc[i[1]]) for i in all_comb]

    TS = [(database['TS, eV'].iloc[i[0]] + database['TS, eV'].iloc[i[1]]) for i in all_comb]

    G = [(database['G, eV'].iloc[i[0]] + database['G, eV'].iloc[i[1]]) for i in all_comb]

    G_ZPE = [(database['G-ZPE, eV'].iloc[i[0]] + database['G-ZPE, eV'].iloc[i[1]]) for i in all_comb]

    dip = ['NaN' for i in enumerate(G_ZPE)]
    mol_vol = ['NaN' for i in enumerate(G_ZPE)]
    path = ['NaN' for i in enumerate(G_ZPE)]

    columns = ['SMILES', 'Molar mass', 'Functional', 'Basis', 'Multiplicity', 'Charge', 'Temperature, K', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'Molecular volume', 'path']

    df = pandas.DataFrame(numpy.array([SMILES, molar_mass, functional, basis, mult, charge, temp, solvent, ZPE, H, TS, G, G_ZPE, dip, mol_vol, path]).T, columns = columns)
    if df.empty == False:
        #df = df.loc[(df['Charge'] == -1) & (df['Multiplicity'] == 2)]
        database = pandas.concat([database, df], ignore_index = True)

    columns = ['Reagent SMILES', 'Reagent multiplicity', 'Reagent charge', 'Reagent energy, eV', 'Product SMILES','Product multiplicity', 'Product charge', 'Product energy, eV', 'Functional', 'Basis', 'Reaction SMILES', 'Solvent', 'Temperature, K', 'Reaction', 'Type', 'Mechanism', 'Free energy, eV']
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
                        'Functional': database['Functional'].iloc[i[1][0][0]],
                        'Basis': database['Basis'].iloc[i[1][0][0]],
                        'Reaction SMILES': database['SMILES'].iloc[i[1][0][0]] + '>>' + database['SMILES'].iloc[i[1][0][1]],
                        'Solvent': database['Solvent'].iloc[i[1][0][0]],
                        'Temperature, K': database['Temperature, K'].iloc[i[1][0][0]],
                        'Reaction': 'thermochem.',
                        'Type': react_type,
                        'Mechanism': mech,
                        'Free energy, eV': float('{:.2f}'.format(database['G, eV'].iloc[i[1][0][1]] - database['G, eV'].iloc[i[1][0][0]]))}
                dataframe.loc[len(dataframe)] = data

    dataframe = dataframe.drop_duplicates(subset = ['Free energy, eV'], keep = 'first')
    dataframe = dataframe[dataframe['Free energy, eV'] < 0].reset_index(drop = True)
    return dataframe

def Reorganization_energy(path1):
    database1 = pandas.read_csv(path1)
    path2 = input('Input path to your database:\t')
    database = pandas.read_csv(path2)
    x = SMILES_analysis(database)
    database2 = database.iloc[x[3]].reset_index(drop = True)
    all_comb = [i for i in combinations(database2.index, 2)]
    for i in enumerate(database2.index):
        all_comb.append((i[1], i[1]))
    all_comb = [(i[0], i[1]) for i in all_comb if (database2['Functional'].iloc[i[0]] == database2['Functional'].iloc[i[1]] and database2['Basis'].iloc[i[0]] == database2['Basis'].iloc[i[1]] and database2['Solvent'].iloc[i[0]] == database2['Solvent'].iloc[i[1]])]
    SMILES = pandas.Series([(database2['SMILES'].iloc[i[0]] + '.' +  database2['SMILES'].iloc[i[1]]) for i in all_comb])
    molar_mass = [(database2['Molar mass'].iloc[i[0]] + database2['Molar mass'].iloc[i[1]]) for i in all_comb]
    functional = [(database2['Functional'].iloc[i[0]]) for i in all_comb]
    basis = [(database2['Basis'].iloc[i[0]]) for i in all_comb]
    spin_1 = numpy.array([(database2['Multiplicity'].iloc[i[0]] - 1)/2 for i in all_comb])
    spin_2 = numpy.array([(database2['Multiplicity'].iloc[i[1]] - 1)/2 for i in all_comb])
    spin = numpy.array([(spin_1[i[0]] + spin_2[i[0]]) if (spin_1[i[0]] + spin_2[i[0]]) != 1 else (spin_1[i[0]] - spin_2[i[0]]) for i in enumerate(spin_1)])
    mult = (2*(spin)+ 1)
    charge = [(database2['Charge'].iloc[i[0]] + database2['Charge'].iloc[i[1]]) for i in all_comb]
    temp = [(database2['Temperature, K'].iloc[i[0]]) for i in all_comb]
    solvent = [(database2['Solvent'].iloc[i[0]]) for i in all_comb]
    ZPE = [(database2['ZPE, eV'].iloc[i[0]] + database2['ZPE, eV'].iloc[i[1]]) for i in all_comb]
    H = [(database2['H, eV'].iloc[i[0]] + database2['H, eV'].iloc[i[1]]) for i in all_comb]
    TS = [(database2['TS, eV'].iloc[i[0]] + database2['TS, eV'].iloc[i[1]]) for i in all_comb]
    G = [(database2['G, eV'].iloc[i[0]] + database2['G, eV'].iloc[i[1]]) for i in all_comb]
    G_ZPE = [(database2['G-ZPE, eV'].iloc[i[0]] + database2['G-ZPE, eV'].iloc[i[1]]) for i in all_comb]
    dip = ['NaN' for i in enumerate(G_ZPE)]
    mol_vol = ['NaN' for i in enumerate(G_ZPE)]
    path = ['NaN' for i in enumerate(G_ZPE)]
    columns = ['SMILES', 'Molar mass', 'Functional', 'Basis', 'Multiplicity', 'Charge', 'Temperature, K', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'Molecular volume', 'path']
    df = pandas.DataFrame(numpy.array([SMILES, molar_mass, functional, basis, mult, charge, temp, solvent, ZPE, H, TS, G, G_ZPE, dip, mol_vol, path]).T, columns = columns)

    #Outer sphere reorganization
    solv_reorg = []
    inner_reorg = []
    tot = []

    for i in enumerate(database1.index):
        SMILES_reagent = str(database1['Reagent SMILES'].iloc[i[0]])
        Charge_reagent = int(database1['Reagent charge'].iloc[i[0]])

        SMILES_product = str(database1['Product SMILES'].iloc[i[0]])
        Charge_product = int(database1['Product charge'].iloc[i[0]])

        if database1['Reduction potential, V'].iloc[i[0]] < 0:
            mol_vol = database['Molecular volume'].loc[(database['SMILES'] == SMILES_reagent) & (database['Charge'] == Charge_reagent)].values[0]
            mol_rad = (3*mol_vol/(4*numpy.pi))**(1/3)
            solvent_reorg_eng = (electron_charge**2)/(8*numpy.pi*vacuum_dielectric_permitivity)*(1/(solvent_refractive_index[database1['Solvent'].iloc[i[0]]]**2)-1/solvent_dielectric_permitivity[database1['Solvent'].iloc[i[0]]])*(1/mol_rad-1/(2*mol_rad))
            solvent_reorg_eng = float('{:.2f}'.format(solvent_reorg_eng*J_to_eV))
            solv_reorg.append(solvent_reorg_eng)
            mol_mass = database['Molar mass'].loc[database['SMILES'] == SMILES_reagent].values[0]
            if mol_mass in database2['Molar mass'].values and database1['Type'].iloc[i[0]] == 'adiab.':
                inner = database2['G, eV'].loc[(database2['Molar mass'] == mol_mass)].values[0] - database1['Reagent energy, eV'].iloc[i[0]]
                inner_reorg.append(inner)
                tot.append(inner + solvent_reorg_eng)
            elif mol_mass in df['Molar mass'].values and database1['Type'].iloc[i[0]] == 'dissoc.':
                inner = df['G, eV'].loc[(df['Molar mass'] == mol_mass)].values[0] - database1['Reagent energy, eV'].iloc[i[0]]
                inner_reorg.append(inner)
                tot.append(inner + solvent_reorg_eng)

        elif database1['Reduction potential, V'].iloc[i[0]] > 0:
            mol_vol = database['Molecular volume'].loc[(database['SMILES'] == SMILES_product) & (database['Charge'] == Charge_product)].values[0]
            mol_rad = (3*mol_vol/(4*numpy.pi))**(1/3)
            solvent_reorg_eng = (electron_charge**2)/(8*numpy.pi*vacuum_dielectric_permitivity)*(1/(solvent_refractive_index[database1['Solvent'].iloc[i[0]]]**2)-1/solvent_dielectric_permitivity[database1['Solvent'].iloc[i[0]]])*(1/mol_rad-1/(2*mol_rad))
            solvent_reorg_eng = float('{:.2f}'.format(solvent_reorg_eng*J_to_eV))
            solv_reorg.append(solvent_reorg_eng)
            mol_mass = database['Molar mass'].loc[database['SMILES'] == SMILES_reagent].values[0]
            if mol_mass in database2['Molar mass'].values and database1['Type'].iloc[i[0]] == 'adiab.':
                inner = database2['G, eV'].loc[(database2['Molar mass'] == mol_mass)].values[0] - database1['Product energy, eV'].iloc[i[0]]
                inner_reorg.append(inner)
                tot.append(inner + solvent_reorg_eng)
            elif mol_mass in df['Molar mass'].values and database1['Type'].iloc[i[0]] == 'dissoc.':
                inner = df['G, eV'].loc[(df['Molar mass'] == mol_mass)].values[0] - database1['Product energy, eV'].iloc[i[0]]
                inner_reorg.append(inner)
                tot.append(inner + solvent_reorg_eng)
        else:
            solv_reorg.append('NaN')
            inner_reorg.append('NaN')
            tot.append('NaN')

    database1['lambda (solv.), eV'] = solv_reorg
    database1['lambda (inner.), eV'] = inner_reorg
    database1['lambda (tot.), eV'] = tot
    database1 = database1.drop(columns = ['Unnamed: 0'])
    database1 = database1.reset_index(drop = True)
    return database1.to_csv(path1)

def Standard_rate(path):
    database = pandas.read_csv(path)
    path1 = input('Input path to your database:\t')
    database1 = pandas.read_csv(path1)
    rate_constant = []
    for i in enumerate(database.index):
        SMILES_reagent = str(database['Reagent SMILES'].iloc[i[0]])
        if database['Reaction'].iloc[i[0]] == 'electrochem.':
            mol_mass = database1['Molar mass'].loc[database1['SMILES'] == SMILES_reagent].values[0]
            Z = numpy.sqrt(gas_constant*database['Temperature, K'].iloc[i[0]]/2/numpy.pi/(mol_mass/1000))*100
            k0 = Z*numpy.exp(-(database['lambda (tot.), eV'].iloc[i[0]]/4)*faraday/gas_constant/database['Temperature, K'].iloc[i[0]])
            rate_constant.append(k0)
        if database['Reaction'].iloc[i[0]] == 'thermochem.':
            if '.' in SMILES_reagent:
                SMILES = SMILES_reagent.split('.')
                D = numpy.sum(numpy.array([database1['Diffusion, cm2/s'].loc[database1['SMILES'] == i].values[0] for i in SMILES])/10000)
                l = numpy.sum(numpy.array([(3*(database1['Molecular volume'].loc[(database1['SMILES'] == i)].values[0])/(4*numpy.pi))**(1/3) for i in SMILES]))
                k0 = 4*numpy.pi*Avogadro*D*l
                rate_constant.append(k0)
            else:
                k0 = k_Boltzman*database['Temperature, K'].iloc[i[0]]/Planck/10000
                rate_constant.append(k0)
    database['k0, cm/s'] = rate_constant
    database = database.drop(columns = ['Unnamed: 0'])
    database = database.reset_index(drop = True)
    return database.to_csv(path)

def Stability(path):
    database = pandas.read_csv(path)
    stability = []
    for i1 in enumerate(database.index):
        file = ORCA_parser(database['path'].iloc[i1[0]])
        stability.append(file.stability())
    database['Stability'] = stability
    database = database.drop(columns = ['Unnamed: 0'])
    return database.to_csv(path)


while True:
    print('''1. Generate reactions.
2. Choose reaction from your database.
3. Delete reactions from your database.
4. Calulate reduction potential.
5. Calculate reorganization energy.
6. Calaulate standard rate constant.
7. Stability analysis.
0. Exit to the main menu''')
    enter = int(input("Input number:\t"))
    if enter == 2:
        path = input('Input path to your database:\t')
        Reaction = Choose_reactions(path)
    elif enter == 1:
        path = input('Input path to your database:\t')
        Newdatabase = pandas.read_csv(path)
        Echem = input("Do you have electhochemical reactions?: (Y/N)\t")
        if Echem == 'Y':
            Electrochem = Electrochem_reactions(Newdatabase)
        elif Echem == 'N':
            Electrochem = pandas.DataFrame()

        Tchem = input("Do you have thermochemical reactions?: (Y/N)\t")
        if Tchem == 'Y':
            Reactions = Thermochem_reactions(Newdatabase)
        elif Tchem == 'N':
            Reactions = pandas.DataFrame()

        if Electrochem.empty == False and Reactions.empty == False:
            Database = pandas.concat([Electrochem, Reactions], ignore_index = True)
        elif Electrochem.empty == True:
            Database = Reactions
        elif Reactions.empty == True:
            Database = Electrochem
        Database.to_csv(input('Output path:\t') + '/' + 'Reaction_data.csv')
    elif enter == 3:
        path = input('Input path to your database:\t')
        Reaction = Drop_reactions(path)
    elif enter == 4:
        path = input('Input path to your database:\t')
        Reduction_potential(path)
    elif enter == 5:
        path = input('Input path to your reaction database:\t')
        Reorganization_energy(path)
    elif enter == 6:
        path = input('Input path to your reaction database:\t')
        #electrode = input('Electrode:\t')
        Standard_rate(path)
    elif enter ==7:
        path = input('Input path to your database:\t')
        Stability(path)
    else:
        print("Back to main menu")
        break
