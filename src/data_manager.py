#imports
import sys
sys.path.append('src')

import pandas
import numpy

def Choose_reactions(path):
    database = pandas.read_csv(path)
    print('Input index number (1 2 3.. n..)')
    Reaction_indexes = [int(i) for i in input("Enter the reaction index numbers:\t").split()]
    database = database.iloc[Reaction_indexes].reset_index()
    database = database.drop(columns = ['Unnamed: 0'])
    return database.to_csv(path)

def Drop_reactions(path):
    database = pandas.read_csv(path)
    print('Input index number (1 2 3.. n..)')
    Reaction_indexes = [int(i) for i in input("Enter the reaction index numbers:\t").split()]
    database = database.drop([Reaction_indexes]).reset_index()
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
    database1 = database1.reset_index()
    database1 = database1.drop(columns = ['Unnamed: 0'])
    return database1.to_csv(path)


while True:
    print('''1. Choose reaction from your database.
2. Delete reactions from your database
3. Calulate reduction potential
0. Exit to the main menu''')
    enter = int(input("Input number:\t"))
    if enter == 1:
        path = input('Input path to your database:\t')
        Reaction = Choose_reactions(path)
    elif enter == 2:
        path = input('Input path to your database:\t')
        Reaction = Drop_reactions(path)
    elif enter == 3:
        path = input('Input path to your database:\t')
        Reduction_potential(path)
    else:
        print("Back to main menu")
        break
