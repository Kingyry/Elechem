#imports
import pandas
import os
from parsers import ORCA_parser

reference_electrode = """Reference electrodes:
    1. Ferrocene
"""
print(reference_electrode)

if int(input("Enter the reference electrode number:\t")) == 1:
    if os.path.isfile('Data/Reference_electrodes/Ferrocene/Energy.csv') is True and os.path.isfile('Data/Reference_electrodes/Ferrocene/Ionization_potential.csv') is True:
        if str(input("Database files exist. Do you want to update? (Y/N)\t")) == "Y":
                input_path = input("Input path to the Database directory:\t")
                Ferrocene_folder = os.listdir(input_path)
                columns = ['Name', 'SMILES', 'Functional', 'Basis', 'Type', 'Charge', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'Dipole, D', 'path']
                dataframe=pandas.DataFrame(columns=columns)
                for item in Ferrocene_folder:
                    if item == 'Ferrocene':
                        path = input_path + '/' + Ferrocene_folder[item == 'Ferrocene']
                        Solvent_folder = os.listdir(path)
                        for item1 in range(len(Solvent_folder)):
                            path1 = path + '/' + Solvent_folder[item1]
                            Functional_folder = os.listdir(path1)
                            for item2 in range(len(Functional_folder)):
                                path2 = path1 + '/' + Functional_folder[item2]
                                Basis_folder = os.listdir(path2)
                                for item3 in range(len(Basis_folder)):
                                    path3 = path2 + '/' + Basis_folder[item3]
                                    File_folder = os.listdir(path3)
                                    for item4 in range(len(File_folder)):
                                        path4 = path3 + '/' + File_folder[item4]
                                        Files = os.listdir(path4)
                                        for item5 in range(len(Files)):
                                            if Files[item5].endswith('.txt'):
                                                path5 = path4 + '/' + Files[item5]
                                                file = ORCA_parser(path5)
                                            if file.thermodynamics()['G, eV'] != 'NaN':
                                                par = file.thermodynamics()['G, eV'] - file.thermodynamics()['ZPE, eV']
                                            else:
                                                par = 'NaN'
                                            if file.info()['charge'] == 0:
                                                name = 'Ferrocene'
                                                ion = 'Neutral'
                                            elif file.info()['charge'] == 1:
                                                    name = 'Ferrocenium'
                                                    ion = 'Cation'
                                            data={'Name': name,
                                                  'SMILES': file.smiles(),
                                                  'Functional': file.info()['functional'],
                                                  'Basis': file.info()['basis'],
                                                  'Type': ion,
                                                  'Charge': file.info()['charge'],
                                                  'Solvent': file.info()['solvation'],
                                                  'ZPE, eV': file.thermodynamics()['ZPE, eV'],
                                                  'H, eV': file.thermodynamics()['H, eV'],
                                                  'TS, eV': file.thermodynamics()['TS, eV'],
                                                  'G, eV': file.thermodynamics()['G, eV'],
                                                  'G-ZPE, eV': par,
                                                  'Dipole, D': file.info()['dipole'],
                                                  'path': path4}
                                            dataframe.loc[len(dataframe)]=data
                for i in range(len(dataframe['SMILES'])):
                            if dataframe['G, eV'].iloc[i] == 'NaN':
                                dataframe['G, eV'].iloc[i] = dataframe['ZPE, eV'].iloc[i] + dataframe[(dataframe['SMILES'] == dataframe['SMILES'].iloc[i]) & (dataframe['Functional'] == 'B3LYP') & (dataframe['Basis'] == 'DEF2-TZVPP') & (dataframe['Charge'] == dataframe['Charge'].iloc[i])]['G-ZPE, eV']
                dataframe.to_csv('Data/Reference_electrodes/Ferrocene/Energy.csv')

                dataframe = pandas.read_csv('Data/Reference_electrodes/Ferrocene/Energy.csv')
                functionals = []
                for i in dataframe['Functional']:
                    if i not in functionals:
                        functionals.append(i)
                basis = []
                for i in dataframe['Basis']:
                    if i not in basis:
                        basis.append(i)
                solvent = []
                for i in dataframe['Solvent']:
                    if i not in solvent:
                        solvent.append(i)
                ferrocene = []
                ferrocenium = []
                functional_dataframe = []
                basis_dataframe = []
                solvent_dataframe = []
                for item in range(len(dataframe['Name'])):
                    for item1 in functionals:
                        for item2 in basis:
                            for item3 in solvent:
                                if (dataframe['Name'][item] == 'Ferrocene' and dataframe['Functional'][item] == item1 and dataframe['Basis'][item] == item2 and dataframe['Solvent'][item] == item3):
                                    ferrocene.append(dataframe['G, eV'][item])
                                elif(dataframe['Name'][item] == 'Ferrocenium' and dataframe['Functional'][item] == item1 and dataframe['Basis'][item] == item2 and dataframe['Solvent'][item] == item3):
                                    ferrocenium.append(dataframe['G, eV'][item])
                                    functional_dataframe.append(dataframe['Functional'][item])
                                    basis_dataframe.append(dataframe['Basis'][item])
                                    solvent_dataframe.append(dataframe['Solvent'][item])
                Ionization_energy = [(ferrocenium[i] - ferrocene[i]) for i in range(len(ferrocene))]
                data = {'Functional': functional_dataframe,
                          'Basis': basis_dataframe,
                          'Solvent': solvent_dataframe,
                          'Ionization energy, eV': Ionization_energy}
                dataframe=pandas.DataFrame(data)
                dataframe.to_csv('Data/Reference_electrodes/Ferrocene/Ionization_energy.csv')

        else:
            exit
    else:
        if os.path.isfile('Data/Reference_electrodes/Ferrocene/Energy.csv') is False and os.path.isfile('Data/Reference_electrodes/Ferrocene/Ionization_potential.csv') is False:
            if str(input("Database file is abcent. Do you want to create? (Y/N)\t")) == "Y":
                input_path = input("Input path to the Database directory:\t")
                Ferrocene_folder = os.listdir(input_path)
                columns=['Name', 'SMILES', 'Functional', 'Basis', 'Type', 'Charge', 'Solvent', 'ZPE, eV', 'H, eV', 'TS, eV', 'G, eV', 'G-ZPE, eV', 'Dipole, D', 'path']
                dataframe=pandas.DataFrame(columns=columns)
                for item in Ferrocene_folder:
                    if item == 'Ferrocene':
                        path = input_path + '/' + Ferrocene_folder[item == 'Ferrocene']
                        Solvent_folder = os.listdir(path)
                        for item1 in range(len(Solvent_folder)):
                            path1 = path + '/' + Solvent_folder[item1]
                            Functional_folder = os.listdir(path1)
                            for item2 in range(len(Functional_folder)):
                                path2 = path1 + '/' + Functional_folder[item2]
                                Basis_folder = os.listdir(path2)
                                for item3 in range(len(Basis_folder)):
                                    path3 = path2 + '/' + Basis_folder[item3]
                                    File_folder = os.listdir(path3)
                                    for item4 in range(len(File_folder)):
                                        path4 = path3 + '/' + File_folder[item4]
                                        Files = os.listdir(path4)
                                        for item5 in range(len(Files)):
                                            if Files[item5].endswith('.txt'):
                                                path5 = path4 + '/' + Files[item5]
                                                file = ORCA_parser(path5)
                                            if file.thermodynamics()['G, eV'] != 'NaN':
                                                par = file.thermodynamics()['G, eV'] - file.thermodynamics()['ZPE, eV']
                                            else:
                                                par = 'NaN'
                                            if file.info()['charge'] == 0:
                                                name = 'Ferrocene'
                                                ion = 'Neutral'
                                            elif file.info()['charge'] == 1:
                                                    name = 'Ferrocenium'
                                                    ion = 'Cation'
                                            data={'Name': name,
                                                  'SMILES': file.smiles(),
                                                  'Functional': file.info()['functional'],
                                                  'Basis': file.info()['basis'],
                                                  'Type': ion,
                                                  'Charge': file.info()['charge'],
                                                  'Solvent': file.info()['solvation'],
                                                  'ZPE, eV': file.thermodynamics()['ZPE, eV'],
                                                  'H, eV': file.thermodynamics()['H, eV'],
                                                  'TS, eV': file.thermodynamics()['TS, eV'],
                                                  'G, eV': file.thermodynamics()['G, eV'],
                                                  'G-ZPE, eV': par,
                                                  'Dipole, D': file.info()['dipole'],
                                                  'path': path4}
                                            dataframe.loc[len(dataframe)]=data
                for i in range(len(dataframe['SMILES'])):
                            if dataframe['G, eV'].iloc[i] == 'NaN':
                                dataframe['G, eV'].iloc[i] = dataframe['ZPE, eV'].iloc[i] + dataframe[(dataframe['SMILES'] == dataframe['SMILES'].iloc[i]) & (dataframe['Functional'] == 'B3LYP') & (dataframe['Basis'] == 'DEF2-TZVPP') & (dataframe['Charge'] == dataframe['Charge'].iloc[i])]['G-ZPE, eV']
                dataframe.to_csv('Data/Reference_electrodes/Ferrocene/Energy.csv')


                dataframe = pandas.read_csv('Data/Reference_electrodes/Ferrocene/Energy.csv')
                functionals = []
                for i in dataframe['Functional']:
                    if i not in functionals:
                        functionals.append(i)
                basis = []
                for i in dataframe['Basis']:
                    if i not in basis:
                        basis.append(i)
                solvent = []
                for i in dataframe['Solvent']:
                    if i not in solvent:
                        solvent.append(i)
                ferrocene = []
                ferrocenium = []
                functional_dataframe = []
                basis_dataframe = []
                solvent_dataframe = []
                for item in range(len(dataframe['Name'])):
                    for item1 in functionals:
                        for item2 in basis:
                            for item3 in solvent:
                                if (dataframe['Name'][item] == 'Ferrocene' and dataframe['Functional'][item] == item1 and dataframe['Basis'][item] == item2 and dataframe['Solvent'][item] == item3):
                                    ferrocene.append(dataframe['G, eV'][item])
                                elif(dataframe['Name'][item] == 'Ferrocenium' and dataframe['Functional'][item] == item1 and dataframe['Basis'][item] == item2 and dataframe['Solvent'][item] == item3):
                                    ferrocenium.append(dataframe['G, eV'][item])
                                    functional_dataframe.append(dataframe['Functional'][item])
                                    basis_dataframe.append(dataframe['Basis'][item])
                                    solvent_dataframe.append(dataframe['Solvent'][item])
                Ionization_energy = [(ferrocenium[i] - ferrocene[i]) for i in range(len(ferrocene))]
                data = {'Functional': functional_dataframe,
                          'Basis': basis_dataframe,
                          'Solvent': solvent_dataframe,
                          'Ionization energy, eV': Ionization_energy}
                dataframe=pandas.DataFrame(data)
                dataframe.to_csv('Data/Reference_electrodes/Ferrocene/Ionization_energy.csv')
            else:
                exit
