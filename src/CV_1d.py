#imports
import sys
sys.path.append('src')

import numpy
import pandas
import os.path
from data import faraday, gas_constant

#Set environment parameters
#import of initial parameters

def Diffusion_constant(Reaction_database, Database):
    '''Function calculates model diffusion parameters for grig construction'''
    elements = []
    for i in enumerate(Reaction_database.index):
        if ('.' not in Reaction_database['Reagent SMILES'].iloc[i[0]]) and (Reaction_database['Reagent SMILES'].iloc[i[0]] not in elements):
            elements.append(Reaction_database['Reagent SMILES'].iloc[i[0]])
        elif ('.' in Reaction_database['Reagent SMILES'].iloc[i[0]]):
            for j in Reaction_database['Reagent SMILES'].iloc[i[0]].split('.'):
                if j not in elements:
                    elements.append(j)
        if ('.' not in Reaction_database['Product SMILES'].iloc[i[0]]) and (Reaction_database['Product SMILES'].iloc[i[0]] not in elements):
            elements.append(Reaction_database['Product SMILES'].iloc[i[0]])
        elif ('.' in Reaction_database['Product SMILES'].iloc[i[0]]):
            for j in Reaction_database['Product SMILES'].iloc[i[0]].split('.'):
                if j not in elements:
                    elements.append(j)
    Diff_coeff = []
    for i in elements:
        Diff_coeff.append(Database['Diffusion, cm2/s'].loc[Database['SMILES'] == i].values[0])
    for i in numpy.linspace(0.5, 0, 1000):
        D = i*numpy.array(Diff_coeff)/Diff_coeff[0]
        if all(D < 0.5):
            break
    data = {}
    for iterations in enumerate(elements):
        data[iterations[1]] = [D[iterations[0]], Diff_coeff[iterations[0]]]
    return data

def SMH(Reactions, iterations):
    '''Function calculates kinetic parameters for Symmetric Marcus-Hush model'''
    reactions = {}
    for i in enumerate(Reactions.index):
        CHARGE = (Reactions['Reagent charge'].iloc[i[0]] + Reactions['Product charge'].iloc[i[0]])
        if Reactions['Reaction'].iloc[i[0]] == 'electrochem.' and Reactions['Reagent charge'].iloc[i[0]] > CHARGE:
            I = numpy.arange(0, iterations, 1)

            Potential_sweep = (Total_scan-Reactions['Reduction potential, V'].iloc[i[0]])*faraday/gas_constant/Reactions['Temperature, K'].iloc[i[0]]
            G = Total_scan-Reactions['Reduction potential, V'].iloc[i[0]]
            alpha_f = 0.5 + (Total_scan - Reactions['Reduction potential, V'].iloc[i[0]])/4/Reactions['lambda (tot.), eV'].iloc[i[0]]
            E_act_f = Reactions['lambda (tot.), eV'].iloc[i[0]]/4 + alpha_f*(Total_scan - Reactions['Reduction potential, V'].iloc[i[0]])
            kf = Reactions['k0, cm/s'].iloc[i[0]]*numpy.exp(-alpha_f*Potential_sweep)
            if '.' in Reactions['Product SMILES'].iloc[i[0]]:
                alpha_b = numpy.ones(iterations)*numpy.nan
                E_act_b = numpy.zeros(iterations)
                kb = numpy.zeros(iterations)
            else:
                alpha_b = (1 - alpha_f)
                E_act_b = Reactions['lambda (tot.), eV'].iloc[i[0]]/4 - alpha_b*(Total_scan -Reactions['Reduction potential, V'].iloc[i[0]])
                kb = Reactions['k0, cm/s'].iloc[i[0]]*numpy.exp(alpha_b*Potential_sweep)

        elif Reactions['Reaction'].iloc[i[0]] == 'electrochem.' and Reactions['Reagent charge'].iloc[i[0]] <= CHARGE:
            I = numpy.arange(0, iterations, 1)

            Potential_sweep = (Total_scan - Reactions['Reduction potential, V'].iloc[i[0]])*faraday/gas_constant/Reactions['Temperature, K'].iloc[i[0]]

            alpha_f = 0.5 + (Total_scan - Reactions['Reduction potential, V'].iloc[i[0]])/4/Reactions['lambda (tot.), eV'].iloc[i[0]]
            alpha_f = 1 - alpha_f
            E_act_f = Reactions['lambda (tot.), eV'].iloc[i[0]]/4 - alpha_f*(Total_scan-Reactions['Reduction potential, V'].iloc[i[0]])
            kf = Reactions['k0, cm/s'].iloc[i[0]]*numpy.exp(alpha_f*Potential_sweep)
            if '.' in Reactions['Reagent SMILES'].iloc[i[0]]:
                alpha_b = numpy.ones(iterations)*numpy.nan
                E_act_b = numpy.zeros(iterations)
                kb = numpy.zeros(iterations)
            else:
                alpha_b = 0.5 + (Total_scan - Reactions['Reduction potential, V'].iloc[i[0]])/4/Reactions['lambda (tot.), eV'].iloc[i[0]]
                E_act_b = Reactions['lambda (tot.), eV'].iloc[i[0]]/4 + alpha_b*(Total_scan-Reactions['Reduction potential, V'].iloc[i[0]])
                kb = Reactions['k0, cm/s'].iloc[i[0]]*numpy.exp(-alpha_b*Potential_sweep)

        data = {'Iterations': I,
                'Potential [V]': Total_scan,
                'dG [eV]': G,
                'alp_{f}': alpha_f,
                'alp_{b}': alpha_b,
                'E_{act. f} [eV]': E_act_f,
                'E_{act. b} [eV]': E_act_b,
                'Potential sweep': Potential_sweep,
                'k_{f}': kf,
                'k_{b}': kb}
        reactions[Reactions['Reaction SMILES'].iloc[i[0]]] = data
    return reactions


print("***Set reactions and environment conditions***\n\n")
path = input('Input path to your reaction database:\t')
Reaction_database = pandas.read_csv(path)

#Separation on Reduction/Oxidation reactions
oxidation_index = []
reduction_index = []
chemical_index = []
for reaction in enumerate(Reaction_database.index):
    CHARGE = (Reaction_database['Reagent charge'].iloc[reaction[0]] + Reaction_database['Product charge'].iloc[reaction[0]])
    if Reaction_database['Reagent charge'].iloc[reaction[0]] > CHARGE and Reaction_database['Reaction'].iloc[reaction[0]] == 'electrochem.':
        reduction_index.append(reaction[0])
    elif Reaction_database['Reagent charge'].iloc[reaction[0]] <= CHARGE and Reaction_database['Reaction'].iloc[reaction[0]] == 'electrochem.':
        oxidation_index.append(reaction[0])
    elif Reaction_database['Reagent charge'].iloc[reaction[0]] == Reaction_database['Product charge'].iloc[reaction[0]] and Reaction_database['Reaction'].iloc[reaction[0]] == 'thermochem.':
        chemical_index.append(reaction[0])

#Generate reduction/oxidation/chemical reaction databases
if len(reduction_index) != 0:
    Reduction_reactions = Reaction_database.iloc[reduction_index].sort_values(by = ['Reduction potential, V'], ascending=False).reset_index(drop = True)
else:
    print('No reduction reactions found')
    Reduction_reactions = None
if len(oxidation_index) != 0:
    Oxidation_reactions = Reaction_database.iloc[oxidation_index].sort_values(by = ['Reduction potential, V'], ascending=True).reset_index(drop = True)
else:
    print('No oxidation reactions found')
    Oxidation_reactions = None
if len(chemical_index) != 0:
    Chemical_reactions = Reaction_database.iloc[chemical_index]
else:
    print('No thermochemical reactions found')
    Chemical_reactions = None
if len(reduction_index) != 0 and len(oxidation_index) != 0:
    raise ValueError('ERROR: Choose one type of the process (reduction/oxidation)')

#Find Independent/dependent electrochemical reactions
#Electroreduction
independent_reaction_index = []
dependent_reaction_index = []
if reduction_index != 0:
    Reduction_reactions_reagents = []
    Reduction_reactions_products = []
    index = []
    for it in enumerate(Reduction_reactions.index):
        if '.' in Reduction_reactions['Reagent SMILES'].iloc[it[0]]:
            for it1 in Reduction_reactions['Reagent SMILES'].iloc[it[0]].split('.'):
                Reduction_reactions_reagents.append(it1)
                index.append(it[0])
        else:
            Reduction_reactions_reagents.append(Reduction_reactions['Reagent SMILES'].iloc[it[0]])
            index.append(it[0])

        if '.' in Reduction_reactions['Product SMILES'].iloc[it[0]]:
            for it1 in Reduction_reactions['Product SMILES'].iloc[it[0]].split('.'):
                Reduction_reactions_products.append(it1)
        else:
            Reduction_reactions_products.append(Reduction_reactions['Product SMILES'].iloc[it[0]])
    for it2 in enumerate(Reduction_reactions_reagents):
        if it2[1] not in Reduction_reactions_products:
            independent_reaction_index.append(index[it2[0]])
        elif index[it2[0]] not in dependent_reaction_index:
            dependent_reaction_index.append(it2[0])
    Independent_reactions = Reduction_reactions.iloc[independent_reaction_index]
    Dependent_reactions = Reduction_reactions.iloc[dependent_reaction_index]

#Electro-oxidation
elif oxidation_index != 0:
    Oxidation_reactions_reagents = []
    Oxidation_reactions_products = []
    index = []
    for it in enumerate(Oxidation_reactions.index):
        if '.' in Oxidation_reactions['Product SMILES'].iloc[it[0]]:
            for it1 in Oxidation_reactions['Product SMILES'].iloc[it[0]].split('.'):
                Oxidation_reactions_reagents.append(it1)
                index.append(it[0])
        else:
            Oxidation_reactions_reagents.append(Oxidation_reactions['Product SMILES'].iloc[it[0]])
            index.append(it[0])

        if '.' in Oxidation_reactions['Reagent SMILES'].iloc[it[0]]:
            for it1 in Oxidation_reactions['Reagent SMILES'].iloc[it[0]].split('.'):
                Oxidation_reactions_products.append(it1)
        else:
            Oxidation_reactions_products.append(Oxidation_reactions['Reagent SMILES'].iloc[it[0]])
    for it2 in enumerate(Oxidation_reactions_reagents):
        if it2[1] not in Oxidation_reactions_products:
            independent_reaction_index.append(index[it2[0]])
        elif index[it2[0]] not in dependent_reaction_index:
            dependent_reaction_index.append(it2[0])
    Independent_reactions = Oxidation_reactions.iloc[independent_reaction_index]
    Dependent_reactions = Oxidation_reactions.iloc[dependent_reaction_index]

path = input('\nInput path to your initial information database:\t')
Database = pandas.read_csv(path)
output_path = input('\nEnter output path:\t')
f = open(output_path + "/" + "OUTPUT.out", "w")

#Model variables
f.write("\t\tElectrochemistry application output file\n\t\t\tCyclic voltammetry module output")
Potential_max = float('{:.2f}'.format(float(input('Potential [V]: max\t'))))
Potential_min = float('{:.2f}'.format(float(input('Potential [V]: min\t'))))
f.write('\n\nPotential [V]:\tmax {0}\n\t\t\t\tmin {1}\n'.format(Potential_max, Potential_min))
Potential_scan = float(input('Potential scan rate [V/s]:\t'))
f.write('Potential scan rate [V/s]: {0}\n'.format(Potential_scan))
Time = abs(Potential_max - Potential_min)*2/Potential_scan
f.write('Simulation total time [s]: {0}\n'.format('{:.2f}'.format(Time)))
dT = 0.1
iterations = int(Time/dT)
f.write('Number of iterations: {0}\n'.format('{:.0f}'.format(iterations)))
boxes = numpy.max(numpy.array([int(6*numpy.sqrt(Diffusion_constant(Reaction_database, Database)[i][0]*iterations)+1) for i in Diffusion_constant(Reaction_database, Database)]))
f.write('Number of boxes: {0}\n'.format('{:.0f}'.format(boxes)))
dX = numpy.max(numpy.array([numpy.sqrt(Diffusion_constant(Reaction_database, Database)[i][1]*dT/Diffusion_constant(Reaction_database, Database)[i][0]) for i in Diffusion_constant(Reaction_database, Database)]))
X = dX*boxes
distance = [0]
for i in range(boxes):
    if i > 0:
        distance.append(distance[i-1] + dX)
f.write('Simulation region [cm]: {0}\n'.format('{:.2f}'.format(X)))

initial_concentrations = {}
for it in enumerate(Independent_reactions.index):
    CHARGE = (Independent_reactions['Reagent charge'].iloc[it[0]] + Independent_reactions['Product charge'].iloc[it[0]])
    if Independent_reactions['Reagent charge'].iloc[it[0]] > CHARGE:
        initial_conc = float(input('Enter initial concentration of {0} [mol/cm3]:\t'.format(Independent_reactions['Reagent SMILES'].iloc[it[0]])))
        f.write('{0} initial concentration [mol/cm3]: {1}\n'.format(Independent_reactions['Reagent SMILES'].iloc[it[0]], float('{:.2e}'.format(initial_conc))))
        initial_concentrations[Independent_reactions['Reagent SMILES'].iloc[it[0]]] = initial_conc
    elif Independent_reactions['Reagent charge'].iloc[it[0]] <= CHARGE:
        initial_conc = float(input('Enter initial concentration of {0} [mol/cm3]:\t'.format(Independent_reactions['Product SMILES'].iloc[it[0]])))
        f.write('{0} initial concentration [mol/cm3]: {1}\n'.format(Independent_reactions['Product SMILES'].iloc[it[0]], float('{:.2e}'.format(initial_conc))))
        initial_concentrations[Independent_reactions['Product SMILES'].iloc[it[0]]] = initial_conc

Forward_scan = numpy.linspace(Potential_max, Potential_min, int(iterations/2), endpoint=False)
Backward_scan = numpy.linspace(Potential_min, Potential_max, int(iterations/2))
Total_scan = numpy.concatenate((Forward_scan, Backward_scan))
EReact = {'Reactions': Reaction_database['Reaction SMILES'].values,
          'Free energy [eV]': Reaction_database['Free energy, eV'].values,
          'Reduction potential [V]': Reaction_database['Reduction potential, V'].values}
EReact = pandas.DataFrame(EReact)
f.write(EReact.to_string(header = True, index = True, justify= 'center'))

#Kinetic model and simulation calculator
kinetic_model = input('Enter electron transfer model name (BV/SMH/AMH):\t')


if kinetic_model in ['SMH', 'smh']:
    #add kinetics to output file
    f.write('\n\nKintetic model is Simmetric Marcus-Hush Model')
    for reactions in SMH(Independent_reactions, iterations):
        f.write('\n\nReaction {0}:\t'.format(reactions) + '\n')
        DATA = pandas.DataFrame(SMH(Independent_reactions, iterations)[reactions])
        f.write(DATA.to_string(header = True, index = False, justify= 'center'))
    for reactions in SMH(Dependent_reactions, iterations):
        f.write('\n\nReaction {0}:\t'.format(reactions) + '\n')
        DATA = pandas.DataFrame(SMH(Dependent_reactions, iterations)[reactions])
        f.write(DATA.to_string(header = True, index = False, justify= 'center'))

    #calculate reaction kinetics and system responces
    kinetics_independent = SMH(Independent_reactions, iterations)
    kinetics_dependent = SMH(Dependent_reactions, iterations)
else:
    raise ValueError('ERROR: Wrong model')

f.write('\n\nCurrent densities from electrochemical reactions J [mA/cm2]')

#Find all reagents for electrochemical and chemical reactions
Reagents = []
Products = []
for reactions in enumerate(Reaction_database.index):
    CHARGE = (Reaction_database['Reagent charge'].iloc[reactions[0]] + Reaction_database['Product charge'].iloc[reactions[0]])
    if Reaction_database['Reagent charge'].iloc[reactions[0]] > CHARGE and Reaction_database['Reaction'].iloc[reactions[0]] == 'electrochem.':
        if '.' in Reaction_database['Reagent SMILES'].iloc[reactions[0]]:
            for reagent in Reaction_database['Reagent SMILES'].iloc[reactions[0]].split('.'):
                if reagent not in Reagents:
                    Reagents.append(reagent)
        else:
            if Reaction_database['Reagent SMILES'].iloc[reactions[0]] not in Reagents:
                Reagents.append(Reaction_database['Reagent SMILES'].iloc[reactions[0]])
        if '.' in Reaction_database['Product SMILES'].iloc[reactions[0]]:
            for product in Reaction_database['Product SMILES'].iloc[reactions[0]].split('.'):
                if product not in Products:
                    Products.append(product)
        else:
            if Reaction_database['Product SMILES'].iloc[reactions[0]] not in Products:
                Products.append(Reaction_database['Product SMILES'].iloc[reactions[0]])
    elif Reaction_database['Reagent charge'].iloc[reactions[0]] <= CHARGE and Reaction_database['Reaction'].iloc[reactions[0]] == 'electrochem.':
        if '.' in Reaction_database['Product SMILES'].iloc[reactions[0]]:
            for reagent in Reaction_database['Product SMILES'].iloc[reactions[0]].split('.'):
                if reagent not in Reagents:
                    Reagents.append(reagent)
        else:
            if Reaction_database['Product SMILES'].iloc[reactions[0]] not in Reagents:
                Reagents.append(Reaction_database['Product SMILES'].iloc[reactions[0]])
        if '.' in Reaction_database['Reagent SMILES'].iloc[reactions[0]]:
            for product in Reaction_database['Reagent SMILES'].iloc[reactions[0]].split('.'):
                if product not in Products:
                    Products.append(product)
        else:
            if Reaction_database['Reagent SMILES'].iloc[reactions[0]] not in Products:
                Products.append(Reaction_database['Reagent SMILES'].iloc[reactions[0]])
    elif Reaction_database['Reagent charge'].iloc[reactions[0]] == Reaction_database['Product charge'].iloc[reactions[0]] and Reaction_database['Reaction'].iloc[reactions[0]] == 'thermochem.':
        if '.' in Reaction_database['Reagent SMILES'].iloc[reactions[0]]:
            for reagent in Reaction_database['Reagent SMILES'].iloc[reactions[0]].split('.'):
                if reagent not in Reagents:
                    Reagents.append(reagent)
        else:
            if Reaction_database['Reagent SMILES'].iloc[reactions[0]] not in Reagents:
                Reagents.append(Reaction_database['Reagent SMILES'].iloc[reactions[0]])
        if '.' in Reaction_database['Product SMILES'].iloc[reactions[0]]:
            for product in Reaction_database['Product SMILES'].iloc[reactions[0]].split('.'):
                if product not in Products:
                    Products.append(product)
        else:
            if Reaction_database['Product SMILES'].iloc[reactions[0]] not in Products:
                Products.append(Reaction_database['Product SMILES'].iloc[reactions[0]])

#Reagents concentration matrixes
Reagents_distribution_with_chem_reactions = {}
for matrix in Reagents:
    if matrix in initial_concentrations:
        Reagents_distribution_with_chem_reactions[matrix] = numpy.ones((iterations, boxes))*initial_concentrations[matrix]
    else:
        Reagents_distribution_with_chem_reactions[matrix] = numpy.zeros((iterations, boxes))
#Products concentration matrixes
Products_distribution_with_chem_reactions = {}
for matrix in Products:
    Products_distribution_with_chem_reactions[matrix] = numpy.zeros((iterations, boxes))

#Reagents fluxes
Reagents_fluxes_with_chem_reaction = {}
for reactions in enumerate(Reaction_database.index):
    CHARGE = (Reaction_database['Reagent charge'].iloc[reactions[0]] + Reaction_database['Product charge'].iloc[reactions[0]])
    if Reaction_database['Reagent charge'].iloc[reactions[0]] > CHARGE and Reaction_database['Reaction'].iloc[reactions[0]] == 'electrochem.':
        if '.' in Reaction_database['Reagent SMILES'].iloc[reactions[0]]:
            for reagent in Reaction_database['Reagent SMILES'].iloc[reactions[0]].split('.'):
                Reagents_fluxes_with_chem_reaction[reagent] = numpy.zeros(iterations)
        else:
            Reagents_fluxes_with_chem_reaction[Reaction_database['Reagent SMILES'].iloc[reactions[0]]] = numpy.zeros(iterations)
    elif Reaction_database['Reagent charge'].iloc[reactions[0]] <= CHARGE and Reaction_database['Reaction'].iloc[reactions[0]] == 'electrochem.':
        if '.' in Reaction_database['Product SMILES'].iloc[reactions[0]]:
            for reagent in Reaction_database['Product SMILES'].iloc[reactions[0]].split('.'):
                Reagents_fluxes_with_chem_reaction[reagent] = numpy.zeros(iterations)
        else:
            Reagents_fluxes_with_chem_reaction[Reaction_database['Product SMILES'].iloc[reactions[0]]] = numpy.zeros(iterations)
#Products fluxes
Products_fluxes_with_chem_reaction = {}
for reactions in enumerate(Reaction_database.index):
    CHARGE = (Reaction_database['Reagent charge'].iloc[reactions[0]] + Reaction_database['Product charge'].iloc[reactions[0]])
    if Reaction_database['Reagent charge'].iloc[reactions[0]] > CHARGE and Reaction_database['Reaction'].iloc[reactions[0]] == 'electrochem.':
        if '.' in Reaction_database['Product SMILES'].iloc[reactions[0]]:
            for product in Reaction_database['Product SMILES'].iloc[reactions[0]].split('.'):
                Products_fluxes_with_chem_reaction[product] = numpy.zeros(iterations)
        else:
            Products_fluxes_with_chem_reaction[Reaction_database['Product SMILES'].iloc[reactions[0]]] = numpy.zeros(iterations)
    elif Reaction_database['Product charge'].iloc[reactions[0]] <= CHARGE and Reaction_database['Reaction'].iloc[reactions[0]] == 'electrochem.':
        if '.' in Reaction_database['Reagent SMILES'].iloc[reactions[0]]:
            for product in Reaction_database['Reagent SMILES'].iloc[reactions[0]].split('.'):
                Products_fluxes_with_chem_reaction[product] = numpy.zeros(iterations)
        else:
            Products_fluxes_with_chem_reaction[Reaction_database['Reagent SMILES'].iloc[reactions[0]]] = numpy.zeros(iterations)

#Chemical reagents fluxes
Chemical_reagents_fluxes = {}
if Chemical_reactions != None:
    for reactions in enumerate(Chemical_reactions.index):
        Chemical_reagents_fluxes[Chemical_reactions['Reaction SMILES'].iloc[reactions[0]]] = numpy.zeros((iterations, boxes))
#Chemical products fluxes
    Chemical_products_fluxes = {}
    for reactions in enumerate(Chemical_reactions.index):
        Chemical_products_fluxes[Chemical_reactions['Reaction SMILES'].iloc[reactions[0]]] = numpy.zeros((iterations, boxes))

#Currents
Reaction_currents = {}
for reactions in enumerate(Reaction_database.index):
    if Reaction_database['Reaction'].iloc[reactions[0]] == 'electrochem.':
        Reaction_currents[Reaction_database['Reaction SMILES'].iloc[reactions[0]]] = numpy.zeros(iterations)

#Calculator electroreduction
for i1 in range(iterations):
    #All reactions iterator
    for reactions in enumerate(Reaction_database.index):
        #Independent reaction iterator
        for indpndt_reactions in enumerate(Independent_reactions.index):
            #compare Reaction SMILES in all reactions and independent reactions
            if Reaction_database['Reaction SMILES'].iloc[reactions[0]] == Independent_reactions['Reaction SMILES'].iloc[indpndt_reactions[0]]:
                #Find chemical reactions with the same reagent as electrochemical reaction
                CHARGE = (Independent_reactions['Reagent charge'].iloc[indpndt_reactions[0]] + Independent_reactions['Product charge'].iloc[indpndt_reactions[0]])
                #Find electrochemical reagents
                if Independent_reactions['Reagent charge'].iloc[indpndt_reactions[0]] > CHARGE:
                #Find electrochemical products electroreduction
                    if '.' in Independent_reactions['Reagent SMILES'].iloc[indpndt_reactions[0]]:
                        Electro_reagents = Independent_reactions['Reagent SMILES'].iloc[indpndt_reactions[0]].split('.')
                    else:
                        Electro_reagents = [Independent_reactions['Reagent SMILES'].iloc[indpndt_reactions[0]]]
                #Find electrochemical products electroreduction
                    if '.' in Independent_reactions['Product SMILES'].iloc[indpndt_reactions[0]]:
                        Electro_products = Independent_reactions['Product SMILES'].iloc[indpndt_reactions[0]].split('.')
                    else:
                        Electro_products = [Independent_reactions['Product SMILES'].iloc[indpndt_reactions[0]]]

                elif Independent_reactions['Reagent charge'].iloc[indpndt_reactions[0]] <= CHARGE:
                #Find electrochemical reagents electro-oxidation
                    if '.' in Independent_reactions['Product SMILES'].iloc[indpndt_reactions[0]]:
                        Electro_reagents = Independent_reactions['Product SMILES'].iloc[indpndt_reactions[0]].split('.')
                    else:
                        Electro_reagents = [Independent_reactions['Product SMILES'].iloc[indpndt_reactions[0]]]
                #Find electrochemical products electro-oxidation
                    if '.' in Independent_reactions['Reagent SMILES'].iloc[indpndt_reactions[0]]:
                        Electro_products = Independent_reactions['Reagent SMILES'].iloc[indpndt_reactions[0]].split('.')
                    else:
                        Electro_products = [Independent_reactions['Reagent SMILES'].iloc[indpndt_reactions[0]]]

                for i2 in range(1, boxes-1):
                    if i1 != 0:
                        for it3 in Electro_reagents:
                            Reagents_distribution_with_chem_reactions[it3][i1, i2] = Reagents_distribution_with_chem_reactions[it3][i1-1, i2] + Diffusion_constant(Reaction_database, Database)[it3][0]*(Reagents_distribution_with_chem_reactions[it3][i1-1, i2+1] - 2*Reagents_distribution_with_chem_reactions[it3][i1-1, i2] + Reagents_distribution_with_chem_reactions[it3][i1-1, i2-1])
                        for it4 in Electro_products:
                            Products_distribution_with_chem_reactions[it4][i1, i2] = Products_distribution_with_chem_reactions[it4][i1-1, i2] + Diffusion_constant(Reaction_database, Database)[it4][0]*(Products_distribution_with_chem_reactions[it4][i1-1, i2+1] - 2*Products_distribution_with_chem_reactions[it4][i1-1, i2] + Products_distribution_with_chem_reactions[it4][i1-1, i2-1])

                for it_rf in Electro_reagents:
                    for it_pf in Electro_products:
                        #fluxes with chemical reactions
                        Reagents_fluxes_with_chem_reaction[it_rf][i1] = -(kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[indpndt_reactions[0]]]['k_{f}'][i1]*Reagents_distribution_with_chem_reactions[it_rf][i1, 1] - kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[indpndt_reactions[0]]]['k_{b}'][i1]*Products_distribution_with_chem_reactions[it_pf][i1, 1])/(1 + kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[indpndt_reactions[0]]]['k_{f}'][i1]*dX/Diffusion_constant(Reaction_database, Database)[it_rf][1] + kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[indpndt_reactions[0]]]['k_{b}'][i1]*dX/Diffusion_constant(Reaction_database, Database)[it_pf][1])
                        Products_fluxes_with_chem_reaction[it_pf][i1] = -Reagents_fluxes_with_chem_reaction[it_rf][i1]

                        #surface concentration change
                        Reagents_distribution_with_chem_reactions[it_rf][i1, 0] = Reagents_distribution_with_chem_reactions[it_rf][i1, 1] + Reagents_fluxes_with_chem_reaction[it_rf][i1]*dX/Diffusion_constant(Reaction_database, Database)[it_rf][1]
                        Products_distribution_with_chem_reactions[it_pf][i1, 0] = Products_distribution_with_chem_reactions[it_pf][i1, 1] + Products_fluxes_with_chem_reaction[it_pf][i1]*dX/Diffusion_constant(Reaction_database, Database)[it_pf][1]

                        #Update Reagents dictionary
                        if it_pf in Reagents_distribution_with_chem_reactions:
                            Reagents_distribution_with_chem_reactions[it_pf][i1, :] = Products_distribution_with_chem_reactions[it_pf][i1, :]

                    #Current densities
                    Reaction_currents[Independent_reactions['Reaction SMILES'].iloc[indpndt_reactions[0]]][i1] = faraday*1000*Reagents_fluxes_with_chem_reaction[it_rf][i1]

    for reactions in enumerate(Reaction_database.index):
        #Dependent reactions iterator
        for dpndt_reactions in enumerate(Dependent_reactions.index):
            #compare Reaction SMILES in all reactions and dependent reactions
            if Reaction_database['Reaction SMILES'].iloc[reactions[0]] == Dependent_reactions['Reaction SMILES'].iloc[dpndt_reactions[0]]:
                #Find chemical reactions with the same reagent as electrochemical reaction
                CHARGE = (Dependent_reactions['Reagent charge'].iloc[dpndt_reactions[0]] + Dependent_reactions['Product charge'].iloc[dpndt_reactions[0]])
                #Find electrochemical reagents
                if Dependent_reactions['Reagent charge'].iloc[dpndt_reactions[0]] > CHARGE:
                #Find electrochemical products electroreduction
                    if '.' in Dependent_reactions['Reagent SMILES'].iloc[dpndt_reactions[0]]:
                        Electro_reagents = Dependent_reactions['Reagent SMILES'].iloc[dpndt_reactions[0]].split('.')
                    else:
                        Electro_reagents = [Dependent_reactions['Reagent SMILES'].iloc[dpndt_reactions[0]]]
                #Find electrochemical products electroreduction
                    if '.' in Dependent_reactions['Product SMILES'].iloc[dpndt_reactions[0]]:
                        Electro_products = Dependent_reactions['Product SMILES'].iloc[dpndt_reactions[0]].split('.')
                    else:
                        Electro_products = [Dependent_reactions['Product SMILES'].iloc[dpndt_reactions[0]]]

                elif Dependent_reactions['Reagent charge'].iloc[dpndt_reactions[0]] <= CHARGE:
                #Find electrochemical reagents electro-oxidation
                    if '.' in Dependent_reactions['Product SMILES'].iloc[dpndt_reactions[0]]:
                        Electro_reagents = Dependent_reactions['Product SMILES'].iloc[dpndt_reactions[0]].split('.')
                    else:
                        Electro_reagents = [Dependent_reactions['Product SMILES'].iloc[dpndt_reactions[0]]]
                #Find electrochemical products electro-oxidation
                    if '.' in Dependent_reactions['Reagent SMILES'].iloc[dpndt_reactions[0]]:
                        Electro_products = Dependent_reactions['Reagent SMILES'].iloc[dpndt_reactions[0]].split('.')
                    else:
                        Electro_products = [Dependent_reactions['Reagent SMILES'].iloc[dpndt_reactions[0]]]

                for i2 in range(1, boxes-1):
                    if i1 != 0:
                        for it3 in Electro_reagents:
                            Reagents_distribution_with_chem_reactions[it3][i1, i2] = Reagents_distribution_with_chem_reactions[it3][i1-1, i2] + Diffusion_constant(Reaction_database, Database)[it3][0]*(Reagents_distribution_with_chem_reactions[it3][i1-1, i2+1] - 2*Reagents_distribution_with_chem_reactions[it3][i1-1, i2] + Reagents_distribution_with_chem_reactions[it3][i1-1, i2-1])
                        for it4 in Electro_products:
                            Products_distribution_with_chem_reactions[it4][i1, i2] = Products_distribution_with_chem_reactions[it4][i1-1, i2] + Diffusion_constant(Reaction_database, Database)[it4][0]*(Products_distribution_with_chem_reactions[it4][i1-1, i2+1] - 2*Products_distribution_with_chem_reactions[it4][i1-1, i2] + Products_distribution_with_chem_reactions[it4][i1-1, i2-1])

                for it_rf in Electro_reagents:
                    for it_pf in Electro_products:
                        #fluxes with chemical reactions
                        Reagents_fluxes_with_chem_reaction[it_rf][i1] = -(kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[dpndt_reactions[0]]]['k_{f}'][i1]*Reagents_distribution_with_chem_reactions[it_rf][i1, 1] - kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[dpndt_reactions[0]]]['k_{b}'][i1]*Products_distribution_with_chem_reactions[it_pf][i1, 1])/(1 + kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[dpndt_reactions[0]]]['k_{f}'][i1]*dX/Diffusion_constant(Reaction_database, Database)[it_rf][1] + kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[dpndt_reactions[0]]]['k_{b}'][i1]*dX/Diffusion_constant(Reaction_database, Database)[it_pf][1])
                        Products_fluxes_with_chem_reaction[it_pf][i1] = -Reagents_fluxes_with_chem_reaction[it_rf][i1]

                        #surface concentration change
                        Reagents_distribution_with_chem_reactions[it_rf][i1, 0] = Reagents_distribution_with_chem_reactions[it_rf][i1, 1] + (Reagents_fluxes_with_chem_reaction[it_rf][i1] + Products_fluxes_with_chem_reaction[it_rf][i1])*dX/Diffusion_constant(Reaction_database, Database)[it_rf][1]
                        Products_distribution_with_chem_reactions[it_pf][i1, 0] = Products_distribution_with_chem_reactions[it_pf][i1, 1] + Products_fluxes_with_chem_reaction[it_pf][i1]*dX/Diffusion_constant(Reaction_database, Database)[it_pf][1]


                        #Update Reagents dictionary
                        if it_pf in Reagents_distribution_with_chem_reactions:
                            Reagents_distribution_with_chem_reactions[it_pf][i1, :] = Products_distribution_with_chem_reactions[it_pf][i1, :]

                    #Current densities
                    Reaction_currents[Dependent_reactions['Reaction SMILES'].iloc[dpndt_reactions[0]]][i1] = faraday*1000*Reagents_fluxes_with_chem_reaction[it_rf][i1]

    if Chemical_reactions != None:
        for reactions in enumerate(Reaction_database.index):
            #Chemical reactions iterator
            for chem_reactions in enumerate(Chemical_reactions.index):
                #compare Reaction SMILES in all reactions and chemical reactions
                if Reaction_database['Reaction SMILES'].iloc[reactions[0]] == Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]:
                    if '.' in Chemical_reactions['Reagent SMILES'].iloc[chem_reactions[0]]:
                        Chemical_reagents = Chemical_reactions['Reagent SMILES'].iloc[chem_reactions[0]].split('.')
                    else:
                        Chemical_reagents = [Chemical_reactions['Reagent SMILES'].iloc[chem_reactions[0]]]
                    #Find chemical products
                    if '.' in Chemical_reactions['Product SMILES'].iloc[chem_reactions[0]]:
                        Chemical_products = Chemical_reactions['Product SMILES'].iloc[chem_reactions[0]].split('.')
                    else:
                        Chemical_products = [Chemical_reactions['Product SMILES'].iloc[chem_reactions[0]]]

                    #Chemical reaction flux
                    if len(Chemical_reagents) != 1:
                        Chemical_reagents_fluxes[Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]][i1, :] = -Chemical_reactions['k0, cm/s'].iloc[chem_reactions[0]]*Reagents_distribution_with_chem_reactions[Chemical_reagents[0]][i1, :]*Reagents_distribution_with_chem_reactions[Chemical_reagents[1]][i1, :]*dT
                        Chemical_products_fluxes[Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]][i1, :] = -Chemical_reagents_fluxes[Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]][i1, :]
                    else:
                        Chemical_reagents_fluxes[Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]][i1, :] = -Chemical_reactions['k0, cm/s'].iloc[chem_reactions[0]]*Reagents_distribution_with_chem_reactions[Chemical_reagents[0]][i1, :]*dT
                        Chemical_products_fluxes[Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]][i1, :] = -Chemical_reagents_fluxes[i1, :]
                    #Update reagents concentration
                    for it3 in Chemical_reagents:
                        Reagents_distribution_with_chem_reactions[it3][i1, :] =  Reagents_distribution_with_chem_reactions[it3][i1, :] + Chemical_reagents_fluxes[Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]][i1, :]

                    for it4 in Chemical_products:
                        Products_distribution_with_chem_reactions[it4][i1, :] = Products_distribution_with_chem_reactions[it4][i1, :] + Chemical_products_fluxes[Chemical_reactions['Reaction SMILES'].iloc[chem_reactions[0]]][i1, :]

Final_distributions = {}
for i in Reagents_distribution_with_chem_reactions:
    if i not in Final_distributions:
        Final_distributions[i] = Reagents_distribution_with_chem_reactions[i]

for i in Products_distribution_with_chem_reactions:
    if i not in Final_distributions:
        Final_distributions[i] = Products_distribution_with_chem_reactions[i]

#Create reaction currents database
Currents = pandas.DataFrame(numpy.array([Reaction_currents[i] for i in Reaction_currents]).T, columns = [i for i in Reaction_currents])
Currents.rename(columns = {'Index': 'Iterations'})
Currents['Total current'] = Currents.sum(axis = 1)
f.write('\n' + Currents.to_string(header = True, index = True, justify= 'center'))
f.close()

conc = open(output_path + "/" + "CONCENTRATIONS.out", "w")
conc.write("Concentration distributions")
conc.write('\n\nPotential [V]:\tmax {0}\n\t\t\t\tmin {1}\n'.format(Potential_max, Potential_min))
conc.write('Potential scan rate [V/s]: {0}\n'.format(Potential_scan))
conc.write('Simulation total time [s]: {0}\n'.format('{:.2f}'.format(Time)))
conc.write('Number of iterations: {0}\n'.format('{:.0f}'.format(iterations)))
conc.write('Number of boxes: {0}\n'.format('{:.0f}'.format(boxes)))
conc.write('Simulation region [cm]: {0}\n'.format('{:.2e}'.format(X)))
conc.write('Simulation region dX [cm]: {0}\n'.format('{:.2e}'.format(dX)))

for chemical_compound in Final_distributions:
    distribution = pandas.DataFrame(Final_distributions[chemical_compound], columns= distance)
    conc.write('\n\nChemical compound:\t' + chemical_compound)
    conc.write('\n' + distribution.to_string(header = True, index = False, justify= 'center', float_format = '{:.6e}'.format))
conc.close()
