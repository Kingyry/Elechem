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
    for i in numpy.linspace(0.5, 0, 50):
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
                'Potential [V]': numpy.around(Total_scan, 2),
                'alp_{f}': numpy.around(alpha_f, 2),
                'alp_{b}': numpy.around(alpha_b, 2),
                'E_{act. f} [eV]': E_act_f,
                'E_{act. b} [eV]': E_act_b,
                'Potential sweep': numpy.around(Potential_sweep, 2),
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
f = open(output_path + "/" + "CV.out", "w")

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
boxes = numpy.max(numpy.array([int(6*numpy.sqrt(Diffusion_constant(Independent_reactions, Database)[i][0]*iterations)+1) for i in Diffusion_constant(Independent_reactions, Database)]))
f.write('Number of boxes: {0}\n'.format('{:.0f}'.format(boxes)))
dX = numpy.max(numpy.array([numpy.sqrt(Diffusion_constant(Independent_reactions, Database)[i][1]*dT/Diffusion_constant(Independent_reactions, Database)[i][0]) for i in Diffusion_constant(Independent_reactions, Database)]))
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
f.write(EReact.to_string(header = True, index = True))

#Kinetic model and simulation calculator
kinetic_model = input('Enter electron transfer model name (BV/SMH/AMH):\t')

Reagents_fluxes = {}
Reagents_fluxes_chem = {}
Reagents_distr = {}
Reagents_distr_chem =  {}
Products_fluxes = {}
Products_fluxes_chem = {}
Products_distr = {}
Products_distr_chem = {}
Current_densities = {}

if kinetic_model in ['SMH', 'smh']:
    #add kinetics to output file
    for reactions in SMH(Independent_reactions, iterations):
        f.write('\n\nReaction {0}:\t'.format(reactions) + '\n')
        DATA = pandas.DataFrame(SMH(Independent_reactions, iterations)[reactions])
        f.write(DATA.to_string(header = True, index = False))
    for reactions in SMH(Dependent_reactions, iterations):
        f.write('\n\nReaction {0}:\t'.format(reactions) + '\n')
        DATA = pandas.DataFrame(SMH(Dependent_reactions, iterations)[reactions])
        f.write(DATA.to_string(header = True, index = False))

    #calculate reaction kinetics and system responces
    kinetics_independent = SMH(Independent_reactions, iterations)
    kinetics_dependent = SMH(Dependent_reactions, iterations)
    #independent reactions
    #Electroreduction
    for reactions in enumerate(Independent_reactions.index):
        Current_densities[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]] = numpy.zeros(iterations)
        if '.' in Independent_reactions['Reagent SMILES'].iloc[reactions[0]]:
            Reagents = Independent_reactions['Reagent SMILES'].iloc[reactions[0]].split('.')
        else:
            Reagents = [Independent_reactions['Reagent SMILES'].iloc[reactions[0]]]
        for it in Reagents:
            Reagents_distr[it] = numpy.ones((iterations, boxes))*initial_concentrations[it]
            Reagents_distr_chem[it] = numpy.ones((iterations, boxes))*initial_concentrations[it]
            Reagents_fluxes[it] = numpy.zeros(iterations)
            Reagents_fluxes_chem[it] = numpy.zeros(iterations)
        if '.' in Independent_reactions['Product SMILES'].iloc[reactions[0]]:
            Products = Independent_reactions['Product SMILES'].iloc[reactions[0]].split('.')
        else:
            Products = [Independent_reactions['Product SMILES'].iloc[reactions[0]]]
        for it in Products:
            Products_distr[it] = numpy.zeros((iterations, boxes))
            Products_distr_chem[it] = numpy.zeros((iterations, boxes))
            Products_fluxes[it] = numpy.zeros(iterations)
            Products_fluxes_chem[it] = numpy.zeros(iterations)
        #calculate distributions
        for i1 in range(iterations):
            for i2 in range(1, boxes-1):
                if i1 != 0:
                    for it_r in Reagents:
                        #distributions without and with chemical reactions
                        Reagents_distr[it_r][i1, i2] = Reagents_distr[it_r][i1-1, i2] + Diffusion_constant(Independent_reactions, Database)[it_r][0]*(Reagents_distr[it_r][i1-1, i2+1] - 2*Reagents_distr[it_r][i1-1, i2] + Reagents_distr[it_r][i1-1, i2-1])
                        Reagents_distr_chem[it_r][i1, i2] = Reagents_distr_chem[it_r][i1-1, i2] + Diffusion_constant(Independent_reactions, Database)[it_r][0]*(Reagents_distr_chem[it_r][i1-1, i2+1] - 2*Reagents_distr_chem[it_r][i1-1, i2] + Reagents_distr_chem[it_r][i1-1, i2-1])
                    for it_p in Products:
                        Products_distr[it_p][i1, i2] = Products_distr[it_p][i1-1, i2] + Diffusion_constant(Independent_reactions, Database)[it_p][0]*(Products_distr[it_p][i1-1, i2+1] - 2*Products_distr[it_p][i1-1, i2] + Products_distr[it_p][i1-1, i2-1])
                        Products_distr_chem[it_p][i1, i2] = Products_distr_chem[it_p][i1-1, i2] + Diffusion_constant(Independent_reactions, Database)[it_p][0]*(Products_distr_chem[it_p][i1-1, i2+1] - 2*Products_distr_chem[it_p][i1-1, i2] + Products_distr_chem[it_p][i1-1, i2-1])
            for it_rf in Reagents:
                for it_pf in Products:
                    #fluxes without chemical reactions
                    Reagents_fluxes[it_rf][i1] = -(kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*Reagents_distr[it_rf][i1, 1] - kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*Products_distr[it_pf][i1, 1])/(1 + kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_rf][1] + kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_pf][1])
                    Products_fluxes[it_pf][i1] = -Reagents_fluxes[it_rf][i1]
                    #fluxes with chemical reactions
                    Reagents_fluxes_chem[it_rf][i1] = -(kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*Reagents_distr_chem[it_rf][i1, 1] - kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*Products_distr_chem[it_pf][i1, 1])/(1 + kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_rf][1] + kinetics_independent[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_pf][1])
                    Products_fluxes_chem[it_pf][i1] = -Reagents_fluxes_chem[it_rf][i1]
                    #surface concentration change
                    Reagents_distr[it_rf][i1, 0] = Reagents_distr[it_rf][i1, 1] + Reagents_fluxes[it_rf][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_rf][1]
                    Reagents_distr_chem[it_rf][i1, 0] = Reagents_distr_chem[it_rf][i1, 1] + Reagents_fluxes_chem[it_rf][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_rf][1]
                    Products_distr[it_pf][i1, 0] = Products_distr[it_pf][i1, 1] + Products_fluxes[it_pf][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_pf][1]
                    Products_distr_chem[it_pf][i1, 0] = Products_distr_chem[it_pf][i1, 1] + Products_fluxes_chem[it_pf][i1]*dX/Diffusion_constant(Independent_reactions, Database)[it_pf][1]

                #Current densities
                Current_densities[Independent_reactions['Reaction SMILES'].iloc[reactions[0]]][i1] = faraday*1000*Reagents_fluxes_chem[it_rf][i1]
    #dependent reactions
    #Electroreduction
    for reactions in enumerate(Dependent_reactions.index):
        Current_densities[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]] = numpy.zeros(iterations)
        if '.' in Dependent_reactions['Reagent SMILES'].iloc[reactions[0]]:
            Reagents = Dependent_reactions['Reagent SMILES'].iloc[reactions[0]].split('.')
        else:
            Reagents = [Dependent_reactions['Reagent SMILES'].iloc[reactions[0]]]
        for it in Reagents:
            Reagents_distr[it] = numpy.zeros((iterations, boxes))
            Reagents_distr_chem[it] = numpy.zeros((iterations, boxes))
            Reagents_fluxes[it] = numpy.zeros(iterations)
            Reagents_fluxes_chem[it] = numpy.zeros(iterations)
        if '.' in Dependent_reactions['Product SMILES'].iloc[reactions[0]]:
            Products = Dependent_reactions['Product SMILES'].iloc[reactions[0]].split('.')
        else:
            Products = [Dependent_reactions['Product SMILES'].iloc[reactions[0]]]
        for it in Products:
            Products_distr[it] = numpy.zeros((iterations, boxes))
            Products_distr_chem[it] = numpy.zeros((iterations, boxes))
            Products_fluxes[it] = numpy.zeros(iterations)
            Products_fluxes_chem[it] = numpy.zeros(iterations)
        #calculate distributions
        for i1 in range(iterations):
            for i2 in range(1, boxes-1):
                if i1 != 0:
                    for it_r in Reagents:
                        #distributions without chemical reactions
                        Reagents_distr[it_r][i1, i2] = Reagents_distr[it_r][i1-1, i2] + Diffusion_constant(Dependent_reactions, Database)[it_r][0]*(Reagents_distr[it_r][i1-1, i2+1] - 2*Reagents_distr[it_r][i1-1, i2] + Reagents_distr[it_r][i1-1, i2-1])
                    for it_p in Products:
                        Products_distr[it_p][i1, i2] = Products_distr[it_p][i1-1, i2] + Diffusion_constant(Dependent_reactions, Database)[it_p][0]*(Products_distr[it_p][i1-1, i2+1] - 2*Products_distr[it_p][i1-1, i2] + Products_distr[it_p][i1-1, i2-1])

                    for it_r in Reagents:
                        #distributions with chemical reactions
                        Reagents_distr_chem[it_r][i1, i2] = Reagents_distr_chem[it_r][i1-1, i2] + Diffusion_constant(Dependent_reactions, Database)[it_r][0]*(Reagents_distr_chem[it_r][i1-1, i2+1] - 2*Reagents_distr_chem[it_r][i1-1, i2] + Reagents_distr_chem[it_r][i1-1, i2-1])
                    for it_p in Products:
                        Products_distr_chem[it_p][i1, i2] = Products_distr_chem[it_p][i1-1, i2] + Diffusion_constant(Dependent_reactions, Database)[it_p][0]*(Products_distr_chem[it_p][i1-1, i2+1] - 2*Products_distr_chem[it_p][i1-1, i2] + Products_distr_chem[it_p][i1-1, i2-1])


            for it_rf in Reagents:
                for it_pf in Products:
                    #fluxes without chemical reactions
                    Reagents_fluxes[it_rf][i1] = -(kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*Reagents_distr[it_rf][i1, 1] - kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*Products_distr[it_pf][i1, 1])/(1 + kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*dX/Diffusion_constant(Dependent_reactions, Database)[it_rf][1] + kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*dX/Diffusion_constant(Dependent_reactions, Database)[it_pf][1])
                    Products_fluxes[it_pf][i1] = -Reagents_fluxes[it_rf][i1]
                    #fluxes with chemical reactions
                    Reagents_fluxes_chem[it_rf][i1] = -(kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*Reagents_distr_chem[it_rf][i1, 1] - kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*Products_distr_chem[it_pf][i1, 1])/(1 + kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{f}'][i1]*dX/Diffusion_constant(Dependent_reactions, Database)[it_rf][1] + kinetics_dependent[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]]['k_{b}'][i1]*dX/Diffusion_constant(Dependent_reactions, Database)[it_pf][1])
                    Products_fluxes_chem[it_pf][i1] = -Reagents_fluxes_chem[it_rf][i1]
                    #surface concentration change
                    Reagents_distr[it_rf][i1, 0] = Reagents_distr[it_rf][i1, 1] + (Reagents_fluxes[it_rf][i1] + Products_fluxes[it_rf][i1])*dX/Diffusion_constant(Dependent_reactions, Database)[it_rf][1]
                    Reagents_distr_chem[it_rf][i1, 0] = Reagents_distr_chem[it_rf][i1, 1] + (Reagents_fluxes_chem[it_rf][i1] + Products_fluxes_chem[it_rf][i1])*dX/Diffusion_constant(Dependent_reactions, Database)[it_rf][1]
                    Products_distr[it_pf][i1, 0] = Products_distr[it_pf][i1, 1] + Products_fluxes[it_pf][i1]*dX/Diffusion_constant(Dependent_reactions, Database)[it_pf][1]
                    Products_distr_chem[it_pf][i1, 0] = Products_distr_chem[it_pf][i1, 1] + Products_fluxes_chem[it_pf][i1]*dX/Diffusion_constant(Dependent_reactions, Database)[it_pf][1]

                #Current densities
                Current_densities[Dependent_reactions['Reaction SMILES'].iloc[reactions[0]]][i1] = faraday*1000*Reagents_fluxes_chem[it_rf][i1]


f.close()


import matplotlib.pyplot as plt

cm = 1/2.54
figure_conc = plt.figure(figsize = (8.3*cm, 8.3*cm), dpi = 600)
ax1 = figure_conc.add_axes([0.15, 0.2, 0.56, 0.45])
for i in Current_densities:
    ax1.plot(Total_scan, Current_densities[i])
plt.show()
