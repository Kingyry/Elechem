#imports
import sys
sys.path.append('src')

import numpy
import pandas
from data import faraday, gas_constant

#Set environment parameters
#import of initial parameters

def Diffusion_constant(Reaction_database, Database):
    elements = []
    Reaction_database = Reaction_database.sort_values(by = ['Reduction potential, V'], ascending=False).reset_index(drop = True)
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
    for i in numpy.linspace(0.5, 0, 10):
        D = i*numpy.array(Diff_coeff)/Diff_coeff[0]
        if all(D < 0.5):
            break
    return elements, Diff_coeff, D

print("***Set reactions and environment conditions***\n\n")
path = input('Input path to your reaction database:\t')
Reaction_database = pandas.read_csv(path)
Reaction_database = Reaction_database.sort_values(by = ['Reduction potential, V'], ascending=False).reset_index(drop = True)
path = input('Input path to your database:\t')
Database = pandas.read_csv(path)

output_path = input('Enter output path:\t')

f = open(output_path + "/" + "CV.out", "w")
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
boxes = int(6*numpy.sqrt(Diffusion_constant(Reaction_database, Database)[2][0]*iterations)+1)
f.write('Number of boxes: {0}\n'.format('{:.0f}'.format(boxes)))
dX = numpy.sqrt(Diffusion_constant(Reaction_database, Database)[1][0]*dT/Diffusion_constant(Reaction_database, Database)[2][0])
X = dX*boxes

distance = [0]

for i in range(boxes):
    if i > 0:
        distance.append(distance[i-1] + dX)

f.write('Simulation region [cm]: {0}\n'.format('{:.2f}'.format(X)))
if Reaction_database['Reduction potential, V'].iloc[0] < 0:
    initial_conc = float(input('Enter initial concentration of {0} [mol/cm3]:\t'.format(Reaction_database['Reagent SMILES'].iloc[0])))
else:
    initial_conc = float(input('Enter initial concentration of {0} [mol/cm3]:\t'.format(Reaction_database['Product SMILES'].iloc[0])))

f.write('{0} initial concentration [mol/cm3]: {1}\n\n'.format(Reaction_database['Reagent SMILES'].iloc[0], float('{:.4f}'.format(initial_conc))))

Forward_scan = numpy.linspace(Potential_max, Potential_min,
                              int(iterations/2), endpoint=False)
Backward_scan = numpy.linspace(Potential_min, Potential_max, int(iterations/2))
Total_scan = numpy.concatenate((Forward_scan, Backward_scan))

f.write("\nElectrochemical reactions:\n")
EReact = Reaction_database.loc[Reaction_database['Reaction'] == 'electrochem.'].drop(columns = ['Reagent SMILES', 'Reagent multiplicity', 'Reagent charge', 'Reagent energy, eV', 'Product SMILES', 'Product multiplicity', 'Product charge', 'Product energy, eV', 'Functional', 'Basis', 'Solvent', 'Temperature, K', 'Reaction', 'Type', 'Mechanism', 'lambda (solv.), eV', 'lambda (inner.), eV', 'lambda (tot.), eV', 'k0, cm/s'])
EReact = EReact.drop(EReact.columns[[0]], axis = 1)
EReact = EReact.sort_values(by = ['Reduction potential, V'], ascending=False).reset_index(drop = True)
f.write(EReact.to_string(header = True, index = True))


Reagent_conc = []
Product_conc = []
Reagent_fluxes = []
Product_fluxes = []
Current_dens = []

for i in enumerate(Reaction_database.loc[Reaction_database['Reaction'] == 'electrochem.'].index):
    if Reaction_database['Reaction'].iloc[i[0]] == 'electrochem.' and Reaction_database['Reduction potential, V'].iloc[0] < 0:
        I = numpy.arange(0, iterations, 1)

        Potential_sweep = (Total_scan-Reaction_database['Reduction potential, V'].iloc[i[0]])*faraday/gas_constant/Reaction_database['Temperature, K'].iloc[i[0]]

        alpha_f = 0.5 + (Total_scan - Reaction_database['Reduction potential, V'].iloc[i[0]])/4/Reaction_database['lambda (tot.), eV'].iloc[i[0]]
        E_act_f = Reaction_database['lambda (tot.), eV'].iloc[i[0]]/4 + alpha_f*(Total_scan - Reaction_database['Reduction potential, V'].iloc[i[0]])
        kf = Reaction_database['k0, cm/s'].iloc[i[0]]*numpy.exp(-alpha_f*Potential_sweep)
        if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
            alpha_b = numpy.ones(iterations)*numpy.nan
            E_act_b = numpy.ones(iterations)*numpy.nan
            kb = numpy.ones(iterations)*numpy.nan
        else:
            alpha_b = (1 - alpha_f)
            E_act_b = Reaction_database['lambda (tot.), eV'].iloc[i[0]]/4 - alpha_b*(Total_scan-Reaction_database['Reduction potential, V'].iloc[i[0]])
            kb = Reaction_database['k0, cm/s'].iloc[i[0]]*numpy.exp(alpha_b*Potential_sweep)
        f.write('\n\nReaction {0}:\t'.format(i[0]) + Reaction_database['Reaction SMILES'].iloc[i[0]] + '\n')
    elif Reaction_database['Reaction'].iloc[i[0]] == 'electrochem.' and Reaction_database['Reduction potential, V'].iloc[0] > 0:
        I = numpy.arange(0, iterations, 1)

        Potential_sweep = (Total_scan-Reaction_database['Reduction potential, V'].iloc[i[0]])*faraday/gas_constant/Reaction_database['Temperature, K'].iloc[i[0]]

        alpha_b = 0.5 + (Total_scan - Reaction_database['Reduction potential, V'].iloc[i[0]])/4/Reaction_database['lambda (tot.), eV'].iloc[i[0]]
        E_act_b = Reaction_database['lambda (tot.), eV'].iloc[i[0]]/4 + alpha_f*(Total_scan-Reaction_database['Reduction potential, V'].iloc[i[0]])
        kb = Reaction_database['k0, cm/s'].iloc[i[0]]*numpy.exp(-alpha_f*Potential_sweep)
        if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
            alpha_f = numpy.ones(iterations)*numpy.nan
            E_act_f = numpy.ones(iterations)*numpy.nan
            kb = numpy.ones(iterations)*numpy.nan
        else:
            alpha_f = (1 - alpha_f)
            E_act_f = Reaction_database['lambda (tot.), eV'].iloc[i[0]]/4 - alpha_b*(Total_scan-Reaction_database['Reduction potential, V'].iloc[i[0]])
            kf = Reaction_database['k0, cm/s'].iloc[i[0]]*numpy.exp(alpha_f*Potential_sweep)
        f.write('\n\nReaction {0}:\t'.format(i[0]) + Reaction_database['Reaction SMILES'].iloc[i[0]] + '\n')
    data = {'Iterations': I,
            'Potential [V]': numpy.around(Total_scan, 2),
            'alp_{f}': numpy.around(alpha_f, 2),
            'alp_{b}': numpy.around(alpha_b, 2),
            'E_{act. f} [eV]': numpy.around(E_act_f, 2),
            'E_{act. b} [eV]': numpy.around(E_act_b, 2),
            'Potential sweep': numpy.around(Potential_sweep, 2),
            'k_{f}': kf,
            'k_{b}': kb}
    database = pandas.DataFrame(data)
    f.write(database.to_string(header = True, index = False))

    if i[0] < 1:
        Current_density = numpy.zeros(iterations)
        Reagent_flux = numpy.zeros(iterations)
        Product_flux = numpy.zeros(iterations)
        Reagent_distr = numpy.ones((iterations, boxes))*initial_conc

        if Reaction_database['Reduction potential, V'].iloc[0] < 0:
            Reagent_SMILES = str(Reaction_database['Reagent SMILES'].iloc[i[0]])
            Diff_reagent = Database['Diffusion, cm2/s'].loc[Database['SMILES'] == Reagent_SMILES].values[0]

            if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
                Products_SMILES = [i for i in Reaction_database['Product SMILES'].iloc[i[0]].split('.')]
            else:
                Products_SMILES = [Reaction_database['Product SMILES'].iloc[i[0]]]
            D_products = []
            Diff_products = []
            for it in enumerate(Diffusion_constant(Reaction_database, Database)[0]):
                if it[1] == Reagent_SMILES:
                    D_reagent = Diffusion_constant(Reaction_database, Database)[2][it[0]]
                    Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])

                for l in Products_SMILES:
                    if it[1] == l:
                        D_products.append(Diffusion_constant(Reaction_database, Database)[2][it[0]])
                        Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])


            if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
                Products_distr = [numpy.zeros((iterations, boxes)) for i in range(len(Reaction_database['Product SMILES'].iloc[i[0]].split('.')))]
            else:
                Products_distr = [numpy.zeros((iterations, boxes))]

        #ELECTROOXYDATION
        elif Reaction_database['Reduction potential, V'].iloc[0] > 0:
            Reagent_SMILES = str(Reaction_database['Product SMILES'].iloc[i[0]])
            Diff_reagent = Database['Diffusion, cm2/s'].loc[Database['SMILES'] == Reagent_SMILES].values[0]

            if '.' in Reaction_database['Reagent SMILES'].iloc[i[0]]:
                Products_SMILES = [i for i in Reaction_database['Reagent SMILES'].iloc[i[0]].split('.')]
            else:
                Products_SMILES = [Reaction_database['Reagent SMILES'].iloc[i[0]]]
            D_products = []
            Diff_products = []
            for it in enumerate(Diffusion_constant(Reaction_database, Database)[0]):
                if it[1] == Reagent_SMILES:
                    D_reagent = Diffusion_constant(Reaction_database, Database)[2][it[0]]
                    Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])

                for l in Products_SMILES:
                    if it[1] == l:
                        D_products.append(Diffusion_constant(Reaction_database, Database)[2][it[0]])
                        Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])


            if '.' in Reaction_database['Reagent SMILES'].iloc[i[0]]:
                Products_distr = [numpy.zeros((iterations, boxes)) for i in range(len(Reaction_database['Reagent SMILES'].iloc[i[0]].split('.')))]
            else:
                Products_distr = [numpy.zeros((iterations, boxes))]














            for i1 in range(iterations):
                for i2 in range(1, boxes-1):
                    if i1 != 0:
                        Reagent_distr[i1, i2] = Reagent_distr[i1-1, i2] + D_reagent*(Reagent_distr[i1-1, i2+1] - 2*Reagent_distr[i1-1, i2] + Reagent_distr[i1-1, i2-1])

                for i2 in range(1, boxes-1):
                    if i1 != 0:
                        for k1 in enumerate(Products_distr):
                            k1[1][i1, i2] = k1[1][i1-1, i2] + D_products[k1[0]]*(k1[1][i1-1, i2+1] - 2*k1[1][i1-1, i2] + k1[1][i1-1, i2-1])

                if len(Products_SMILES) > 1:
                    Reagent_flux[i1] = -(kf[i1]*Reagent_distr[i1, 1])/(1 + kf[i1]*dX/Diff_reagent)
                else:
                    Reagent_flux[i1] = -(kf[i1]*Reagent_distr[i1, 1]-kb[i1]*Products_distr[0][i1, 1])/(1 + kf[i1]*dX/Diff_reagent + kb[i1]*dX/Diff_products[0])

                Product_flux[i1] = -Reagent_flux[i1]
                Reagent_distr[i1, 0] = Reagent_distr[i1, 1] + Reagent_flux[i1]*dX/Diff_reagent
                for k1 in enumerate(Products_distr):
                    k1[1][i1, 0] = k1[1][i1, 1] + Product_flux[i1]*dX/Diff_products[k1[0]]
                Current_density[i1] = faraday*1000*Reagent_flux[i1]


        Reagent_database = pandas.DataFrame(Reagent_distr, columns=distance)
        Reagent_database.to_csv(output_path + "/" + Reagent_SMILES + ".csv")
        f.write('\nReaction {0} output:\t'.format(int(i[0])) + output_path + "/" + Reagent_SMILES + ".csv\n")
        for j in enumerate(Products_distr):
            Products_database = pandas.DataFrame(j[1], columns=distance)
            Products_database.to_csv(output_path + "/" + Products_SMILES[j[0]] + ".csv")
            f.write('Reaction {0} output:\t'.format(int(i[0])) + output_path + "/" + Products_SMILES[j[0]] + ".csv\n")

        Reagent_fluxes.append(Reagent_flux)
        Product_fluxes.append(Product_flux)
        Current_dens.append(Current_density)


    else:
        if '.' in Reaction_database['Product SMILES'].iloc[i[0]-1]:
            SMILES = Reaction_database['Product SMILES'].iloc[i[0]-1].split('.')
        else:
            SMILES = Reaction_database['Product SMILES'].iloc[i[0]-1]

        if Reaction_database['Reagent SMILES'].iloc[i[0]] in SMILES:
            Current_density = numpy.zeros(iterations)
            Reagent_flux = numpy.zeros(iterations)
            Product_flux = numpy.zeros(iterations)
            Reagent_distr = numpy.zeros((iterations, boxes))


            if Reaction_database['Reduction potential, V'].iloc[i[0]] < 0:
                Reagent_SMILES = str(Reaction_database['Reagent SMILES'].iloc[i[0]])
                Diff_reagent = Database['Diffusion, cm2/s'].loc[Database['SMILES'] == Reagent_SMILES].values[0]

                if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
                    Products_SMILES = [i for i in Reaction_database['Product SMILES'].iloc[i[0]].split('.')]
                else:
                    Products_SMILES = [Reaction_database['Product SMILES'].iloc[i[0]]]
                D_products = []
                Diff_products = []
                for it in enumerate(Diffusion_constant(Reaction_database, Database)[0]):
                    if it[1] == Reagent_SMILES:
                        D_reagent = Diffusion_constant(Reaction_database, Database)[2][it[0]]
                        Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])

                    for l in Products_SMILES:
                        if it[1] == l:
                            D_products.append(Diffusion_constant(Reaction_database, Database)[2][it[0]])
                            Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])


                if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
                    Products_distr = [numpy.zeros((iterations, boxes)) for i in range(len(Reaction_database['Product SMILES'].iloc[i[0]].split('.')))]
                else:
                    Products_distr = [numpy.zeros((iterations, boxes))]

                for i1 in range(iterations):
                    Reagent_distr[i1, 0] = Reagent_distr[i1, 1] + (Reagent_flux[i1] + Product_fluxes[i[0]-1][i1])*dX/Diff_reagent

                    for i2 in range(1, boxes-1):
                        if i1 != 0:
                            Reagent_distr[i1, i2] = Reagent_distr[i1-1, i2] + D_reagent*(Reagent_distr[i1-1, i2+1] - 2*Reagent_distr[i1-1, i2] + Reagent_distr[i1-1, i2-1])

                    if len(Products_SMILES) > 1:
                        Reagent_flux[i1] = -(kf[i1]*Reagent_distr[i1, 1])/(1 + kf[i1]*dX/Diff_reagent)
                    else:
                        Reagent_flux[i1] = -(kf[i1]*Reagent_distr[i1, 1]-kb[i1]*Products_distr[0][i1, 1])/(1 + kf[i1]*dX/Diff_reagent + kb[i1]*dX/Diff_products[0])
                    Product_flux[i1] = -Reagent_flux[i1]
                    Current_density[i1] = faraday*1000*Reagent_flux[i1]


                    for k1 in enumerate(Products_distr):
                        k1[1][i1, 0] = k1[1][i1, 1] + Product_flux[i1]*dX/Diff_products[k1[0]]
                    for i2 in range(1, boxes-1):
                        if i1 != 0:
                            for k1 in enumerate(Products_distr):
                                k1[1][i1, i2] = k1[1][i1-1, i2] + D_products[k1[0]]*(k1[1][i1-1, i2+1] - 2*k1[1][i1-1, i2] + k1[1][i1-1, i2-1])



            Reagent_database = pandas.DataFrame(Reagent_distr, columns=distance)
            Reagent_database.to_csv(output_path + "/" + Reagent_SMILES + ".csv")
            f.write('\nReaction {0} output:\t'.format(int(i[0])) + output_path + "/" + Reagent_SMILES + ".csv\n")
            for j in enumerate(Products_distr):
                Products_database = pandas.DataFrame(j[1], columns=distance)
                Products_database.to_csv(output_path + "/" + Products_SMILES[j[0]] + ".csv")
                f.write('Reaction {0} output:\t'.format(int(i[0])) + output_path + "/" + Products_SMILES[j[0]] + ".csv\n")

            Reagent_fluxes.append(Reagent_flux)
            Product_fluxes.append(Product_flux)
            Current_dens.append(Current_density)


        else:
            Current_density = numpy.zeros(iterations)
            Reagent_flux = numpy.zeros(iterations)
            Product_flux = numpy.zeros(iterations)
            Reagent_distr = numpy.ones((iterations, boxes))*initial_conc

            if Reaction_database['Reduction potential, V'].iloc[i[0]] < 0:
                Reagent_SMILES = str(Reaction_database['Reagent SMILES'].iloc[i[0]])
                Diff_reagent = Database['Diffusion, cm2/s'].loc[Database['SMILES'] == Reagent_SMILES].values[0]

                if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
                    Products_SMILES = [i for i in Reaction_database['Product SMILES'].iloc[i[0]].split('.')]
                else:
                    Products_SMILES = [Reaction_database['Product SMILES'].iloc[i[0]]]
                D_products = []
                Diff_products = []
                for it in enumerate(Diffusion_constant(Reaction_database, Database)[0]):
                    if it[1] == Reagent_SMILES:
                        D_reagent = Diffusion_constant(Reaction_database, Database)[2][it[0]]
                        Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])

                    for l in Products_SMILES:
                        if it[1] == l:
                            D_products.append(Diffusion_constant(Reaction_database, Database)[2][it[0]])
                            Diff_products.append(Diffusion_constant(Reaction_database, Database)[1][it[0]])


                if '.' in Reaction_database['Product SMILES'].iloc[i[0]]:
                    Products_distr = [numpy.zeros((iterations, boxes)) for i in range(len(Reaction_database['Product SMILES'].iloc[i[0]].split('.')))]
                else:
                    Products_distr = [numpy.zeros((iterations, boxes))]

                for i1 in range(iterations):
                    for i2 in range(1, boxes-1):
                        if i1 != 0:
                            Reagent_distr[i1, i2] = Reagent_distr[i1-1, i2] + D_reagent*(Reagent_distr[i1-1, i2+1] - 2*Reagent_distr[i1-1, i2] + Reagent_distr[i1-1, i2-1])

                    for i2 in range(1, boxes-1):
                        if i1 != 0:
                            for k1 in enumerate(Products_distr):
                                k1[1][i1, i2] = k1[1][i1-1, i2] + D_products[k1[0]]*(k1[1][i1-1, i2+1] - 2*k1[1][i1-1, i2] + k1[1][i1-1, i2-1])

                    if len(Products_SMILES) > 1:
                        Reagent_flux[i1] = -(kf[i1]*Reagent_distr[i1, 1])/(1 + kf[i1]*dX/Diff_reagent)
                    else:
                        Reagent_flux[i1] = -(kf[i1]*Reagent_distr[i1, 1]-kb[i1]*Products_distr[0][i1, 1])/(1 + kf[i1]*dX/Diff_reagent + kb[i1]*dX/Diff_products[0])

                    Product_flux[i1] = -Reagent_flux[i1]
                    Reagent_distr[i1, 0] = Reagent_distr[i1, 1] + Reagent_flux[i1]*dX/Diff_reagent
                    for k1 in enumerate(Products_distr):
                        k1[1][i1, 0] = k1[1][i1, 1] + Product_flux[i1]*dX/Diff_products[k1[0]]
                    Current_density[i1] = faraday*1000*Reagent_flux[i1]

            Reagent_database = pandas.DataFrame(Reagent_distr, columns=distance)
            Reagent_database.to_csv(output_path + "/" + Reagent_SMILES + ".csv")
            f.write('\nReaction {0} output:\t'.format(int(i[0])) + output_path + "/" + Reagent_SMILES + ".csv\n")
            for j in enumerate(Products_distr):
                Products_database = pandas.DataFrame(j[1], columns=distance)
                Products_database.to_csv(output_path + "/" + Products_SMILES[j[0]] + ".csv")
                f.write('Reaction {0} output:\t'.format(int(i[0])) + output_path + "/" + Products_SMILES[j[0]] + ".csv\n")

            Reagent_fluxes.append(Reagent_flux)
            Product_fluxes.append(Product_flux)
            Current_dens.append(Current_density)



Current_dens = numpy.array([elem for elem in Current_dens if numpy.sum(elem) != 0]).T

column = ['Reaction {0}'.format(int(i[0])) for i in enumerate(Reaction_database['Reaction'].loc[Reaction_database['Reaction'] == 'electrochem.'].index)]

column = ['Iterations', 'Potential [V]'] + column

Currents = numpy.concatenate((numpy.array([numpy.arange(0, iterations, 1), [float('{:.2f}'.format(i)) for i in Total_scan]]).T, Current_dens), axis = 1)
Currents = pandas.DataFrame(Currents, columns = column)
Currents['Iterations'] = numpy.arange(0, iterations, 1)
Currents['Total'] = Currents.iloc[:,2:].sum(axis = 1)
f.write('\nCurrent densities from electrochemical reactions J [mA/cm2]')
f.write('\n\n' + Currents.to_string(header = True, index = False))

f.close()
