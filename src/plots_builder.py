# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 18:40:23 2023

@author: Aleksandr
"""

#imports
import sys
sys.path.append('src')

import matplotlib.pyplot as plt
import pandas

#functions
def Ionization_energy(path):
    database = pandas.read_csv(path)
    cm = 1/2.54
    figure_conc = plt.figure(figsize = (8.3*cm, 8.3*cm), dpi = 600)
    ax1 = figure_conc.add_axes([0.15, 0.2, 0.56, 0.45])
    ax1.bar(database['Functional'], database['Free energy, eV'], color = 'white', edgecolor='black', linewidth = 0.5, zorder = 2, hatch = '//'*2)
    ax1.grid(True, color = 'black', linestyle = '-', linewidth = 0.3, zorder = 0, alpha = 0.3)
    ax1.set_ylim(database['Free energy, eV'].max()*1.25, database['Free energy, eV'].min()*0.75)
    ax1.set_ylabel('Ionization energy [eV]', fontsize=8)
    ax1.set_xlabel('Functional', fontsize=8)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        label.set_fontsize(6)
    plt.xticks(rotation = 45, ha = 'right')
    plt.savefig(input('Output path:\t') + '/' + 'Ionization_energy.png')

welcome_window = """Choose plot
    1. Ionization energy plot.
"""

print(welcome_window)
while True:
    enter = int(input("Input number:\t"))
    if enter == 1:
        path = input('Input path to your reaction database:\t')
        Ionization_energy(path)
    else:
        print("Back to main menu")
        break
