# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 11:18:13 2023

@author: Aleksandr
"""
#imports

welcome_window="""Choose an option
    1. Create/Refresh database.
    2. Database manager.
    3. Cyclic voltammetry.
    0. Exit.
           """
#Database building
while True:
    print(welcome_window)
    enter = int(input("Input number:\t"))
    if enter == 1:
        with open("src/database_builder.py") as f:
            exec(f.read())
    elif enter == 2:
        with open("src/data_manager.py") as f:
            exec(f.read())
    elif enter == 3:
        with open("src/CV.py") as f:
            exec(f.read())
    else:
        print("Exit")
        break
