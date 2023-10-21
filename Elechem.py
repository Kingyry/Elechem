# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 11:18:13 2023

@author: Aleksandr
"""
#imports

from output import welcome_window
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
    else:
        print("Exit")
        break
