import f90nml
import time
import os
import pandas as pd
import numpy as np 

def writeParameterToNML(filename, name, value):
    nml = f90nml.read(filename)
    nml['model_parameters'][name] = value
    os.remove(filename)
    nml.write(filename)

# Defining the grid
T = 20
min = 0
max = 2
stepSize = (max-min)/(T)
alph_array = np.arange(0,2+stepSize,stepSize)

ensembleResults = []
ensembleParameterValues = []

print("Compiling preditor_prey")
os.system('gfortran -fopenmp preditor_prey.f95')

i = 0
for alpha in alph_array:
    # Adjust the NML file for the current value of alpha
    writeParameterToNML('/home/faiz/SS_2020/Ocean/cs_in_ocean_and_climate_research_python_fortran/exercies/10/predatorprey.nml','alpha',alpha)

    # Perform the simulation
    print("Starting ensemble run for alpha: ", alpha)
    os.system('./a.out')
