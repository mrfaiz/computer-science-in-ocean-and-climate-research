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
alph_array = np.arange(0,2,stepSize)

print(alph_array)

# ensembleResults = []
# ensembleParameterValues = []

# # Removing old files from the output folder
# os.system('rm output/*')

# alphaGrid = [x*stepSize for x in range(0,n)]
# # print("Grid for alpha:")
# # print([round(x,2) for x in alphaGrid])

# # Compiling the source
# print("Compiling mod_precision")
# os.system('gfortran -c mod_precision.f90')
# print("Compiling twoThreadsSections")
# os.system('gfortran -fopenmp twoThreadsSections.f95')

# i = 0
# for alpha in alphaGrid:
#     # Adjust the NML file for the current value of alpha
#     writeParameterToNML('predatorprey.nml','alpha',alpha)

#     # Perform the simulation
#     print("Starting ensemble run for alpha: ", alpha)
#     os.system('./a.out')

#     # Save selected results
#     df = pd.read_csv('results.csv')
#     headers = list(df)
#     prey = df[headers[2]]
#     ensembleResults += [prey[len(prey)-1]]
#     ensembleParameterValues += [alpha]

#     # Move the output file to subdirectory
#     os.system('mv results.csv output/results_'+str(i)+'.csv')
#     i+=1


# # Save the results to a csv file
# df = pd.DataFrame({'alpha': ensembleParameterValues,
#                    'prey_population': ensembleResults})
# df.to_csv('ensembleResults.csv', index=False)
