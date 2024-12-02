import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import re
import pandas as pd

#load the yaml file
filename = 'fabm.yaml'


#flag to add nutrients to adjacency matrix
flagN = True

with open(filename, 'r') as file:
    texto = file.readlines()

total_FTs = [i for i, line in enumerate(texto) if re.search("model: ogs/(Phyto|PelBac|MicroZoo|MesoZoo)", line)]
bloque = [texto[i - 3] for i in total_FTs] + [texto[i - 2] for i in total_FTs] + [texto[i - 1] for i in total_FTs]

nombres_FTs = [line.strip() for line in bloque if re.search(r"[A-Z][0-9]:|[A-Z][0-9]_", line)]
nombres_FTs = [line.replace(":", "") for line in nombres_FTs]

# Get 'clasesG' (exponents for size classes)
inicio = [line.find('_') + 1 for line in nombres_FTs]
clasesG = [int(line[start:start + 2].replace('m', '-')) for line, start in zip(nombres_FTs, inicio)]

# Phyto large to small, B, Zoo small to large
nombres_ordenados = (
    sorted([name for name in nombres_FTs if "P4" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P1" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P5" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P2" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P8" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P7" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P3" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P9" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "P6" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "B1" in name], key=lambda x: clasesG[nombres_FTs.index(x)], reverse=True) +
    sorted([name for name in nombres_FTs if "Z6" in name], key=lambda x: clasesG[nombres_FTs.index(x)]) +
    sorted([name for name in nombres_FTs if "Z5" in name], key=lambda x: clasesG[nombres_FTs.index(x)]) +
    sorted([name for name in nombres_FTs if "Z4" in name], key=lambda x: clasesG[nombres_FTs.index(x)]) +
    sorted([name for name in nombres_FTs if "Z3" in name], key=lambda x: clasesG[nombres_FTs.index(x)])
)

# Create adjacency matrix
n = len(total_FTs)
adj_matrix = np.zeros((n, n))
row_col_names = {i: name for i, name in enumerate(nombres_ordenados)}


# Define color scheme
#related to the size classes
# P4*10, P1*10, P5*10, P2*10, P8*6, P7*6, P3*3, P9*3, P6*2, B1*1, Z6*7, Z5*10, Z4*10, Z3*10
colorinchis = ["#FF69B4"] * 10 + ["#FFB6C1"] * 10 + ["#B8860B"] * 10 + ["#FFD700"] * 10 + \
              ["#008000"] * 6 + ["#6B8E23"] * 6 + ["#8FBC8F"] * 3 + ["#00BFFF"] * 3 + \
              ["#9932CC"] * 2 + ["#000000"] * 1 + ["#333333"] * 7 + ["#666666"] * 10 + \
              ["#999999"] * 10 + ["#CCCCCC"] * 10

# Fill the matrix for MicroZoo (and Nanoflagellates)
instances = [i for i, line in enumerate(texto) if re.search("model: ogs/MicroZoo", line)]
bloque = [texto[i - 3] for i in instances] + [texto[i - 2] for i in instances] + [texto[i - 1] for i in instances]
stop = [i for i, line in enumerate(texto) if re.search("model: ogs/MesoZoo", line)][0]

cazadores = [line.strip() for line in bloque if re.search(r"[A-Z][0-9]:|[A-Z][0-9]_", line)]
cazadores = [line.replace(":", "") for line in cazadores]

#adjacency matrix fill
for i in range(len(instances)):
    cazador = cazadores[i]
    #largo is need to search in between a specific FT parameters block
    if i != len(instances) - 1:
        largo = instances[i + 1] - instances[i]
    else:
        largo = stop - instances[i]

    presas = [line for line in texto[instances[i]:instances[i] + largo] if re.search(r" prey\d+:", line)]
    preferencias = [line for line in texto[instances[i]:instances[i] + largo] if "suprey" in line]
    
    nombritos = [line.split(":")[1].strip() for line in presas]
    prefis    = [re.search(r':\s*([+-]?\d*\.?\d+)', line).group(1) for line in preferencias] #get the preference values
    
    cuales = [list(row_col_names.values()).index(nom) for nom in nombritos]
    cazador_idx = list(row_col_names.values()).index(cazador)
    adj_matrix[cuales, cazador_idx] = prefis

# Fill the matrix for MesoZoo
instances = [i for i, line in enumerate(texto) if re.search("model: ogs/MesoZoo", line)]
bloque = [texto[i - 3] for i in instances] + [texto[i - 2] for i in instances] + [texto[i - 1] for i in instances]

cazadores = [line.strip() for line in bloque if re.search(r"[A-Z][0-9]:|[A-Z][0-9]_", line)]
cazadores = [line.replace(":", "") for line in cazadores]

#adjacency matrix fill
for i in range(len(instances)):
    cazador = cazadores[i]
    #largo is need to search in between a specific FT parameters block
    if i != len(instances) - 1:
        largo = instances[i + 1] - instances[i]
    else:
        largo = len(texto) - instances[i]
    
    presas = [line for line in texto[instances[i]:instances[i] + largo] if re.search(r" prey\d+:", line)]
    preferencias = [line for line in texto[instances[i]:instances[i] + largo] if "suprey" in line]
    
    nombritos = [line.split(":")[1].strip() for line in presas]
    prefis    = [re.search(r':\s*([+-]?\d*\.?\d+)', line).group(1) for line in preferencias] #get the preference values
    
    cuales = [list(row_col_names.values()).index(nom) for nom in nombritos]
    cazador_idx = list(row_col_names.values()).index(cazador)
    adj_matrix[cuales, cazador_idx] = prefis


#if nutrients flag is True add a row and column for nutrients as first FT
if flagN:
    adj_matrix = np.insert(adj_matrix, 0, 0, axis=0)
    adj_matrix = np.insert(adj_matrix, 0, 0, axis=1)
    nombres_ordenados = ['N'] + nombres_ordenados
    for iname,name in enumerate(nombres_ordenados):
        if 'B' or 'Z' in name:      #bacteria and zooplankton eat nutrients and excrete them
            adj_matrix[0, iname] = 1
            adj_matrix[iname, 0] = 1
        if 'P' in name:             #phytoplankton die in nutrients
            adj_matrix[0, iname] = 0
            adj_matrix[iname, 0] = 1
    adj_matrix[0, 0] = 0
    #add nutrient to size classes and colors
    clasesG = [np.mean(clasesG)] + clasesG
    colorinchis = ["#FF4500"] + colorinchis


#save the adjacency matrix in a csv file with as header and index the FT names
df = pd. DataFrame(adj_matrix, index=nombres_ordenados, columns=nombres_ordenados)
if flagN:
    df.to_csv('adj_matrix_nutrients.csv')
else:
    df.to_csv('adj_matrix.csv')

# Create the graph
G = nx.from_numpy_array(adj_matrix, create_using=nx.DiGraph())

# Plotting
def plot_graph(network, layout, title):
    pos = layout(network)
    pesos = [adj_matrix[i, j] for i, j in network.edges]
    node_sizes = [np.log10(2 ** g) + 4 for g in clasesG]
    
    nx.draw(network, pos, node_color=colorinchis, edge_color=pesos, width=pesos, 
            node_size=node_sizes, with_labels=True, font_color='white')
    plt.title(title)
    plt.show()

# Plot circle layout
plot_graph(G, nx.circular_layout, "Circle Layout")

# Plot Fruchterman-Reingold layout
plot_graph(G, nx.spring_layout, "Fruchterman-Reingold Layout")

# Plot LGL layout
plot_graph(G, nx.shell_layout, "Large Graph Layout")
