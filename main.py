import classic_simulation
import json
import os
import networkx as nx

if __name__ == "__main__":
    
    data_file  = open("data.json", "r")
    data = json.load(data_file)
    data_file.close()

    # parameters
    names = ['four_phen', 'musculus', 'neurospora', 'arabidopsis']
    gspace_name = names[int(os.environ['GSPACE'])]
    gspace = nx.read_gml(data[gspace_name]['g-space'], label='id')

    phenotypes_file = open(data[gspace_name]["phenotypes"], "r")
    phenotypes = phenotypes_file.read().splitlines()
    phenotypes_file.close()

    genotype_networks_file = open(data[gspace_name]["g-networks"], "r")
    genotype_networks = json.load(genotype_networks_file)
    genotype_networks_file.close()

    initial_genotype = 200
    max_simulation_time = 10**(8)
    gamma_c = 10**(-7)
    
    classic_simulation.simulation_cw(gspace, gspace_name, phenotypes, initial_genotype, max_simulation_time, gamma_c)