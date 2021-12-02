import classic_simulation
import json
import networkx as nx
import pandas as pd

if __name__ == "__main__":
    
    data_file  = open("data.json", "r")
    data = json.load(data_file)
    data_file.close()

    parameters_file = open("classical-parameters.json", "r")
    classical_parameters = json.load(parameters_file)
    parameters_file.close()

    # parameters
    names = ['four_phen', 'musculus', 'neurospora', 'arabidopsis']

    job_index = 0
    num_simulations = 1000

    if job_index < 10000:
        gspace_name = names[1]
        initial_node_index = job_index//1000
    elif job_index < 20000:
        gspace_name = names[2]
        initial_node_index = (job_index-10000)//1000
    elif job_index < 30000:
        gspace_name = names[3]
        initial_node_index = (job_index-20000)//1000
    else:
        gspace_name = names[0]
        initial_node_index = (job_index-30000)//1000

    gspace = nx.read_gml("data/"+data[gspace_name]['g-space'], label='id')

    phenotypes_file = open("data/"+data[gspace_name]["phenotypes"])
    phenotypes = phenotypes_file.read().splitlines()
    phenotypes_file.close()

    genotype_networks_file = open("data/"+data[gspace_name]["g-networks"])
    genotype_networks = json.load(genotype_networks_file)
    genotype_networks_file.close()

    initial_nodes = classical_parameters['initial_nodes'][gspace_name]

    initial_genotype = initial_nodes[initial_node_index]
    max_simulation_time = classical_parameters['max_time'] 
    gamma_c = classical_parameters['transition_rate']

    df = pd.DataFrame()
    
    for i in range(num_simulations):
        print(f"Simulation {i}")
        simulation = classic_simulation.simulation_cw(gspace, phenotypes, initial_genotype, max_simulation_time, gamma_c)
        df = df.append(simulation, ignore_index=True)

    df_past = pd.read_csv(f"results/{gspace_name}-cw.csv")

    df = pd.concat([df_past, df], ignore_index=True)
    df.to_csv(f"results/{gspace_name}-cw.csv", index=False)