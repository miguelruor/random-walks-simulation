import numpy as np
import networkx as nx
from scipy.stats import expon
import os
import time as timing
from datetime import datetime

def infinitesimalGenerator(G, gamma):
  # gamma is the mutation rate 

  A = nx.laplacian_matrix(G)
  H = -gamma* A

  return H

def nextState(states, i, Q):
  # choose next state to jump to from the i-th row of the jump process transition matrix
  # states are the set of state of the random walk
  # i is the actual state
  # Q is the infinitesimal generator of the continuous random walk
  probs = [-Q[i, j]/Q[i, i] for j in states]
  probs[i] = 0 

  return np.random.choice(states, p = probs)

def simulation_cw(gspace, phenotypes, initial_genotype, max_simulation_time, gamma):
  start = timing.time()
  date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

  tau = {} # tau (hitting time) for every phenotype: time it takes quantum walk to find given phenotype
  N = {} # number of jumps the random walk takes to find a genotype with a new phenotype

  H = infinitesimalGenerator(gspace, gamma)

  jump = 0
  time = 0
  
  #initialization
  for phen in phenotypes:
    tau[phen] = -1
    N[phen] = -1

  actual_state = initial_genotype 

  # phenotypes of actual state
  phenotypes_actual_state = gspace.nodes[actual_state]['phenotypeName']

  # start of simulation
  while time < max_simulation_time:
    for phen in phenotypes_actual_state:
      if tau[phen] < 0: # update hitting times of novel phenotypes
        tau[phen] = time
        N[phen] = jump

    T = expon.rvs(scale = -1/H[actual_state, actual_state], size=1)[0] # holding time

    time += T
    jump += 1
    
    # choice next random state
    actual_state = nextState(list(gspace.nodes), actual_state, H)

    # phenotypes of actual state
    phenotypes_actual_state = gspace.nodes[actual_state]['phenotypeName']
  
  # end of simulation
  end = timing.time()

  simulation = {}

  simulation["initial_gen_index"] = initial_genotype
  simulation["initial_gen"] = gspace.nodes[initial_genotype]['sequence']
  simulation["initial_phen"] = gspace.nodes[initial_genotype]['phenotypeName'][0]
  simulation["transition_rate"] = gamma
  simulation["max_simulation_time"] = max_simulation_time
  simulation["total_mutations"] = jump
  simulation["computing_time"] = end-start
  simulation["simulation_time"] = time
  simulation["date"] = date

  for phen in phenotypes:
    simulation['tau_'+phen] = tau[phen] if tau[phen] >= 0 else time
    simulation['mutations_'+phen] = N[phen]

  return simulation

  