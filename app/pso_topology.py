from pyswarms.single.global_best import GlobalBestPSO
from pyswarms.utils.plotters import plot_cost_history
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from app.systems import Systems
from scipy.integrate import odeint
from itertools import product
from random import uniform
import numpy as np
import sys

class PSO_topology:
    
    """gene circuit design automation particle swarm object
    args:
        particles - number of particles in swarm
        dimensions - number of dimensions for particle positions
        options - hyperparam dict
            cognitive, social and inertia params
        bounds - bounds tuple
            list containing lower bounds, list containing upper bounds
        loops - number of iterations to run optimization
        init_cond - starting quantities for each entity
        time - number of time points to run ODE models over
        target - time series of target entity
        inputs - time series of input entities
        genes - number of genes to use in circuit design
        param_known - constant parameters
            promoter leakiness, transcription rate, rna degredation rate, translation rate, protein degredation rate per gene
            5 * gene_count in length
        processes - number of threads over which to run swarm in parallel
        param_penalty - penalty factor for edge elimination
            added to mean squared error for each edge in gene circuit graph
    """
    
    def __init__(self, particles, dimensions, options, bounds, loops, init_cond, time, target, inputs, genes, param_known, processes, param_penalty):
        self.particles = particles
        self.dimensions = dimensions
        self.options = options
        self.bounds = bounds
        self.loops = loops
        self.init_cond = init_cond
        self.time = time
        self.target = target
        self.inputs = inputs
        self.genes = genes
        self.param_known = param_known
        self.processes = processes
        self.param_penalty = param_penalty
    
    #decide initial particle positions such that each particle starts representing a different circuit topology
    def initial_particle_positions(self):
        all_tops = [i for i in product([-110., 0., 110.], repeat=(self.dimensions - self.genes))]
        all_tops.sort(key=lambda top: 1 - len(set(top)))
        init_tops = all_tops[:self.particles]
        init_pos = [[uniform(self.bounds[0][0] + 0.1, self.bounds[1][0] - 0.1) for i in range(self.genes)] + list(top) for top in init_tops]
        return np.array(init_pos)
    
    #run ODE model corresponding to single particle and calculate corresponding cost
    def particle_cost(self, particle):
        model = Systems(self.param_known + list(particle), self.genes, self.inputs)
        result_sys = odeint(model.sys, self.init_cond, self.time)
        num_active_edges = len([e for e in particle if e > 100. or e < -100.])
        return sum([mean_squared_error(self.target[t], result_sys[:, t]) for t in range(len(self.target))]) + (self.param_penalty * num_active_edges)
    
    #calculate cost for each particle in swarm and serialize positions
    def swarm_cost(self, particle_swarm):
        f = open('positions.txt', 'a')
        f.write(str(particle_swarm) + '\n')
        f.close()
        return [self.particle_cost(p) for p in particle_swarm]
    
    #carry out particle swarm optimization and plot cost of swarm best per iteration
    def swarm(self):
        circuit_swarm = GlobalBestPSO(n_particles=self.particles, dimensions=self.dimensions, options=self.options, bounds=self.bounds)
        cost, pos = circuit_swarm.optimize(self.swarm_cost, iters=self.loops, n_processes=self.processes)
        plot_cost_history(circuit_swarm.cost_history)
        plt.show()
        return cost, pos
    
        
    
