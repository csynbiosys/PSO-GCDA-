from app.pso_topology import PSO_topology
from app.systems import Systems
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sys
from app.rates import Rates
from pickle import load

"""calculate number of dimensions for particle positions
(hill coeff per gene) + ((number of genes)^2 activation/repression coeffs) + (activation/repression coeff per input per gene)
dimensions = n+n^2+i*n
"""
def get_dimensions(gene_count, input_count):
    return ((int(gene_count)) + (int(gene_count)**2) + (int(input_count)*int(gene_count)))

#run ODE model for designed circuit
def run_sys(param_position, gene_count, inputs, init):
    init_cond_sys = init
    time_sys = list(range(1000))
    model = Systems(param_position, gene_count, inputs)
    return odeint(model.sys, init_cond_sys, time_sys)

#plot time series for each in entity in ODE model of designed circuit
def plot(results, out_file, target, gene_count):
    for i in range(gene_count*2):    
        fig = plt.figure(figsize=(10, 10))    
        res = results[:, i]
        plt.plot(list(range(1000)), res)
        if i in range(len(target)):
            plt.plot(list(range(1000)), target[i])
        plt.savefig(out_file + str(i) + '.png')

#draw gene circuit graph for designed circuit        
def draw_graph(position, genes):
    dc = draw_circuit.Draw_Circuit(position, genes)
    dc.build_circuit_graph()
    dc.draw_circuit_graph()

#read pickle generated from image to time series code - used to get targets and inputs    
def get_time_series(circuit_name, nodes, folder):
    time_series = []
    for t in range(int(nodes)):
        f = open('app/static/' + folder +'/extracted/' + circuit_name + str(t) + '.pickle', 'rb')
        time_series.append(load(f))
        f.close()
    return [[max(0, v) for v in ts] for ts in time_series]
    
def build_circuit(circuit_name, gene_nodes, output_nodes, input_nodes, time_span):
    targets = get_time_series(circuit_name, output_nodes, 'outputs')
    inputs = get_time_series(circuit_name, input_nodes, 'inputs')
    
    #object for rates and bounds
    rate = Rates(int(gene_nodes), get_dimensions(gene_nodes, input_nodes))
    
    """gene circuit design automation particle swarm object
    args:
        0 - number of particles in swarm
        1 - number of dimensions for particle positions
        2 - hyperparam dict
            cognitive, social and inertia params
        3 - bounds tuple
            list containing lower bounds, list containing upper bounds
        4 - number of iterations to run optimization
        5 - starting quantities for each entity
        6 - number of time points to run ODE models over
        7 - time series of target entity
        8 - time series of input entities
        9 - number of genes to use in circuit design
        10 - constant parameters
            promoter leakiness, transcription rate, rna degredation rate, translation rate, protein degredation rate per gene
            5 * gene_count in length
        11 - number of threads over which to run swarm in parallel
        12 - penalty factor for edge elimination
            added to mean squared error for each edge in gene circuit graph
    """
    opt = PSO_topology(1, 
                       get_dimensions(gene_nodes, input_nodes), 
                       {'c1': 1.3, 'c2': 2.7, 'w':0.3}, 
                       rate.get_bounds(), 
                       1, 
                       rate.get_init(), 
                       list(range(int(time_span))), 
                       targets, 
                       inputs, 
                       int(gene_nodes), 
                       rate.get_constant_params(), 
                       1, 
                       2750.)
    
    #best cost found by swarm and corresponding position in n + n^2 + i*n dimensional space
    cost, pos = opt.swarm()
    
    #run ODE model for designed circuit and generate output figures
    res = run_sys(rate.get_constant_params() + list(pos), gene_nodes, inputs, rate.get_init())
    plot(res, 'app/static/builds/' + circuit_name, targets, int(gene_nodes))
    
    return pos

