import igraph as ig
import matplotlib.pyplot as plt
from decimal import Decimal, getcontext
    
def get_coeff(coeff):
    getcontext().prec = 5
    if coeff > 100.:
        return float(Decimal(coeff)-Decimal(100.))
    elif coeff < -100.:
        return float(Decimal(coeff)+Decimal(100.))
    
def draw(circuit_name, gene_nodes, output_nodes, input_nodes, circuit):
    #get nodes and edges
    gene_labels = [chr(l) for l in range(65, 65 + int(gene_nodes))]
    inputs_labels = ['I' + str(l) for l in range(input_nodes)]
    for o in range(output_nodes):
        gene_labels[o] = gene_labels[o] + '_out'
    node_labels = gene_labels + inputs_labels
    edges = [(v, n) for n in node_labels for v in node_labels]
    filtered_edges = [e for e in edges if e[1][0] != 'I']
    edges_with_coeffs = [e for e in zip(filtered_edges, circuit[gene_nodes:])]
    pruned = [((le[0][0], le[0][1]), get_coeff(le[1])) for le in edges_with_coeffs if le[1] > 100. or le[1] < -100.]
    
    circuit_g = ig.Graph(directed=True)
    
    for node in node_labels:
        circuit_g.add_vertices(node)
    circuit_g.vs['hill'] = [h[1] + ' hill: ' + str(h[0]) for h in zip(circuit[:gene_nodes]+['n/a' for i in range(len(inputs_labels))], node_labels)]
    for edge, label in pruned:
        circuit_g.add_edges([(edge[0], edge[1])])
    circuit_g.es['coeff'] = [l[1] for l in pruned]
    
    layout = circuit_g.layout('circle')
    visual_style = {}
    visual_style['vertex_label'] = circuit_g.vs['hill']
    visual_style['vertex_size'] = 150
    visual_style['vertex_color'] = 'green'
    visual_style['edge_label'] = circuit_g.es['coeff']
    visual_style['bbox'] = (1000, 1000)
    visual_style['margin'] = 200
    out = ig.plot(circuit_g, layout=layout, **visual_style)
    
    out.save(fname='/static/builds/' + circuit_name + '_viz.png')
    
  
