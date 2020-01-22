from pyeda.inter import *
import re
#from opt_topology import get_dimensions
from random import uniform

def mini_expression(tt, variables):
    x = ttvars('x', variables)
    o = truthtable(x, tt)
    return espresso_tts(o)

def extract_nots(expressions):
    nots = [set(re.findall('~x' + re.escape('[') + '\d+', exp)) for exp in expressions]
    return [i[-1] for inner in nots for i in inner]

def extract_ands(expressions):
    and_split = [e.split('And')[1:] for e in expressions]
    flat_and_split = [a for inner in and_split for a in inner]
    grab = lambda and_str: ''.join([c for c in and_str if c.isdigit() or c=='~' or c==','])
    extract = [grab(a) for a in flat_and_split]
    filter_trailing = lambda ex: ex[:-1] if ex[-1] == ',' else ex
    return [filter_trailing(a).split(',') for a in extract]

def extract_ors(ands, expressions):
    return [exp.count('And') for exp in expressions]

def logic_synth_position(nots, ands, ors, dimensions, inputs, outputs, genes):
    position = [0. for i in range(dimensions)]
    get_nodes = lambda o, g, i : ['o' + str(o_count) for o_count in range(outputs)] + ['g' + str(g_count) for g_count in range(genes)] + ['i' + str(i_count) for i_count in range(inputs)]
    nodes = get_nodes(outputs, genes, inputs)
    position_labels = ['hill o' + str(i) for i in range(outputs)] + ['hill g' + str(i) for i in range(genes)] + [i for inner in [[(i, j) for i in nodes] for j in nodes] for i in inner if i[1] not in ['i' + str(i_count) for i_count in range(inputs)]]
    gate_labels = ['' for i in range(len(nodes))]
    
    and_assigned = ''
    #add ors
    for or_gate in range(len(ors)):
        #add gate node
        or_out_node = 'o' + str(or_gate)
        node_idx = nodes.index(or_out_node)
        gate_labels[node_idx] = 'or gate' + str(or_gate)
        #add edges
        self_edge_idx = position_labels.index(('o' + str(or_gate), 'o' + str(or_gate)))
        position[self_edge_idx] = -295.
        gate_inputs = ors[or_gate]
        for i in range(gate_inputs):
            if and_assigned == '':
                node_idx = nodes.index('g0')
                gate_labels[node_idx] = 'and gate0'
                and_edge_idx = position_labels.index(('g0', or_out_node))
                position[and_edge_idx] = 200.
                and_assigned = 'g0'
            else:
                node_idx = nodes.index(and_assigned)
                gate_labels[node_idx] = 'and gate' + str(int(and_assigned[-1]) + 1)
                and_edge_idx = position_labels.index(('g' + str(int(and_assigned[-1]) + 1), or_out_node))
                position[and_edge_idx] = 200.
                and_assigned = 'g' + str(int(and_assigned[-1]) + 1)
    
    print(and_assigned)
    print(position)
    
    #add nots
    count = 0
    for i in range(int(and_assigned[-1]) + 1, genes):
        node_idx = nodes.index('g' + str(i))
        gate_labels[node_idx] = 'not gate' + nots[count]
        self_edge_idx = position_labels.index(('g' + str(i), 'g' + str(i)))
        position[self_edge_idx] = -299.5
        not_edge_idx = position_labels.index(('i' + nots[count], 'g' + str(i)))
        position[not_edge_idx] = -143.2
        count = count + 1
    
    #add ands
    for i in range(int(and_assigned[-1]) + 1):
        and_self_edge_idx = position_labels.index(('g' + str(i), 'g' + str(i)))
        position[and_self_edge_idx] = -294.5
        for and_input in ands[i]:
            if and_input[0] != '~':
                and_edge_idx = position_labels.index(('i' + and_input, 'g' + str(i)))
                position[and_edge_idx] = 140.
            else:
                not_gate = gate_labels.index('not gate' + and_input[1:])
                and_edge_idx = position_labels.index((nodes[not_gate], 'g' + str(i)))
                position[and_edge_idx] = 140.
    
    return position

def mutate(circuit_supergraph, genes, outputs):
    mod_circuit = []
    for h in range(genes + outputs):
        mod_circuit.append(uniform(1.0, 5.0))
    for e in range(genes + outputs, len(circuit_supergraph)):
        if circuit_supergraph[e] == 1.:
            if uniform(0, 1) > 0.5:
                mod_circuit.append(uniform(-300, 300))
            else:
                mod_circuit.append(1.)
        else:
            if uniform(0, 1) > 0.5:
                potential_km = circuit_supergraph[e] + uniform(-30, 30)
                if potential_km < 300. and potential_km > -300.:
                    mod_circuit.append(potential_km)
                else:
                    mod_circuit.append(circuit_supergraph[e])
            else:
                mod_circuit.append(circuit_supergraph[e])
    return mod_circuit

def starting_positions(best_guess, swarm_size, split, genes, outputs, inputs):
    from_logic_synth = int(swarm_size*split)
    from_rnd = swarm_size - from_logic_synth
    
    set_1s = lambda d: d if d != 0. else 1.
    logic_synth = [set_1s(d) for d in best_guess]
    logic_synth_starts = []
    for start in range(from_logic_synth - 1):
        add = mutate(logic_synth, genes, outputs)
        logic_synth_starts.append(add)
        
    rnd_starts = []
    #edge_count = get_dimensions(genes + outputs, inputs) - (genes + outputs)
    edge_count = 192
    for start in range(from_rnd):
        rnd_starts.append([uniform(1.0, 5.0) for h in range(genes + outputs)] + [uniform(-300., 300.) for e in range(edge_count)])
    
    return [logic_synth] + logic_synth_starts + rnd_starts
    
if __name__ == '__main__':
    tables = ['00010111', '01101001']
    expressions = [str(mini_expression(tt, 3)) for tt in tables]
    nots = extract_nots(expressions)
    print(sorted(nots))
    ands = extract_ands(expressions)
    print(ands)
    ors = extract_ors(ands, expressions)
    print(ors)
    #position = logic_synth_position(sorted(nots), ands, ors, get_dimensions(12, 3), 3, 2, 10)
    position = logic_synth_position(sorted(nots), ands, ors, 192, 3, 2, 10)
    
    print('################################################')
    print(position)
    print('################################################')
    
    sp = starting_positions(position, 100, 0.5, 10, 2, 3)
    
    lt = []
    print('################################################')
    for p in sp:
        print(p)
        lt.append(len(p))
    print('################################################')
    
    print(set(lt))
    
    print((len(sp)))
    
    print(sp[0])
    print(len(sp[0]))
    
    #check in bounds
    for p in sp:
        for dim in p:
            if dim >= 300. or dim <= -300.:
                print(p)
