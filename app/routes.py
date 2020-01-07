from app import app, db
from flask import render_template, redirect, url_for, request
from app.forms import SpecificationForm
from app.models import GeneCircuit, Input_Output
import sys
from app.image_handling import decode_io, image_handling
from app.opt_topology import build_circuit
from app.draw_circuit import draw
from decimal import Decimal, getcontext
 
@app.route('/')
@app.route('/home')
def home():
    return render_template('home.html')

@app.route('/database')
def database_view():
    q = GeneCircuit.query.all()
    database_sum = []
    for c in q:
        database_sum.append((c.circuit_name,
                             c.gene_nodes,
                             c.input_nodes, 
                             c.output_nodes))
    return render_template('database_summary.html',
                           database_sum=database_sum)

@app.route('/circuit_summary/<circuit_name>')
def circuit_summary(circuit_name):
    q = GeneCircuit.query.all()
    for c in q:
        print(str(c), file=sys.stderr)
        if str(c.circuit_name) == circuit_name:
            print('trig', file=sys.stderr)
            top = c.topology.split()
            topology = eval(''.join([top[0]] + [s + ',' for s in top[1:-1]] + [top[-1]]))
            gene_nodes = c.gene_nodes
            input_nodes = c.input_nodes
            output_nodes = c.output_nodes
            break
        
    gene_labels = [chr(l) for l in range(65, 65 + int(gene_nodes))]
    inputs_labels = ['I' + str(l) for l in range(int(input_nodes))]
    for o in range(int(output_nodes)):
        gene_labels[o] = gene_labels[o] + '_out'
    node_labels = gene_labels + inputs_labels
    edges = [(v, n) for n in node_labels for v in node_labels]
    filtered_edges = [e for e in edges if e[1][0] != 'I']
    edges_with_coeffs = [e for e in zip(filtered_edges, topology[int(gene_nodes):])]
    pruned = [((le[0][0], le[0][1]), get_coeff(le[1])) for le in edges_with_coeffs if le[1] > 100. or le[1] < -100.]
        
    edge_starts = [e[0][0] for e in pruned]
    edge_ends = [e[0][1] for e in pruned]
    kms = [e[1] for e in pruned]
    
    return render_template('circuit_summary.html', 
                           circuit_name=circuit_name, 
                           hills=[h for h in topology[:int(gene_nodes)]], 
                           edge_starts=edge_starts, 
                           edge_ends=edge_ends, 
                           kms=kms, 
                           gene_labels=gene_labels)

@app.route('/start', methods=['GET', 'POST'])
def start():
    form = SpecificationForm()
    if form.validate_on_submit():
        circuit_name = form.circuit_name.data
        gene_nodes = form.gene_nodes.data
        output_nodes = form.output_nodes.data
        input_nodes = form.input_nodes.data
        max_expression = form.max_expression.data
        time_span = form.time_span.data
        circuit = GeneCircuit(circuit_name=circuit_name,
                              gene_nodes=gene_nodes,
                              output_nodes=output_nodes,
                              input_nodes=input_nodes,
                              max_expression=max_expression,
                              time_span=time_span)
        db.session.add(circuit)
        db.session.commit()
        return redirect(url_for('out_time_series', 
                                circuit_name=circuit_name, 
                                gene_nodes=gene_nodes, 
                                output_nodes=output_nodes, 
                                input_nodes=input_nodes,
                                max_expression=max_expression,
                                time_span=time_span))
    return render_template('start.html', form=form)

@app.route('/outputs/<circuit_name>/<gene_nodes>/<output_nodes>/<input_nodes>/<max_expression>/<time_span>', 
           methods=['GET', 'POST'])
def out_time_series(circuit_name, gene_nodes, output_nodes, input_nodes, max_expression, time_span):
    if request.method == 'GET':
        return render_template('time_series.html', 
                               circuit_name=circuit_name, 
                               gene_nodes=int(gene_nodes), 
                               output_nodes=int(output_nodes), 
                               input_nodes=int(input_nodes), 
                               out_iter=list(range(int(output_nodes))), 
                               in_iter=list(range(int(input_nodes))), 
                               max_expression=int(max_expression), 
                               time_span=int(time_span))
    if request.method == 'POST':
        data = request.form['saved_canvas']
        #print(data[22:], file=sys.stderr)
        with open('app/static/outputs/' + str(circuit_name), 'w') as outfile:
            outfile.write(data[22:])
        outfile.close()
        io = Input_Output(circuit=circuit_name)
        db.session.add(io)
        db.session.commit()
        decode_io(circuit_name, 'outputs')
        return 'output time series'

@app.route('/inputs/<circuit_name>/<gene_nodes>/<output_nodes>/<input_nodes>/<max_expression>/<time_span>', 
           methods=['GET', 'POST'])
def in_time_series(circuit_name, gene_nodes, output_nodes, input_nodes, max_expression, time_span):
    if request.method == 'GET':
        return render_template('in_time_series.html', 
                               circuit_name=circuit_name, 
                               gene_nodes=int(gene_nodes), 
                               output_nodes=int(output_nodes), 
                               input_nodes=int(input_nodes), 
                               out_iter=list(range(int(output_nodes))), 
                               in_iter=list(range(int(input_nodes))), 
                               max_expression=int(max_expression), 
                               time_span=int(time_span))
    if request.method == 'POST':
        data = request.form['saved_canvas']
        #print(data[22:], file=sys.stderr)
        with open('app/static/inputs/' + str(circuit_name), 'w') as outfile:
            outfile.write(data[22:])
        outfile.close()
        decode_io(circuit_name, 'inputs')
        return 'input time series'
        
@app.route('/check_build/<circuit_name>/<gene_nodes>/<output_nodes>/<input_nodes>/<max_expression>/<time_span>', methods=['GET'])
def check_build(circuit_name, gene_nodes, output_nodes, input_nodes, max_expression, time_span):
    if request.method == 'GET':
        image_handling(circuit_name, max_expression, time_span, output_nodes, 'outputs')
        image_handling(circuit_name, max_expression, time_span, input_nodes, 'inputs')
        return render_template('check_build.html',
                               circuit_name=circuit_name,
                               gene_nodes=gene_nodes,
                               output_nodes=int(output_nodes),
                               input_nodes=int(input_nodes),
                               max_expression=max_expression,
                               time_span=time_span)
    
@app.route('/building/<circuit_name>/<gene_nodes>/<output_nodes>/<input_nodes>/<max_expression>/<time_span>', methods=['GET', 'POST'])
def building(circuit_name, gene_nodes, output_nodes, input_nodes, max_expression, time_span):
    if request.method == 'GET':
        return render_template('building.html',
                               circuit_name=circuit_name,
                               gene_nodes=gene_nodes,
                               output_nodes=int(output_nodes),
                               input_nodes=int(input_nodes),
                               max_expression=max_expression,
                               time_span=time_span)
    
def get_coeff(coeff):
    getcontext().prec = 5
    if coeff > 100.:
        return float(Decimal(coeff)-Decimal(100.))
    elif coeff < -100.:
        return float(Decimal(coeff)+Decimal(100.))
    
@app.route('/circuit_build/<circuit_name>/<gene_nodes>/<output_nodes>/<input_nodes>/<max_expression>/<time_span>', methods=['GET'])
def circuit_build(circuit_name, gene_nodes, output_nodes, input_nodes, max_expression, time_span):
    if request.method == 'GET':
        circuit = build_circuit(circuit_name,
                                gene_nodes,
                                output_nodes,
                                input_nodes,
                                time_span)
        #draw(circuit_name, int(gene_nodes), int(output_nodes), int(input_nodes), circuit)
        query = GeneCircuit.query.filter_by(circuit_name=circuit_name).first()
        query.topology = str(circuit)
        db.session.commit()
        
        gene_labels = [chr(l) for l in range(65, 65 + int(gene_nodes))]
        inputs_labels = ['I' + str(l) for l in range(int(input_nodes))]
        for o in range(int(output_nodes)):
            gene_labels[o] = gene_labels[o] + '_out'
        node_labels = gene_labels + inputs_labels
        edges = [(v, n) for n in node_labels for v in node_labels]
        filtered_edges = [e for e in edges if e[1][0] != 'I']
        edges_with_coeffs = [e for e in zip(filtered_edges, circuit[int(gene_nodes):])]
        pruned = [((le[0][0], le[0][1]), get_coeff(le[1])) for le in edges_with_coeffs if le[1] > 100. or le[1] < -100.]
        
        edge_starts = [e[0][0] for e in pruned]
        edge_ends = [e[0][1] for e in pruned]
        kms = [e[1] for e in pruned]
    
        return render_template('circuit_build.html', 
                               circuit_name=circuit_name,
                               hills=[h for h in circuit[:int(gene_nodes)]],
                               edge_starts=edge_starts,
                               edge_ends=edge_ends,
                               kms=kms,
                               gene_labels=gene_labels)
    

