from app import db

class GeneCircuit(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    circuit_name = db.Column(db.String(32), index=True)
    gene_nodes = db.Column(db.Integer, index=True)
    output_nodes = db.Column(db.Integer, index=True)
    input_nodes = db.Column(db.Integer, index=True)
    topology = db.Column(db.String(10000), index=True)
    max_expression = db.Column(db.Integer, index=True)
    time_span = db.Column(db.Integer, index=True)
    
class Input_Output(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    circuit = db.ForeignKey('gene_circuit.id')
    
    
