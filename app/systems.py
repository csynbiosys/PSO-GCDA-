import math

class Systems:
    
    """ODE system object representing a circuit design - created and evaluated for every particle in each iteration
    args:
        params - known/constant params + particle position in n+n^2+in dimensional space
            promoter leakiness, transcription rate, rna degredation rate, translation rate, protein degredation rate per gene - 5 * gene_count in length 
            hill coeff per gene - gene_count in length
            activation/repression coeffs - gene_count^2 + input_count * gene_count in length - coeff for each edge in gene circuit supergraph, +ve for activation, -ve for repression
        gene_count - number of genes in circuit design
        inputs - time series for input entities
    """

    def __init__(self, params, genes, inputs):
        self.gene_count = int(genes)
        self.params = params[:5*int(genes)]
        self.hills = params[5*int(genes):5*int(genes) + int(genes)]
        self.inputs = inputs
        #nested list of activation/repression coeffs - each inner list contains coeffs for in edges of gene at same index as that inner list - length gene_count + input_count
        self.reg_coeffs = [params[5*int(genes) + int(genes):][i:i+self.gene_count+len(self.inputs)] for i in range(0, len(params[5*int(genes) + int(genes):]), self.gene_count+len(self.inputs))]
    
    #rate of change for protein entity at position pos
    def protein_ode(self, protein, rna, pos): 
        return (self.params[1 + 5*pos] * rna) - (self.params[3 + 5*pos] * protein)
    
    #calculate sum of hill function values for all in edges to a particular node (representing a gene entity at position pos) in gene circuit super graph 
    def regulation(self, entities, pos):
        entity_reg_coeffs = self.reg_coeffs[pos]
        running = 1.
        for count, coeff in enumerate(entity_reg_coeffs):
            if coeff < -100.:
                running = running * (((abs(coeff)-100)**self.hills[pos])/((abs(coeff)-100)**self.hills[pos] + entities[count]**self.hills[pos]))
            elif coeff > 100.:
                running = running * ((entities[count]**self.hills[pos])/((abs(coeff)-100)**self.hills[pos] + entities[count]**self.hills[pos]))
            else:
                running = running * 1.
        return running
    
    #rate of change for rna entity at position pos
    def rna_ode(self, entities_and_inputs, entities, pos):
        return (self.params[5*pos] + self.params[2 + 5*pos] * self.regulation(entities_and_inputs, pos)) - (self.params[4 + 5*pos] * entities[pos + self.gene_count])
    
    #used by odeint to evaluate ODE model, quantities of each entity and time point in and rate of change of each entity out
    def sys(self, y, t):
        constraint = lambda i: 0. if y[i] < 0. else y[i]
        entities = [constraint(i) for i in range(len(y))]
        if t>999.:
            time = 999
        else:
            time = int(t)
        #proteins + inputs
        entities_and_inputs = entities[:int(len(entities)/2)] + [input_entity[time] for input_entity in self.inputs]
        rna_odes = [self.rna_ode(entities_and_inputs, entities, i) for i in range(len(entities[len(entities)//2:]))]
        protein_odes = [self.protein_ode(entities[i], entities[i + self.gene_count], i) for i in range(len(entities[:len(entities)//2]))]
        return [ode for ode in protein_odes + rna_odes]
