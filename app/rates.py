
class Rates:
    
    def __init__(self, gene_count, dimensions):
        self.gene_count = gene_count
        self.dimensions = dimensions
        
    # boundaries of hill coeff, activation/repression coeff per gene
    def get_bounds(self):
        bounds = ([1. for i in range(self.gene_count)] + [-300. for i in range(self.dimensions - self.gene_count)], [5. for i in range(self.gene_count)] + [300. for i in range(self.dimensions  - self.gene_count)])
        return bounds

    # promoter leakiness, transcription rate, rna degredation rate, translation rate, protein degredation rate per gene
    def get_constant_params(self):
        constant_params = [0.03, 6.93, 29.97, 0.069, 0.347]*self.gene_count
        return constant_params

    # initial quantities of entities
    def get_init(self):
        init = [0. for i in range(self.gene_count)] + [float(i) for i in range(self.gene_count)]
        return init
