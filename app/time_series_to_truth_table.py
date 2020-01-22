import math
import random
import numpy as np
from scipy.ndimage import gaussian_filter1d 
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.optimize import curve_fit
from statistics import median

def partitions(time_series):
    f = gaussian_filter1d([i[1] for i in time_series], 1.5)
    grads = np.gradient(f)
    tps = [(i, grads[i]) for i in range(1, len(grads)) if (grads[i] > 0. and grads[i-1] < 0.) or (grads[i] < 0. and grads[i-1] > 0.)]
    return [tp[0] for tp in tps]

def tt_rows(inputs, all_time_series):
    expected_partitions = 2**inputs - 1
    p = [partitions(ts) for ts in all_time_series]
    p_flat = [[i] for inner in p for i in inner]
    kmc = KMeans(n_clusters=expected_partitions).fit(p_flat)
    return [int(round(i)) for inner in kmc.cluster_centers_ for i in inner]

def tt_values(tt_rows, all_time_series, max_expression):
    tt_regions = [0] + sorted(tt_rows) + ['end']
    tt_regions = [(tt_regions[i-1], tt_regions[i]) for i in range(1, len(tt_regions))]
    exponential = lambda t, a, b, c: a * np.exp(b * t) + c
    on_off = lambda t: 0 if median(t) < 0.5*max_expression else 1
    tt = []
    for ts in all_time_series:
        col = []
        for row in tt_regions:
            to_fit = []
            if row[1] != 'end':
                to_fit = [i[1] for i in ts[row[0]:row[1]]]
            else:
                to_fit = [i[1] for i in ts[row[0]:]]
            col.append(on_off(to_fit))
        tt.append(col)
    return tt

if __name__ == '__main__':
    #test0 = [math.sin(50*i) + random.uniform(0.01, 0.3) for i in range(100)]
    #test1 = [math.sin(50*i+5) + random.uniform(0.01, 0.3) for i in range(100)]
    #test2 = [math.sin(50*i-5) + random.uniform(0.01, 0.3) for i in range(100)]
    #test3 = [math.sin(50*i+10) + random.uniform(0.01, 0.3) for i in range(100)]
    test0 = [0. + random.uniform(0., 100.) for i in range(500)] + [2000. + random.uniform(-100., 100.) for i in range(500)]
    test1 = [0. + random.uniform(0., 100.) for i in range(250)] + [2000. + random.uniform(-100., 100.) for i in range(250)] + [0. + random.uniform(0., 100.) for i in range(250)] + [2000. + random.uniform(-100., 100.) for i in range(250)]
    test2 = [0. + random.uniform(0., 100.) for i in range(125)] + [2000. + random.uniform(-100., 100.) for i in range(125)] + [0. + random.uniform(0., 100.) for i in range(125)] + [2000. + random.uniform(-100., 100.) for i in range(125)] + [0. + random.uniform(0., 100.) for i in range(125)] + [2000. + random.uniform(-100., 100.) for i in range(125)] + [0. + random.uniform(0., 100.) for i in range(125)] + [2000. + random.uniform(-100., 100.) for i in range(125)] 
    test3 = [0. + random.uniform(0., 100.) for i in range(375)] + [2000. + random.uniform(-100., 100.) for i in range(125)] + [0. + random.uniform(0., 100.) for i in range(125)] + [2000. + random.uniform(-100., 100.) for i in range(375)]
    
    test0 = list(enumerate(test0))
    test1 = list(enumerate(test1))
    test2 = list(enumerate(test2))
    test3 = list(enumerate(test3))
    
    p = [partitions(ts) for ts in [test0, test1, test2, test3]]
    
    pf = tt_rows(3, [test0, test1, test2, test3])
    
    tt = tt_values(pf, [test0, test1, test2, test3], 2000.)
    
    print(tt)
    
    plt.figure()
    plt.subplot(411)
    plt.plot([t[0] for t in test0], [t[1] for t in test0])
    for part in pf:
        plt.axvline(part, color='r')
        
    plt.subplot(412)
    plt.plot([t[0] for t in test1], [t[1] for t in test1])
    for part in pf:
        plt.axvline(part, color='r')
        
    plt.subplot(413)
    plt.plot([t[0] for t in test2], [t[1] for t in test2])
    for part in pf:
        plt.axvline(part, color='r')
        
    plt.subplot(414)
    plt.plot([t[0] for t in test3], [t[1] for t in test3])
    for part in pf:
        plt.axvline(part, color='r')
        
    plt.show()
