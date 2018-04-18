import sys
import numpy as np
import matplotlib.pyplot as plt

input_file = sys.argv[1]
method = 'scatter'
if len(sys.argv) >= 3:
    method = sys.argv[2]

with open(input_file, 'r') as fin:
    lines = fin.readlines()
    maxx = -1e10
    minx = 1e10
    maxy = -1e10
    miny = 1e10
    for line in lines:
        tokens = np.asarray([float(s) for s in line.strip().split(' ')])
        x = tokens[0::2]
        maxx = max(maxx, max(x))
        minx = min(minx, min(x))
        
        y = tokens[1::2]
        maxy = max(maxy, max(y))
        miny = min(miny, min(y))
        
        if method == 'scatter':
            plt.scatter(x, y, s=0.01)
        else:
            plt.plot(x, y, 'k-+', markersize=0.01, linewidth=0.1)
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)
    plt.title(input_file)
    plt.show()
