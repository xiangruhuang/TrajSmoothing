import sys
from util import *
import matplotlib.pyplot as plt
import numpy

input_file = sys.argv[1]
trajs = read_traj(input_file)
acc = 0.0
down = 0.0

with open(sys.argv[2], 'r') as fin:
    lines = fin.readlines()
    for line in lines:
        tokens = line.strip().split(' ')
        i = int(tokens[0])
        fi = int(tokens[1])
        j = int(tokens[2])
        fj = int(tokens[3])
        w = float(tokens[4])
        ti = numpy.asarray(trajs[i])
        tj = numpy.asarray(trajs[j])
        g = 5
        #if (ti[fi, 0] <= 1.0) and (tj[fj, 0] <= 1.0):
        #    acc += 1.0
        #elif (ti[fi, 0] >= 1.0) and (tj[fj, 0] >= 1.0):
        #    acc += 1.0
        #else:
        si = max(fi-g, 0)
        sj = max(fj-g, 0)
        ei = min(fi+g, len(ti))
        ej = min(fj+g, len(tj))
        plt.plot(ti[si:ei, 0], ti[si:ei, 1], 'b-')
        plt.plot(tj[sj:ej, 0], tj[sj:ej, 1], 'r-')
        plt.scatter(ti[fi, 0], ti[fi, 1], s=10)
        plt.scatter(tj[fj, 0], tj[fj, 1], s=10)
        plt.title(str(w))
        plt.show()
