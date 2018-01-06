from util import *
import matplotlib.pyplot as plt
import numpy

with open('matches.txt', 'r') as fin:
    lines = fin.readlines()
    for count, line in enumerate(lines):
        if count % 3 == 0:
            trace_i = numpy.asarray([float(x) for x in line.strip().split(' ')])
            x_i = trace_i[0::2]
            y_i = trace_i[1::2]
        elif count % 3 == 1:
            trace_j = numpy.asarray([float(x) for x in line.strip().split(' ')])
            x_j = trace_j[0::2]
            y_j = trace_j[1::2]
        else:
            if len(line) < 2:
                continue
            match = numpy.asarray([int(x) for x in line.strip().split(' ')])
            if len(match) == 0:
                continue
            print match.shape
            assert match.shape[0] % 2 == 0
            match_i = match[0::2]
            match_j = match[1::2]
            plt.plot(x_i, y_i)
            plt.plot(x_j, y_j)
            for i, j in zip(match_i, match_j):
                plt.plot([x_i[i], x_j[j]], [y_i[i], y_j[j]], '--')
            plt.xlim(0, 10)
            plt.ylim(0, 10)
            plt.show()
