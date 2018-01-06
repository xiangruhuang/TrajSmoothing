import numpy
import itertools
import matplotlib.pyplot as plt
import sys

#points = numpy.random.randint(0, 100, (1000, 2))

X = 20
Y = 20

def generate_traces_on_grids(X, Y):

    points = [[a,b] for a, b in itertools.product(range(X), range(Y))]
    num_trace = 10
    p = 0.9
    dirx = [0, 1, 0, -1]
    diry = [1, 0, -1, 0]

    traces = []
    for t in range(num_trace):
        trace = []
        trace.append([numpy.random.randint(0, X-1), numpy.random.randint(0, Y-1)])
        while (numpy.random.uniform(0.0, 1.0) <= p) or (len(trace) < 20):
            x, y = trace[-1]
            dir_idx = numpy.random.randint(0, 3)
            new_x = x + dirx[dir_idx]
            new_y = y + diry[dir_idx]
            while not (new_x < X and new_x >= 0 and new_y < Y and new_y >= 0):
                dir_idx = numpy.random.randint(0, 3)
                new_x = x + dirx[dir_idx]
                new_y = y + diry[dir_idx]
            trace.append([new_x, new_y])
        trace = numpy.asarray(trace)
        traces.append(trace)
        #plt.plot(trace[:, 0], trace[:, 1])

    #plt.show()
    return traces

def create_samples(output_file, traces):
    with open(output_file, 'w') as fout:
        h = plt.figure()
        for t in range(100):
            for i, trace in enumerate(traces):
                T = len(trace)
                #print(T, trace.shape)
                time_stamps = numpy.random.uniform(0, T-1, size=(T, 1))
                #print(time_stamps.shape)
                time_stamps = sorted(time_stamps)
                sample = []
                for count, stamp in enumerate(time_stamps):
                    idx = int(numpy.trunc(stamp))
                    frac = stamp - idx
                    point = []
                    noise = numpy.random.uniform(-0.1, 0.1, size=(2))
                    point.append(frac*trace[idx, 0] + (1.0-frac)*trace[idx+1, 0] + noise[0])
                    point.append(frac*trace[idx, 1] + (1.0-frac)*trace[idx+1, 1] + noise[1])
                    if count != 0:
                        fout.write(' ')
                    fout.write('%f %f' % (point[0], point[1]))
                    sample.append(point)
                fout.write('\n')
                sample = numpy.asarray(sample)
                
                plt.plot(sample[:, 0], sample[:, 1], '.', markersize=0.1)
        plt.xlim(0, X)
        plt.ylim(0, Y)
        plt.savefig('sample.eps')

traces = generate_traces_on_grids(X, Y)

create_samples(sys.argv[1], traces)
