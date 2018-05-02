import numpy
import itertools
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
from util import repeat

#points = numpy.random.randint(0, 100, (1000, 2))

def generate_colors(total):

    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.Normalize(vmin=-20, vmax=10)
    cmap = cm.hot

    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = []
    for i in range(total):
        x = i * 1.0 / total
        colors.append(cm.hot(x))
    return colors

X = 10
Y = 10

def generate_traces_on_grids(X, Y):

    points = [[a,b] for a, b in itertools.product(range(X), range(Y))]
    num_trace = 1
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

def generate_traces_on_streets(X, Y):
    street_width = 0.01
    street_length = 1.0-street_width * 2.0
    mx = numpy.asarray([1, -1, -1, 1]) * street_width
    my = numpy.asarray([1,  1, -1, -1]) * street_width
    dx = numpy.asarray([0, -1, 0 , 1]) * street_length
    dy = numpy.asarray([1,  0, -1, 0]) * street_length
    new_ori = [3, 0, 1, 2]
    candidate_ori = [0, 1, 2, 3] 

    p = 0.8
    num_trace = 10
    traces = []
    for t in range(num_trace):
        x = numpy.random.randint(4, X-4)
        y = numpy.random.randint(4, Y-4)
        trace = []
        ori = 0
        #print t
        while ((numpy.random.uniform(0.0, 1.0) <= p) or (len(trace) < 20)):
            old_ori = ori

            ori = numpy.random.randint(0, len(candidate_ori))
            ori = candidate_ori[ori]

            while abs(old_ori - ori) == 2:
                ori = numpy.random.randint(0, len(candidate_ori))
                ori = candidate_ori[ori]
                
            #print x, y, ori
            x0 = x
            y0 = y
            x = x + mx[ori]
            y = y + my[ori]
            trace.append([x, y, ori])
            x = x + dx[ori]
            y = y + dy[ori]
            #trace.append([x, y, ori])
            x = x - mx[new_ori[ori]]
            y = y - my[new_ori[ori]]
            if ((x < 0.5) or (x > X-0.5) or (y < 0.5) or (y > Y-0.5)):
                trace = trace[:-1]
                #x, y, ori = trace[-1]
                x = x0
                y = y0
                #x = x - mx[new_ori[ori]]
                #y = y - my[new_ori[ori]]
        if len(trace) < 20:
             continue
        trace = numpy.asarray(trace)
        trace = trace[:, :-1]
        #print trace
        
        traces.append(trace)
        #plt.plot(trace[:, 0], trace[:, 1], 'k-', linewidth=0.01)
        #plt.xlim(0, X)
        #plt.ylim(0, Y)
    #plt.show()
    
    print '#trace=', len(traces)

    return traces

def generate_traces_on_circle_async(X, Y, r):
    angles = numpy.random.uniform(0.0, 2.0*numpy.pi, size=(repeat))
    noise_level = numpy.exp(-numpy.random.uniform(1.0, 1.0, size=(100)))
    points = numpy.asarray([[r*numpy.sin(angle), r*numpy.cos(angle)] for angle in angles])
    num_trace = repeat
    traces = []
    for t in range(num_trace):
        trace = []
        angle = angles[t]
        T = numpy.random.randint(10, 100)
        for i in range(T):
            angle += 0.02*numpy.pi
            p = [r * numpy.sin(angle), r*numpy.cos(angle)] 
            trace.append(p)
        trace = numpy.asarray(trace)
        noise = numpy.random.randn(len(trace), 2) * 0.03
        trace = trace + noise
        traces.append(trace)
        plt.plot(trace[:, 0], trace[:, 1], 'r-')

    plt.show()

    with open('circle.txt', 'w') as fout:
        for trace in traces:
            for i in range(len(trace)):
                if i != 0:
                    fout.write(' ')
                fout.write('%f %f' % (trace[i, 0], trace[i, 1]))
            fout.write('\n')
    
    return traces

def generate_traces_on_line(T):
    traces = []
    sigma = 0.001
    start_points = numpy.random.uniform(0.0, 2*sigma, size=(repeat))
    for t in range(repeat):
        p = start_points[t]
        trace = []
        for i in range(T):
            pi = [p + i * sigma, 0.0]
            trace.append(pi)
        trace = numpy.asarray(trace)
        noise = numpy.random.uniform(-sigma, sigma, size=(len(trace), 2))
        trace = trace + noise
        traces.append(trace)
        plt.plot(trace[:, 0], trace[:, 1], 'r-')

    plt.show()

    with open('line.txt', 'w') as fout:
        for trace in traces:
            for i in range(len(trace)):
                if i != 0:
                    fout.write(' ')
                fout.write('%f %f' % (trace[i, 0], trace[i, 1]))
            fout.write('\n')
    
    return traces
     

def generate_traces_on_circle(X, Y, r, noisy=False):
    angles = numpy.linspace(0.0, 3.0*numpy.pi, 100)
    noise_level = numpy.exp(-numpy.random.uniform(1.0, 5.0, size=(100)))
    points = numpy.asarray([[r*numpy.sin(angle), r*numpy.cos(angle)] for angle in angles])
    num_trace = repeat
    traces = []
    for t in range(num_trace):
        trace = []
        if not noisy:
            a = numpy.random.randint(0, 99)
            b = numpy.random.randint(0, 99)
            while (a + 10 >= b):
                a = numpy.random.randint(0, 99)
                b = numpy.random.randint(0, 99)
            noise = numpy.random.randn(b-a+1, 2)*0.1
            trace = points[a:b+1, :]
            trace = trace + noise
        else:
            a = numpy.random.randint(0, 99)
            b = numpy.random.randint(0, 99)
            while (a + 10 >= b):
                a = numpy.random.randint(0, 99)
                b = numpy.random.randint(0, 99)
            noise = numpy.random.uniform(-1.0, 1.0, size=(b-a+1, 2))
            for i in range(a, b+1):
                noise[i - a, :] = noise[i - a, :] * noise_level[i]
            #noise = numpy.random.randn(b-a+1, 2)*0.1
            trace = points[a:b+1, :]
            trace = trace + noise
        traces.append(trace)
        plt.plot(trace[:, 0], trace[:, 1], 'r-')

    plt.show()

    #with open('circle.txt', 'w') as fout:
    #    for trace in traces:
    #        for i in range(len(trace)):
    #            if i != 0:
    #                fout.write(' ')
    #            fout.write('%f %f' % (trace[i, 0], trace[i, 1]))
    #        fout.write('\n')
   
    with open('circle.gt', 'w') as fout:
        points = []
        N = 200
        for t in range(N+1):
            angle = t * 1.0 / N * 2.0 * numpy.pi
            x = r * numpy.sin(angle)
            y = r * numpy.cos(angle)
            if t != 0:
                fout.write(' ')
            fout.write('%f %f' % (x, y))
        fout.write('\n')
        
    return traces

def read_traces(filename):
    traces = []
    with open(filename, 'r') as fin:
        lines = fin.readlines()
        for line in lines:
            numbers = [float(token) for token in line.strip().split(' ')]
            x = numbers[0::2]
            y = numbers[1::2]
            trace = [[xi, yi] for xi, yi in zip(x, y)]
            trace = numpy.array(trace)
            traces.append(trace)
    return traces

def create_samples(output_file, traces, repeat, sigma):
    colors = generate_colors(len(traces))
    with open(output_file+'.txt', 'w') as fout:
        h = plt.figure()
        for i, trace in enumerate(traces):
            for t in range(repeat):
                T = len(trace)
                #print(T, trace.shape)
                time_stamps = numpy.random.uniform(0, (T-1), size=(T / 20))
                #time_stamps = [i + 1e-5 for i in range(T - 1)]
                #print(time_stamps.shape)
                time_stamps = sorted(time_stamps)
                sample = []
                for count, stamp in enumerate(time_stamps):
                    idx = int(numpy.trunc(stamp))
                    frac = stamp - idx
                    point = []
                    noise = numpy.random.uniform(-sigma, sigma, size=(2))
                    px = (1.0-frac)*trace[idx, 0] + frac*trace[idx+1, 0] + noise[0]
                    py = (1.0-frac)*trace[idx, 1] + frac*trace[idx+1, 1] + noise[1]
                    point = [px, py]
                    if count != 0:
                        fout.write(' ')
                    fout.write('%f %f' % (point[0], point[1]))
                    sample.append(point)
                fout.write('\n')
                sample = numpy.asarray(sample)
                plt.plot(sample[:, 0], sample[:, 1], markersize=1,
                        color=colors[i])
        plt.savefig(output_file + '.png')

if sys.argv[1] == 'grids':
    traces = generate_traces_on_grids(X, Y)
    create_samples('grids', traces, 1000, 0.1)
elif sys.argv[1] == 'circle':
    traces = generate_traces_on_circle(0.0, 0.0, 1.0)
    #create_samples('circle', traces, 1, 0.0)
elif sys.argv[1] == 'street':
    traces = generate_traces_on_streets(X, Y)
    create_samples('street', traces, 100, 0.05)
elif sys.argv[1] == 'line':
    traces = generate_traces_on_line(10)
else:
    """ read traces """
    traces = read_traces(sys.argv[1])
    create_samples('read', traces, repeat, 0.05)
