import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util import repeat
import os

false_cor = 0
average_w = 0.0

records = []
if os.path.isfile(sys.argv[1]):
    with open(sys.argv[1], 'r') as fin:
        lines = fin.readlines()
        for line in lines:
            tokens = line.strip().split(' ')
            i = int(tokens[0])
            fi = int(tokens[1])
            j = int(tokens[2])
            fj = int(tokens[3])
            w = float(tokens[4])
            if i / repeat != j / repeat:
                false_cor += 1
                average_w += w
            record = [i,fi,j,fj,w]
            records.append(record)

print 'matching file=', sys.argv[1]
print '#false correspondences= %d / %d' % (false_cor, len(records))
if false_cor > 0:
    print 'average false weight=%f' % (average_w / false_cor)

input_file = sys.argv[2]
print 'input file=', input_file

maxx = -1e10
minx = 1e10
maxy = -1e10
miny = 1e10

if input_file.startswith('-'):
    input_file = input_file[1:]
    method = 'plot'
else:
    method = 'scatter'
traces = []
with open(input_file, 'r') as fin:
    lines = fin.readlines()
    if len(lines) % repeat == 0:
        num_trace = len(lines) / repeat
    else:
        repeat = len(lines)
        num_trace = 1
    print '#trace=%d' % num_trace
    n = len(lines)
    for count, line in enumerate(lines):
        tokens = np.asarray([float(s) for s in line.strip().split(' ')])
        x = tokens[0::2]
        maxx = max(maxx, max(x))
        minx = min(minx, min(x))
        
        y = tokens[1::2]
        maxy = max(maxy, max(y))
        miny = min(miny, min(y))
        trace = np.array([[xi, yi] for xi, yi in zip(x, y)])
        traces.append(trace)

        #c = cm.hot((count / 100)*1.0/num_trace)
        #if method == 'scatter':
        #    plt.scatter(x, y, s=0.5, color=c)
        #else:
        #    plt.plot(x, y, color=c, marker='o', markersize=0.2, linewidth=0.1)

if num_trace > 1:
    plot_i = np.random.randint(repeat*2, n)
else:
    plot_i = np.random.randint(0, n)

for count, trace in enumerate(traces):
    if count != plot_i:
        c = cm.hot((count / repeat)*1.0 / num_trace)
        plt.scatter(trace[:, 0], trace[:, 1], s=1.0, color=c)
    else:
        c = cm.hot((count / repeat)*1.0 / (num_trace))
        plt.plot(trace[:, 0], trace[:, 1], color=c, marker='o', markersize=1.0, linewidth=1.0)

#maxx = 2
#minx = -2
#maxy = 2
#miny = -2
drawlines = []
count = 0
for I in range(repeat*2, n):
    i = plot_i
    for record in records:
        if i == int(record[0]):
            fi = int(record[1])
            j = int(record[2])
            fj = int(record[3])
            w = float(record[4])
            #if (j / repeat == i / repeat) or w < 0.5:
            #    continue
            count+=1
            p1 = traces[i][fi, :]
            p2 = traces[j][fj, :]
            p = np.array([p1, p2])
            drawlines.append(p[:, 0])
            drawlines.append(p[:, 1])
            if j / repeat == i / repeat:
                drawlines.append('g-')
            else:
                drawlines.append('g-')
            #plt.plot(p[:, 0], p[:, 1], 'g-', linewidth=0.5)
    if count > 0:
        break
#drawlines.append(linewidth=0.5)
plt.plot(linewidth=0.5, *drawlines)
plt.axis('equal')
plt.xlim(minx, maxx)
plt.ylim(miny, maxy)
plt.savefig('%s.png' % input_file)
#plt.show()
