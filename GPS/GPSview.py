import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util import repeat

input_files = sys.argv

print input_files

plt.figure(figsize=(5, 5))

maxx = -1e10
minx = 1e10
maxy = -1e10
miny = 1e10
name_concat = ''
draw = []
for input_file in input_files[1:]:
    if input_file.startswith('-'):
        input_file = input_file[1:]
        method = 'plot'
    else:
        method = 'scatter'
    print input_file
    name_concat = name_concat + '_' + input_file
    with open(input_file, 'r') as fin:
        lines = fin.readlines()
        num_trace = len(lines) / repeat
        if num_trace == 0:
            num_trace = 1
        print '#trace=%d' % num_trace
        for count, line in enumerate(lines):
            tokens = np.asarray([float(s) for s in line.strip().split(' ')])
            x = tokens[0::2]
            maxx = max(maxx, max(x))
            minx = min(minx, min(x))
            
            y = tokens[1::2]
            maxy = max(maxy, max(y))
            miny = min(miny, min(y))
            
            c = cm.hot((count / repeat)*1.0/num_trace)
            c = 'k'
            for xi, yi in zip(x, y):
                draw.append([xi, yi])
            #if method == 'scatter':
            #    plt.scatter(x, y, s=0.02, color=c)
            #else:
            #    plt.plot(x, y, color=c, markersize=0.02,
            #            linewidth=0.01)

draw = np.array(draw)
if method == 'scatter':
    plt.scatter(draw[:, 0], draw[:, 1], s=0.02, color=c)
else:
    plt.plot(draw[:, 0], draw[:, 1], color=c, markersize=0.02, linewidth=0.01)

plt.xlim(minx, maxx)
plt.ylim(miny, maxy)
plt.title('Raw Data')
plt.xticks([], visible=False)
plt.yticks([], visible=False)
plt.axis('equal')
plt.show()
