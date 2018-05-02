import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

input_file = sys.argv[1]
output_file = sys.argv[2]
method = 'scatter'
if 'graph' in input_file:
    method = 'plot'

#if len(sys.argv) > 2 and sys.argv[2] == 'zoom':
#    zoom = True
#else:
#    zoom = False

zoom = True
    

print 'Going to %s file: %s' % (method, input_file)

if len(sys.argv) > 3:
    title = sys.argv[3]
    plt.rc('text', usetex=True)
else:
    title = input_file

p2w = {}

fig, ax = plt.subplots(figsize=(15,15))
if zoom:
    if 'graph' in input_file:
        axins = zoomed_inset_axes(ax, 24, loc=1)
    else:
        axins = zoomed_inset_axes(ax, 3.1, loc=1)

if ('ours_graph' in input_file) or ('dds_graph' in input_file):
    maxx = -1e10
    minx = 1e10
    maxy = -1e10
    miny = 1e10
    draw = []
    scale = 3e-7
    with open(input_file, 'r') as fin:
        lines = fin.readlines()[1:]
        num_point = 0
        for count, line in enumerate(lines):
            if line.startswith('Edges'):
                num_point = count
                break
        
        """ Points: """
        points = []
        for c, line in enumerate(lines[:num_point]):
            p = [float(token) for token in line.strip().split(' ')]
            points.append(p)
        points = np.array(points)
        print points.shape
        ax.scatter(points[:, 0], points[:, 1], marker='o', s=2, color='k',
                facecolor='none')
        if zoom:
            axins.scatter(points[:, 0], points[:, 1], marker='o', s=2, color='k',
                    facecolor='none')
        
        """ Edges: """
        for count, line in enumerate(lines[num_point+1:]):
            a, b, w = [float(s) for s in line.strip().split(' ')]
            a = int(a)
            b = int(b)
            if w >= -1:
                p2w[(a,b)] = w

        for key, val in p2w.items():
            a, b = key
            op_val = p2w.get((b,a), 0.0)
            #if val < op_val * 1.0 / r:
            #    continue
            pa = points[a, :]
            pb = points[b, :]

            x = [pa[0], pb[0]]
            maxx = max(maxx, max(x))
            minx = min(minx, min(x))
            
            y = [pa[1], pb[1]]
            maxy = max(maxy, max(y))
            miny = min(miny, min(y))
            
            c = 'k'
            if a % 3 == 0:
                ax.arrow(pa[0], pa[1], pb[0]-pa[0], pb[1]-pa[1], color=c,
                        linewidth=scale, width=scale)
                if zoom:
                    axins.arrow(pa[0], pa[1], pb[0]-pa[0], pb[1]-pa[1], color=c,
                            linewidth=scale, width=scale)
            else:
                ax.plot(x, y, color=c, linewidth=0.1)
                if zoom:
                    axins.plot(x, y, color=c, linewidth=0.1)
            
    #ax.plot(color=c, linewidth=scale, *draw)
    #ax.arrow(color=c, linewidth=1.0, *draw)
    #if zoom:
    #    axins.plot(color=c, linewidth=scale, *draw)
    print minx, maxx, miny, maxy
    #ax.axis('equal')
    #axins.axis('equal')
    #plt.xticks(fontsize=0.01)
    #plt.yticks(fontsize=0.01)
    #plt.xlabel('x', fontsize=0.01)
    #plt.ylabel('y', fontsize=0.01)
    ##plt.title(, fontsize=40)
    #plt.xlim((minx, maxx))
    #plt.ylim((miny, maxy))
    #ax.axis('off')
    #plt.tight_layout(pad=1.0)
    #plt.axes().set_aspect('auto')
    #plt.show()
    #plt.savefig(output_file) 
else:
    if input_file == 'raw_graph':
        lw = 3.0
    else:
        lw = 0.1
    ms = 0.2
    draw = []
    if method == 'scatter':
        drawx = []
        drawy = []
    with open(input_file, 'r') as fin:
        lines = fin.readlines()
        maxx = -1e5
        minx = 1e5
        maxy = -1e5
        miny = 1e5
        for count, line in enumerate(lines):
            tokens = np.asarray([float(s) for s in line.strip().split(' ')])
            x = tokens[0::2]
            maxx = max(maxx, max(x))
            minx = min(minx, min(x))
            
            y = tokens[1::2]
            maxy = max(maxy, max(y))
            miny = min(miny, min(y))
            c = 'k'
            if method == 'scatter':
                #plt.scatter(x, y, color=c, s=0.5)
                drawx = drawx + list(x)
                drawy = drawy + list(y)
            else:
                draw.append(x)
                draw.append(y)
                #plt.plot(x, y, color=c, marker='.', markersize=0.2, linewidth=lw)
    if method == 'scatter':
        ax.scatter(drawx, drawy, color='k', s=0.5)
        if zoom:
            axins.scatter(drawx, drawy, color='k', s=0.5)
    else:
        #ax.plot(linewidth=lw, markersize=ms, color='k', *draw)
        ax.arrow(linewidth=lw, markersize=ms, color='b', *draw)
        if zoom:
            axins.plot(linewidth=lw, markersize=ms, color='k', *draw)
    #print minx, maxx, miny, maxy
    #plt.xticks([], visible=False, fontsize=0.01)
    #plt.yticks([], visible=False, fontsize=0.01)
    #plt.tight_layout(pad=0.0)
    #plt.show()
    #plt.savefig(output_file)

#ax.scatter(draw[:, 0], draw[:, 1], s=0.007, color='k')
ax.axis('equal')
if zoom:
    axins.axis('equal')
ax.axis('off')
ax.set_xlim([minx, maxx])
ax.set_ylim([miny, maxy])

if zoom:
    if 'graph' in input_file:
        axins.set_xlim([116.2865, 116.2870])
        axins.set_ylim([39.9976, 39.9981])
        mark_inset(ax, axins, loc1=4, loc2=2, fc="none", ec="0.0")
    else:
        axins.set_xlim([116.2815, 116.2865])
        axins.set_ylim([39.995, 39.998])
        mark_inset(ax, axins, loc1=4, loc2=2, fc="none", ec="0.0")

plt.xlabel('x', fontsize=0.01)
plt.ylabel('y', fontsize=0.01)
plt.xticks([], visible=False, fontsize=0.01)
plt.yticks([], visible=False, fontsize=0.01)
plt.tight_layout()
#plt.show()
plt.savefig(output_file)
