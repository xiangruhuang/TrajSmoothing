import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

def filter(minx, maxx, miny, maxy, x, y):
    new_x = []
    new_y = []
    flag = False
    for xi, yi in zip(x, y):
        if (xi < minx or xi > maxx):
            continue
        if (yi < miny or yi > maxy):
            continue
        flag = True
        break
    if flag:
        return x, y
    else:
        return new_x, new_y

input_file = sys.argv[1]
output_file = sys.argv[2]
method = 'scatter'
if 'graph' in input_file:
    method = 'plot'

if 'raw_points' in input_file:
    method = 'plot'

#if len(sys.argv) > 2 and sys.argv[2] == 'zoom':
#    zoom = True
#else:
#    zoom = False

zoom = True

if zoom:
    if 'graph' in input_file:
        zoom_maxx = 116.2870
        zoom_minx = 116.2865
        zoom_maxy = 39.9981
        zoom_miny = 39.9976
    else:
        #zoom_maxx = 116.2865
        #zoom_minx = 116.2815
        #zoom_maxy = 39.998
        #zoom_miny = 39.995
        zoom_maxx = 116.2872
        zoom_minx = 116.2860
        zoom_maxy = 39.9985
        zoom_miny = 39.9973

print 'Going to %s file: %s' % (method, input_file)

if len(sys.argv) > 3:
    title = sys.argv[3]
    plt.rc('text', usetex=True)
else:
    title = input_file

p2w = {}

zoomlines =[]
if 'dds_points' in input_file:
    """ load dds_edges"""
    with open('dds_edges', 'r') as fin:
        lines = fin.readlines()
        for line in lines:
            tokens = [float(token) for token in line.strip().split(' ')]
            x = tokens[0::2]
            y = tokens[1::2]
            zoomlines.append(x)
            zoomlines.append(y)

if 'ours_points' in input_file:
    """ load ours_edges"""
    zoomlines =[]
    with open('ours_edges', 'r') as fin:
        lines = fin.readlines()
        for line in lines:
            tokens = [float(token) for token in line.strip().split(' ')]
            x = tokens[0::2]
            y = tokens[1::2]
            zoomlines.append(x)
            zoomlines.append(y)

fig, ax = plt.subplots(figsize=(15,15))
if zoom:
    if 'graph' in input_file:
        axins = zoomed_inset_axes(ax, 24, loc=1)
    else:
        axins = zoomed_inset_axes(ax, 10, loc=1)
    for axis in ['top','bottom','left','right']:
        axins.spines[axis].set_linewidth(3)

if ('ours_graph' in input_file) or ('dds_graph' in input_file):
    maxx = -1e10
    minx = 1e10
    maxy = -1e10
    miny = 1e10
    draw = []
    scale = 1e-6
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
        ax.scatter(points[:, 0], points[:, 1], marker='o', s=2, color='r',
                facecolor='none')
        if zoom:
            axins.scatter(points[:, 0], points[:, 1], marker='o', s=2, color='r',
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
            
            c = 'r'
            if a % 3 == 0:
                ax.arrow(pa[0], pa[1], pb[0]-pa[0], pb[1]-pa[1], color=c,
                        linewidth=scale, width=scale)
                if zoom:
                    axins.arrow(pa[0], pa[1], pb[0]-pa[0], pb[1]-pa[1], color=c,
                            linewidth=scale, width=scale)
            else:
                ax.plot(x, y, color=c, linewidth=0.5)
                if zoom:
                    axins.plot(x, y, color=c, linewidth=2.0)
            
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
    if 'raw_points' in input_file:
        lw = 5e-1
        ms = 0.2
    else:
        lw = 1e-2
        ms = 0.2
    draw = []
    drawx = []
    drawy = []
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
            c = 'r'
            if method == 'scatter':
                #plt.scatter(x, y, color=c, s=0.5)
                drawx = drawx + list(x)
                drawy = drawy + list(y)
            else:
                #if count % 5 == 0:
                #    ax.arrow(x[0], y[0], x[1]-x[0], y[1]- y[0], lw=lw, width=lw,
                #            color='k')
                #    if zoom:
                #        axins.arrow(x[0], y[0], x[1]-x[0], y[1]- y[0], lw=lw,
                #                width=lw, color='k')
                #else:
                drawx = drawx + list(x)
                drawy = drawy + list(y)
                #if zoom:
                #x, y = filter(zoom_minx, zoom_maxx, zoom_miny, zoom_maxy, x, y)
                if count % 100 == 0:
                    draw.append(x)
                    draw.append(y)
                #plt.plot(x, y, color=c, marker='.', markersize=0.2, linewidth=lw)
    if method == 'scatter':
        ax.scatter(drawx, drawy, color='r', s=0.5)
        if zoom:
            axins.scatter(drawx, drawy, color='r', s=0.5)
            if len(zoomlines) > 0:
                axins.plot(lw=lw, color='r', *zoomlines)
    else:
        ax.plot(linewidth=lw, markersize=ms, color='r', *draw)
        #ax.plot(linewidth=lw, markersize=ms, color='k', *draw)
        if zoom:
            axins.plot(linewidth=lw, markersize=ms, color='r', *draw)
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
        axins.set_xlim([zoom_minx, zoom_maxx])
        axins.set_ylim([zoom_miny, zoom_maxy])
        mark_inset(ax, axins, loc1=4, loc2=2, fc="none", ec="0.0", lw=2)
    else:
        axins.set_xlim([zoom_minx, zoom_maxx])
        axins.set_ylim([zoom_miny, zoom_maxy])
        mark_inset(ax, axins, loc1=4, loc2=2, fc="none", ec="0.0", lw=2)

plt.xlabel('x', fontsize=0.01)
plt.ylabel('y', fontsize=0.01)
plt.xticks([], visible=False, fontsize=0.01)
plt.yticks([], visible=False, fontsize=0.01)
plt.tight_layout()
#plt.show()
plt.savefig(output_file)
