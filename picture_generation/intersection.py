import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull
from scipy.interpolate import interp1d

def distance2set(x, y, a, b):
    min_dist = 1e10
    for xi, yi in zip(x, y):
        if (xi-a) ** 2+ (yi-b) ** 2 < min_dist:
            min_dist = (xi-a) ** 2+ (yi-b) ** 2
    return min_dist

def interpolate(x, y, prec=1e-3):
    p = zip(x, y)
    #print type(p), len(p), len(p[0])
    p = np.array(p)
    new_p = [p[0]]
    for pi in p[1:]:
        last_pi = new_p[-1]
        num_inter = int(np.linalg.norm(pi - last_pi, 2) / prec) + 1
        for n in range(num_inter):
            ratio = (n + 1.0) / (num_inter + 1.0)
            next_pi = last_pi * (1.0-ratio) + pi * ratio
            new_p.append(next_pi)
        new_p.append(pi)
    new_p =  np.array(new_p)
    return new_p[:, 0], new_p[:, 1]

def fat(x0, y0):
    eps = 2.5e-2
    points = []
    x, y = interpolate(x0, y0, 3e-3)
    for xi, yi in zip(x, y):
        for t in np.linspace(0, np.pi*2.0, 200):
            points.append([xi + eps*np.cos(t), yi + eps * np.sin(t)])
        
    #xs = np.linspace(min(x)-eps, max(x) + eps, 1000)
    #ys = np.linspace(min(y)-eps, max(y) + eps, 1000)
    #points = []
    #for xi in xs:
    #    for yi in ys:
    #        if distance2set(x, y, xi, yi) < eps:
    #            points.append([xi, yi])
    points = np.array(points)
    #plt.scatter(points[:, 0], points[:, 1], s=0.01, facecolors='r', edgecolors='none')
    #hull = ConvexHull(points)
    #poly = []
    #for simplex in hull.simplices:
    #    poly.append([points[simplex, 0], points[simplex, 1]])
        #print simplex
        #plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
    #poly = np.array(poly)
    #print poly.shape
    return points

y1 = np.linspace(-2.0, 1.5, 200)
x1 = np.linspace(0.0, 0.0, 200)
plt.plot(x1, y1, 'b-')
points1 = fat(x1, y1)

y2 = np.linspace(0.0, 0.0, 200)
x2 = np.linspace(-2.0, 1.0, 200)
plt.plot(x2, y2, 'b-', label='Road')
points2 = fat(x2, y2)
points12 = np.concatenate((points1, points2), axis=0)
#plt.scatter(points12[:, 0], points12[:, 1], s=0.01, facecolors='b',
#        edgecolors='none')

with open('bezier.out', 'r') as fin:
    lines = fin.readlines()
    x3 = []
    y3 = []
    for line in lines:
        xi, yi = [float(token) for token in line.strip().split(' ')]
        x3.append(xi)
        y3.append(yi)
    x3 = np.array(x3)
    y3 = np.array(y3)
    x3 = x3 - x3[0]
    ratio = 2.0/(np.max(x3) - np.min(x3))
    x3 = x3 * ratio
    y3 = y3 - y3[0]
    ratio = 2.0/(np.max(y3) - np.min(y3))
    y3 = y3 * ratio
    y3 = y3 - y3[-1]
    
    points3 = fat(x3, y3)
    #plt.scatter(points3[:, 0], points3[:, 1], s=0.01, facecolors='r',
    #        edgecolors='none')
    plt.plot(x3, y3, 'r-', label='Road')

plt.scatter([0, 0, x3[-1]], [0, y3[0], 0], label='Intersections',
        facecolors='none',
        edgecolors='k', marker='o', s=80)

plt.scatter([0, 0.6544], [0.3385, 0.0], label='Fake Intersections', facecolors='none',
        edgecolors='g', marker='v', s=80)

lgd = plt.legend(loc='lower left', fontsize=18)
lgd.legendHandles[2]._sizes = [50]
lgd.legendHandles[3]._sizes = [50]
plt.axis('equal')
plt.axis('off')
plt.show()
