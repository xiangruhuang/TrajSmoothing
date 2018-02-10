from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import sys
import os
sys.path.append('../')
import BVH
import matplotlib.colors as colors

try:
    list_file = sys.argv[1]
except:
    print "python visualize.py <list_file>"
    exit(1)

def fullpath(path):
    if path.startswith('/'):
        return path
    cwd = os.getcwd()
    return cwd + '/' + path

with open(list_file, 'r') as fin:
    bvh_list = [fullpath(line.strip()) for line in fin.readlines() if
            line.strip().split('/')[-1] != 'rest.bvh']

motions = []
X = []
order = 3
used_joints = np.array([
    0,
    2, 3, 4,
    7, 8, 9,
    12, 13, 15, 16,
    18, 19, 20,
    25, 26, 27])
for count, f in enumerate(bvh_list):
    motion = BVH.load(f)
    motions.append(motion)
    print '%d / %d' % (count, len(bvh_list))
    el = motion.rotations.euler()
    el = el[:, used_joints, 0:2]
    el = np.reshape(el, [el.shape[0], -1])
    #X_c = []
    #c1 = 0.0
    #c2 = 0.0
    #for i in range(order, el.shape[0]-order):
    #    p_i = []
    #    for offset in range(-order, order):
    #        diff = el[i+offset+1] - el[i+offset]
    #        c1 += np.linalg.norm(diff[used_joints], 2) ** 2
    #        c2 += np.linalg.norm(diff, 2) ** 2
    #        diff = diff[used_joints]
    #        diff = diff / np.linalg.norm(diff, 2)
    #        p_i.append(diff)
    #    p_i = np.concatenate(p_i)
    #    p_i = p_i / np.linalg.norm(p_i)
    #    X_c.append(p_i)
    #print c1, c2
    #el = np.array(X_c)
    #el = el[:, [i for i in range(el.shape[1]) if i % 17 <= 2]]
    X.append(el)

X = np.concatenate(X)

print X.shape

X_embedded = TSNE(n_components=2, perplexity=300, verbose=1, angle=0.2).fit_transform(X)
print X_embedded.shape
fig = plt.figure()
ax = fig.add_subplot(111)#, projection='3d')

l = 0
r = 0
acolors = ['violet', 'steelblue', 'coral']
print acolors
offset = []
for count, motion in enumerate(motions):
    offset.append(l)
    l = r
    r = r + motion.rotations.euler().shape[0]
    ax.plot(X_embedded[l:r, 0], X_embedded[l:r, 1],
            color=acolors[count])

#with open('../smooth/matchings', 'r') as fin:
#    lines = fin.readlines()
#    for line in lines:
#        m1, f1, m2, f2, d = [float(token) for token in line.strip().split(' ')]
#        m1 = int(m1)
#        f1 = int(f1)
#        m2 = int(m2)
#        f2 = int(f2)
#        id1 = offset[m1] + f1
#        id2 = offset[m2] + f2
#        ax.plot(X_embedded[[id1, id2], 0], X_embedded[[id1, id2], 1],
#                X_embedded[[id1, id2], 2], 'r-',
#                linewidth=d)

#ax.scatter(X_embedded[:, 0], X_embedded[:, 1], X_embedded[:, 2], c=acolors)
plt.show()
