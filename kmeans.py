from sklearn.cluster import KMeans

import numpy as np
import matplotlib.pyplot as plt
import os
import cPickle

#xmin = 116.0
#xmax = 116.4
#ymin= 39.6
#ymax = 40.0

n=1000

#filename = 'truncated_traj.txt'
filename = 'synthetic/trajectories.txt'
model_path = './synthetic/kmeans.save%d' % n

if not os.path.exists(model_path):
    xy = []
    with open(filename, 'r') as fin:
        lines = fin.readlines()
        for count, line in enumerate(lines):
            token = line.strip().split(' ')
            
            x = token[0::2]
            y = token[1::2]
            for x_i, y_i in zip(x, y):
                xy.append([x_i, y_i])

    kmeans = KMeans(n_clusters=n, random_state=0).fit(xy)

    with open(model_path, 'wb') as fout:
        cPickle.dump(kmeans, fout)
else:
    with open(model_path, 'rb') as fin:
        kmeans = cPickle.load(fin)

with open(filename, 'r') as fin:
    plt.figure()
    lines = fin.readlines()
    for count, line in enumerate(lines):
        token = line.strip().split(' ')
        
        x = token[0::2]
        y = token[1::2]

        a = []
        b = []
        for x_i, y_i in zip(x, y):
            label = kmeans.predict([[float(x_i), float(y_i)]])
            ab = kmeans.cluster_centers_[label][0]
            a.append(ab[0])
            b.append(ab[1])
        plt.plot(a, b, '.', markersize=0.1)
    #plt.xlim(xmin, xmax)
    #plt.ylim(ymin, ymax)
    plt.savefig('kmeans.eps')
