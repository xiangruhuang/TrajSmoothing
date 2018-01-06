import sys
import matplotlib.pyplot as plt
import os
from sklearn.cluster import KMeans
import cPickle
import numpy as np

#xmin = 116.0
#xmax = 116.8
#ymin= 39.6
#ymax = 40.2

def plot_from_file(input_file):
    with open(input_file, 'r') as fin:
        lines = fin.readlines()
        for line in lines:
            tokens = np.asarray([float(s) for s in line.strip().split(' ')])
            x = tokens[0::2]
            y = tokens[1::2]
            plt.plot(x, y, '.', markersize=0.1)
        #plt.show()
        plt.xlim(-0.5, 20.5)
        plt.ylim(-0.5, 20.5)
        plt.title(input_file)
        output_file = '/'.join(input_file.split('/')[:-1] +
                [input_file.split('/')[-1].split('.')[0]+'.png'])
        print output_file
        plt.savefig(output_file)

def read_traj(input_file):
    trace = []
    with open(input_file, 'r') as fin:
        lines = fin.readlines()
        for line in lines:
            tokens = np.asarray([float(s) for s in line.strip().split(' ')])
            x = tokens[0::2]
            y = tokens[1::2]
            xy = [[x_i, y_i] for x_i, y_i in zip(x, y)]
            trace.append(xy)

def train_kmeans(n, input_file):

    if not os.path.exists('./kmeans.save%d' % n):
        with open(input_file, 'r') as fin:
            xy = []
            lines = fin.readlines()
            for count, line in enumerate(lines):
                token = line.strip().split(' ')
                x = token[0::2]
                y = token[1::2]
                for x_i, y_i in zip(x, y):
                    xy.append([x_i, y_i])
        print len(xy)
        kmeans = KMeans(n_clusters=n, random_state=0, n_jobs=-1).fit(xy)
        with open('kmeans.save%d' % n, 'wb') as fout:
            cPickle.dump(kmeans, fout)
    else:
        with open('kmeans.save%d' % n, 'rb') as fin:
            kmeans = cPickle.load(fin)
    return kmeans

def predict_kmeans(kmeans, input_file):
    with open(input_file, 'r') as fin:
        lines = fin.readlines()
        ab = None
        for count, line in enumerate(lines):
            token = line.strip().split(' ')
            print count
            x = token[0::2]
            y = token[1::2]
            xy = []
            for x_i, y_i in zip(x, y):
                xy.append([float(x_i), float(y_i)])
            labels = kmeans.predict(xy)
            ab = kmeans.cluster_centers_[labels]
            plt.plot(ab[:, 0], ab[:, 1])
        #plt.xlim(xmin, xmax)
        #plt.ylim(ymin, ymax)
        plt.show()

def plot_traj(input_file, output_file = None, xmin = 116.0, xmax = 116.4, ymin = 39.6, ymax = 40.0):
    if output_file is not None:
        fout = open(output_file, 'w')
    else:
        fout = None
    with open(input_file, 'r') as fin:
        lines = fin.readlines()
        for count, line in enumerate(lines):
            token = line.strip().split(' ')
            x = token[1::2]
            y = token[0::2]
            if len(x) <= 10 or len(y) <= 10:
                continue
            x = [float(x_i) for x_i in x]
            y = [float(y_i) for y_i in y]
            plt.plot(x, y, '+')
            if fout is not None:
                _count = 0
                for x_i, y_i in zip(x, y):
                    if _count != 0:
                        fout.write(' ')
                    fout.write('%f %f' % (x_i, y_i))
                    _count += 1
                fout.write('\n')
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.show()
    if fout is not None:
        fout.close()

def in_range(x, y, xmin, xmax, ymin, ymax):
    if (x <= xmax) and (x >= xmin):
        if (y <= ymax) and (y >= ymin):
            return True
    return False

def truncate_traj(input_file, output_file, xmin, xmax, ymin, ymax):
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        lines = fin.readlines()
        total_count = len(lines)
        plot_count = 0
        for count, line in enumerate(lines):
            token = line.strip().split(' ')
            print('%d/%d' % (count, total_count))
            
            x = token[0::2]
            y = token[1::2]
            if len(x) <= 10:
                continue
            x = [float(x_i) for x_i in x]
            y = [float(y_i) for y_i in y]
            
            workx = []
            worky = []
            for x_i, y_i in zip(x, y):
                if in_range(x_i, y_i, xmin, xmax, ymin, ymax):
                    workx.append(x_i)
                    worky.append(y_i)
                else:
                    if len(workx) > 4:
                        plot_count += 1
                        plt.plot(workx, worky, '-+')
                        _count = 0
                        for x_i, y_i in zip(workx, worky):
                            if _count != 0:
                                fout.write(' ')
                            fout.write('%f %f' % (x_i, y_i))
                            _count += 1
                        fout.write('\n')
                    workx = []
                    worky = []
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            if plot_count >= 10:
                break
            
        plt.show()
