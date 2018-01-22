from util import *

import sys
import numpy as np

from util import animation_plot
from util import process_file
import scipy.io as io

if len(sys.argv) < 5:
    print "python view_align.py <a in list> <b in list> <threshold> <matchings> "
    exit(1)

with open('walk_remove_odd/lists.txt', 'r') as fin:
    filelist = [line.strip().split('/')[-1] for line in fin.readlines()]

a = int(sys.argv[1])
b = int(sys.argv[2])
threhold = float(sys.argv[3])

print 'aligning %s and %s with threshold %f' % (filelist[a], filelist[b],
        threshold)

with open(sys.argv[4], 'r') as fin:
    lines = fin.readlines()
    pairs = lines[0::2]
    matchings = lines[1::2]
    for pair, match in zip(pairs, matchings):
        i, j, T = [int(token) for token in pair.strip().split(' ')]
        if T == 0:
            continue
        tuples = [int(token) for count, token in
                enumerate(match.strip().split(' ')) if (count % 3 != 2)]
        assert T*2 == len(tuples)
        
        if not (i == a and j == b):
            continue
        for a, b in zip(tuples[0::2], tuples[1::2]):
            compare(i, j, a+3, b+3)

clips = []

clips += process_file(sys.argv[1])
clips += process_file(sys.argv[2])

clips = np.array(clips)
np.savez_compressed('data/compare', clips=clips)
io.savemat('data/compare.mat', dict(clips=clips), do_compression=True)

database = np.load('data/compare.npz')['clips']

animation_plot([database[0], database[1]], H36=False, repeat=True)

