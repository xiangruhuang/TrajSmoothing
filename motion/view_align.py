from util import *

import sys
import numpy as np
import BVH
import scipy.io as io
import random

try:
    filelist_path = sys.argv[1]
    matching_file = sys.argv[2]
except:
    print "python view_align.py <filelist_path> <matching_file>"
    exit(1)

motions = []

with open(filelist_path, 'r') as fin:
    filelist = [line.strip() for line in fin.readlines()]
    motions = [BVH.load(file) for file in filelist]

animations = []
offsets = []

grid_size = 1

with open(matching_file, 'r') as fin:
    matching = []
    count = 0
    lines = fin.readlines()
    random.seed(816)
    random.shuffle(lines)
    for line in lines:
        tokens = line.strip().split(' ')
        m1,f1,m2,f2 = [int(token) for token in tokens[:-1]]
        dist = float(tokens[-1])
        if dist < float(sys.argv[3]) or dist > float(sys.argv[4]):
            continue
        print dist
        a1 = motions[m1].copy()
        a2 = motions[m2].copy()
        l = 100
        ll = min(f1, f2, l)
        rr = min(len(a1.rotations)-1-f1, len(a2.rotations)-1-f2, l)
        a1.rotations = a1.rotations[f1-ll:f1+rr]
        a1.positions = a1.positions[f1-ll:f1+rr]
        a2.rotations = a2.rotations[f2-ll:f2+rr]
        a2.positions = a2.positions[f2-ll:f2+rr]
        animations.append(a1)
        animations.append(a2)

        c = count % (grid_size ** 2)
        r = 15
        offsets.append([r * (c / grid_size - grid_size / 2 + 0), r * (c %
            grid_size - grid_size / 2 + 2)])
        offsets.append([r * (c / grid_size - grid_size / 2 + 0), r * (c %
            grid_size - grid_size / 2 + 2)])
        if (count + 1) %  grid_size ** 2 == 0:
            clips = []

            for anim in animations:
                clips += process_animation(anim)

            clips = np.array(clips)
            output_video = 'videos/view_align%d.mp4' % count
            animation_plot_with_offset([clips[i] for i in range(len(animations))],
                    output_video = None, _offsets = offsets, repeat=True,
                    no_root_motion=True)
            animations = []
            offsets = []
        count += 1

#        
#
#a = int(sys.argv[1])
#b = int(sys.argv[2])
#threhold = float(sys.argv[3])
#
#print 'aligning %s and %s with threshold %f' % (filelist[a], filelist[b],
#        threshold)
#
#with open(sys.argv[4], 'r') as fin:
#    lines = fin.readlines()
#    pairs = lines[0::2]
#    matchings = lines[1::2]
#    for pair, match in zip(pairs, matchings):
#        i, j, T = [int(token) for token in pair.strip().split(' ')]
#        if T == 0:
#            continue
#        tuples = [int(token) for count, token in
#                enumerate(match.strip().split(' ')) if (count % 3 != 2)]
#        assert T*2 == len(tuples)
#        
#        if not (i == a and j == b):
#            continue
#        for a, b in zip(tuples[0::2], tuples[1::2]):
#            compare(i, j, a+3, b+3)
#
#clips = []
#
#clips += process_file(sys.argv[1])
#clips += process_file(sys.argv[2])
#
#clips = np.array(clips)
#np.savez_compressed('data/compare', clips=clips)
#io.savemat('data/compare.mat', dict(clips=clips), do_compression=True)
#
#database = np.load('data/compare.npz')['clips']
#
#animation_plot([database[0], database[1]], H36=False, repeat=True)

