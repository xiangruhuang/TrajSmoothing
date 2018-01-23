import sys
import BVH
import numpy as np
from util import *

try:
    folder = sys.argv[1]
    I = int(sys.argv[2])
except:
    print "python align.py <folder> <id>"
    print "By default assuming <folder>/list.txt exists"
    print "Aligning file with name shown in <id>-th line in <folder>/list.txt"
    exit(1)

try:
    with open(folder+'/list.txt', 'r') as fin:
        filelist = [line.strip() for line in fin.readlines()]
except:
    print "Error: file %s/list.txt doesn't exists" % folder
    print "Call `python bvh2traj.py %s`" % folder
    exit(1)

print "aligning %s" % filelist[I]
""" Align other bvh files with this specific file"""

rest = [f for count, f in enumerate(filelist) if I != count]

src = BVH.load(filelist[I])

#print type(src)

#plot([src])

threshold = 1.0

tgts = []
src_euler = src.rotations.euler()
align = [[] for src_frame in src_euler]
min_dist = [threshold for src_frame in src_euler]
for tnum, filename in enumerate(rest):
    tgt = BVH.load(filename)
    tgts.append(tgt)
    """ Aligning tgt and src """
    tgt_euler = tgt.rotations.euler()
    for i, src_frame in enumerate(src_euler):
        #print i, src_frame.shape
        min_d = 1e10
        min_index = -1
        dists = np.linalg.norm(tgt_euler - src_frame, 2, axis=(1,2))
        for j, tgt_frame in enumerate(tgt_euler):
            dist = dists[j]
            if dist < min_d:
                min_d = dist
                min_index = j
        if min_d < threshold:
            align[i].append((tnum, j))

align_rate = 0.0
for i in range(len(align)):
    if len(align[i]) > 0:
        align_rate += len(align[i])
align_rate /= len(align)

print 'align_rate=%f' % align_rate

aligned_plot(src, tgts, align)
