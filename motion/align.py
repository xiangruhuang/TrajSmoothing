import sys
import BVH
import numpy as np
from util import *
import os

try:
    list_file = sys.argv[1]
except:
    print "python align.py <list_file> (nickname)"
    exit(1)

if len(sys.argv) > 2:
    nickname = sys.argv[2]
else:
    for i in range(10000):
        if os.path.exists('/home/xiangru/Projects/Qixing/TrajSmoothing/motion/videos/%d.mp4'):
            continue
        else:
            nickname = '/home/xiangru/Projects/Qixing/TrajSmoothing/motion/videos/%d' % i
            break

folder = nickname.split('/')[-1]
with open(list_file, 'r') as fin:
    bvh_list = ['/'+line.strip() for line in fin.readlines() if
            line.strip().split('/')[-1] != 'rest.bvh']

tgt_files = [f for f in bvh_list if os.path.exists(f+'.aligned')]
src_file = [f for f in bvh_list if not (f in tgt_files)][0]
print 'src=', src_file
#print 'tgts=', tgt_files

print "aligning %s" % src_file
""" Align other bvh files with this specific file"""

src, src_rotations = BVH.load(src_file, with_org_rotations=True)

print src_rotations.shape

""" new implementation, with velocity """
#tgts = []
#
#align_rate = 0.0
#
#for tnum, filename in enumerate(tgt_files):
#    aligned_ids = []
#    with open('%s.aligned' % filename, 'r') as fin:
#        aligned_ids = [int(line.strip()) for line in fin.readlines()]
#    if len(aligned_ids) == 0:
#        continue
#    print "reading %s" % filename
#    align_rate += len(aligned_ids)
#    tgt = BVH.joint_load(filename, src, src_rotations, aligned_ids)
#    tgts.append(tgt)
#
#align_rate /= len(tgts)
#
#print 'align_rate=%f' % align_rate
#
#print 'ploting and saving as %s.mp4' % nickname
#
##plot([src]+tgts, output_video='%s.mp4' % nickname)
#plot([src]+tgts, output_video=None)

""" old implementation, no velocity """
tgts = []
align = [[] for src_frame in range(src.rotations.shape[0])]

act_tnum = 0

for tnum, filename in enumerate(tgt_files):
    print "reading %s" % filename
    tgt = BVH.load(filename)
    if len(tgt.rotations) == 0:
        continue
    tgts.append(tgt)
    with open('%s.aligned' % filename, 'r') as fin:
        lines = fin.readlines()
        for fnum, line in enumerate(lines):
            aligned_id = int(line.strip())
            align[aligned_id].append((act_tnum, fnum))
            #print np.linalg.norm(tgt.positions[fnum]- src.positions[aligned_id])
    
    act_tnum += 1

align_rate = 0.0
for i in range(len(align)):
    if len(align[i]) > 0:
        align_rate += len(align[i])
align_rate /= len(align)

print 'align_rate=%f' % align_rate

aligned_plot(src, tgts, align)
