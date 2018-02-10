import sys
import BVH
import numpy as np
from util import *
import os

try:
    list_file = sys.argv[1]
    src_file = sys.argv[2]
except:
    print "python align.py <list_file> <source file> (nickname)"
    exit(1)

def fullpath(path):
    if path.startswith('/'):
        return path
    cwd = os.getcwd()
    return cwd + '/' + path

if len(sys.argv) > 3:
    nickname = sys.argv[3]
else:
    for i in range(10000):
        if os.path.exists('/home/xiangru/Projects/Qixing/TrajSmoothing/motion/videos/%d.mp4'):
            continue
        else:
            nickname = '/home/xiangru/Projects/Qixing/TrajSmoothing/motion/videos/%d' % i
            break

folder = nickname
with open(list_file, 'r') as fin:
    bvh_list = [fullpath(line.strip())+'.project' for line in fin.readlines() if
            line.strip().split('/')[-1] != 'rest.bvh']

print bvh_list

tgt_files = [f for f in bvh_list if not src_file in f]
print tgt_files
src_file = [f for f in bvh_list if not (f in tgt_files)][0]
print 'src=', src_file
#print 'tgts=', tgt_files

print "aligning %s" % src_file
""" Align other bvh files with this specific file"""

src, src_rotations = BVH.load(src_file, with_org_rotations=True)

#print src_rotations.shape

""" new implementation, with velocity """
tgts = []

align_rate = 0.0

for tnum, filename in enumerate(tgt_files):
    aligned_ids = []
    with open(filename.replace('.project', '.aligned'), 'r') as fin:
        aligned_ids = [int(line.strip()) for line in fin.readlines()]
    if len(aligned_ids) == 0:
        continue
    print "reading %s" % filename
    align_rate += len(aligned_ids)
    tgt = BVH.joint_load(filename, src, src_rotations, aligned_ids)
    tgts.append(tgt)

if len(tgts) > 0:
    align_rate /= len(tgts)
    print 'align_rate=%f' % align_rate

print 'ploting and saving as %s.mp4' % nickname

#plot([src]+tgts, output_video='%s.mp4' % nickname)
plot([src]+tgts, output_video=None)

""" old implementation, no velocity """
#tgts = []
#align = [[] for src_frame in range(src.rotations.shape[0])]
#
#act_tnum = 0
#
#for tnum, filename in enumerate(tgt_files):
#    print "reading %s ( %d / %d )" % (filename, tnum, len(tgt_files))
#    tgt = BVH.load(filename)
#    if len(tgt.rotations) == 0:
#        continue
#    tgts.append(tgt)
#    with open('%s.aligned' % filename, 'r') as fin:
#        lines = fin.readlines()
#        for fnum, line in enumerate(lines):
#            aligned_id = int(line.strip())
#            align[aligned_id].append((act_tnum, fnum))
#            #print np.linalg.norm(tgt.positions[fnum]- src.positions[aligned_id])
#    print "aligned = %s " % str(tgt.rotations.shape)
#    
#    act_tnum += 1
#
#align_rate = 0.0
#for i in range(len(align)):
#    if len(align[i]) > 0:
#        align_rate += len(align[i])
#align_rate /= len(align)

#print 'align_rate=%f' % align_rate

# aligned_plot(src, tgts, align)
