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

plot([src])

#for filename in rest:
#    tgt = BVH.load(filename)
#    """ Aligning tgt and src """
#    print dir(src.rotations.interpolate)
#    #print src.rotations.euler()/np.pi * 180.0
#    break
#    #for src_rot in src.rotations:
#        #print src_rot.shape
