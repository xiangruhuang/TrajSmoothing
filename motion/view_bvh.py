import sys
import numpy as np

from util import animation_plot
from util import process_file
import scipy.io as io

if len(sys.argv) < 2:
    print "python <bvh file> (start) (end)"
    exit(1)

clips = []

clips += process_file(sys.argv[1])

clips = np.array(clips)
np.savez_compressed('/home/xiangru/Projects/Qixing/TrajSmoothing/motion/data/tmp', clips=clips)
io.savemat('/home/xiangru/Projects/Qixing/TrajSmoothing/motion/data/tmp.mat', dict(clips=clips), do_compression=True)

database = np.load('/home/xiangru/Projects/Qixing/TrajSmoothing/motion/data/tmp.npz')['clips']

print database.shape

if len(sys.argv) >= 4:
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    animation_plot([database[0][start:end+1]], H36=False, repeat=True)
else:
    animation_plot([database[0]], H36=False, repeat=True)

