import sys
import numpy as np

from util import animation_plot
from util import process_file
import scipy.io as io

if len(sys.argv) < 3:
    print "python <bvh file1> <bvh file2>"
    exit(1)

clips = []

clips += process_file(sys.argv[1])
clips += process_file(sys.argv[2])

clips = np.array(clips)
print clips[0].shape, type(clips[0])
print clips.shape, type(clips)
np.savez_compressed('data/compare', clips=clips)
io.savemat('data/compare.mat', dict(clips=clips), do_compression=True)

database = np.load('data/compare.npz')['clips']


#print type(database[0])
#print type(clips)
#print clips.shape

database = clips
print database[0].shape, type(database[0])
print database.shape, type(database)
print clips[0].shape, type(clips[0])
print clips.shape, type(clips)

animation_plot([database[0], database[1]], H36=False, repeat=True)


