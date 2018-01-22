import os
import sys
import numpy as np
import scipy.io as io
import scipy.ndimage.filters as filters

sys.path.append('D:/Projects/motion/motion')

import BVH
import Animation
from Quaternions import Quaternions


def process_file(filename, window=240, window_step=120, export_trajectory=False, H36 = False):
    # for CMU and Edinburg
    # anim, names, frametime = BVH.load(filename)
    anim = BVH.load(filename, bvh = True)


    """ Subsample to 60 fps """
    #anim = anim[::2]
    print(anim.shape)

    if ~H36:
        used_joints = np.array([
            0,
            2, 3, 4,
            7, 8, 9,
            12, 13, 15, 16,
            18, 19, 20,
            25, 26, 27])
    else:
        used_joints = np.array([
            0,
            2, 3, 4,
            7, 8, 9,
            12, 13, 15, 17,
            18, 19, 20,
            25, 26, 27])

    positions = anim[:, used_joints]
    window=positions.shape[0]
    window_step = window
    print window


    """ Slide over windows """
    windows = []

    for j in range(0, len(positions), window_step):

        slice = positions[j:j + window]
        if len(slice) < window:
            left = slice[:1].repeat((window - len(slice)) // 2 + (window - len(slice)) % 2, axis=0)
            left[:, -7:-4] = 0.0
            right = slice[-1:].repeat((window - len(slice)) // 2, axis=0)
            right[:, -7:-4] = 0.0
            slice = np.concatenate([left, slice, right], axis=0)

        if len(slice) != window: raise Exception()
        windows.append(slice)

    return windows



#for database in [
#    'cmu', 'hdm05', 'mhad',
#    'edin_locomotion', 'edin_kinect',
#    'edin_xsens', 'edin_misc',
#    'edin_fight']:

for database in [sys.argv[1]]:
    path = ''
    files = [os.path.join(path + database, f) for f in sorted(list(os.listdir(path + database)))
             if os.path.isfile(os.path.join(path + database, f))
             and f.endswith('.bvh') and f != 'rest.bvh']

    clips = []

    for i, item in enumerate(files):
        print('Database %s Processing %i of %i (%s)' % (database, i, len(files), item))
        clips += process_file(item, H36 = False)

    clips = np.array(clips)
    np.savez_compressed('data/' + database, clips=clips)
    io.savemat('data/' + database + '.mat', dict(clips=clips), do_compression=True)

