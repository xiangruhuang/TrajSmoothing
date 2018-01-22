import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.colors as colors
from matplotlib.animation import ArtistAnimation
import matplotlib.patheffects as pe
from Quaternions import Quaternions
import scipy.io as io
import scipy.ndimage.filters as filters

import os
import BVH
import Animation
from Quaternions import Quaternions

def animation_plot(animations, interval=8.33, H36 = False, repeat=False):
    
    footsteps = []
    for ai in range(len(animations)):
        anim = animations[ai].copy()
        
        joints, root_x, root_z, root_r = anim[:,:-3], anim[:,-3], anim[:,-2], anim[:,-1]
        
        joints = joints.reshape((len(joints), -1, 3))
        
        rotation = Quaternions.id(1)
        offsets = []
        translation = np.array([[0,0,0]])
        
        for i in range(len(joints)):
            joints[i,:,:] = rotation * joints[i]
            joints[i,:,0] = joints[i,:,0] + translation[0,0]
            joints[i,:,2] = joints[i,:,2] + translation[0,2]
            rotation = Quaternions.from_angle_axis(-root_r[i], np.array([0,1,0])) * rotation
            offsets.append(rotation * np.array([0,0,1]))
            translation = translation + rotation * np.array([root_x[i], 0, root_z[i]])
        
        print joints.shape, animations[ai].shape
        animations[ai] = joints
        footsteps.append(anim[:,-4:])
        
    footsteps = np.array(footsteps)
    
    #scale = 1.25*((len(animations))/2)
    scale = 1.25#*((len(animations)))
    #print(animation.writers.list())
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(-scale*30, scale*30)
    ax.set_zlim3d( 0, scale*60)
    ax.set_ylim3d(-scale*30, scale*30)
    ax.set_xticks([], [])
    ax.set_yticks([], [])
    ax.set_zticks([], [])
    ax.set_aspect('equal')
    
    acolors = ['r', 'b']
    #list(sorted(colors.cnames.keys()))[::-1]
    #print acolors
    lines = []
    if ~H36:
        parents = np.array([-1, 0, 1, 2, 0, 4, 5, 0, 7, 8, 9, 9, 11, 12, 9, 14, 15])
    else:
        parents = np.array([-1, 0, 1, 2, 0, 4, 5, 0, 7, 8, 9, 8, 11, 12, 8, 14, 15])
    for ai, anim in enumerate(animations):
        lines.append([plt.plot([0,0], [0,0], [0,0], color=acolors[ai], 
            lw=2, path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])[0] for _ in range(anim.shape[1])])
    
    def animate(i):
        changed = []
        #print i
        #print np.linalg.norm(np.reshape(animations[0], [-1]) - np.reshape(animations[1],
        #    [-1]), 2)
        for ai in range(len(animations)):
            #offset = 25*(ai-((len(animations))/2))
            offset = 0.0 #25*(ai-((len(animations)))/2.0)
            num_frame = animations[ai].shape[0]
            ii = i % num_frame
            for j in range(len(parents)):
                if parents[j] != -1:
                    if i % 2 == 0 and ai == 1:
                        lines[ai][j].set_data([], [])
                        lines[ai][j].set_3d_properties([])
                    else:
                        lines[ai][j].set_data(
                            [ animations[ai][ii,j,0]+offset, animations[ai][ii,parents[j],0]+offset],
                            [-animations[ai][ii,j,2],       -animations[ai][ii,parents[j],2]])
                        lines[ai][j].set_3d_properties(
                            [ animations[ai][ii,j,1],        animations[ai][ii,parents[j],1]])
            changed += lines
            
        #if i == 3:
        #    plt.pause(50)
        return changed
    
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    plt.tight_layout()
    num_frames = np.arange(max([len(ai) for ai in animations]))
    
    ani = animation.FuncAnimation(fig, 
        animate, num_frames,
        interval=interval, repeat=repeat)
    #global count
    #ani.save('%s/videos/%d.mp4' % (sys.argv[1], count), writer=writer)
    #count += 1
    plt.show()

def process_file(filename, window=240, window_step=120, export_trajectory=False, H36 = False ):

    #for CMU and Edinburg
    #anim, names, frametime = BVH.load(filename)

    anim= BVH.load(filename)
    #print(anim.shape)

    """ Subsample to 60 fps """
    #anim = anim[::1]

    """ Do FK """
    global_positions = Animation.positions_global(anim)
    #print global_positions[0,:,:]
    #print(global_positions.shape)

    window=global_positions.shape[0]
    window_step = window

    """ Remove Uneeded Joints """
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
    positions = global_positions[:,used_joints]
    

    """ Put on Floor """
    positions[:,:,1] -= positions[:,:,1].min()
    '''
    """ Get Foot Contacts """
    velfactor, heightfactor = np.array([0.05,0.05]), np.array([3.0, 2.0])
    
    fid_l, fid_r = np.array([3,4]), np.array([7,8])
    feet_l_x = (positions[1:,fid_l,0] - positions[:-1,fid_l,0])**2
    feet_l_y = (positions[1:,fid_l,1] - positions[:-1,fid_l,1])**2
    feet_l_z = (positions[1:,fid_l,2] - positions[:-1,fid_l,2])**2
    feet_l_h = positions[:-1,fid_l,1]
    feet_l = (((feet_l_x + feet_l_y + feet_l_z) < velfactor) & (feet_l_h < heightfactor)).astype(np.float)
    
    feet_r_x = (positions[1:,fid_r,0] - positions[:-1,fid_r,0])**2
    feet_r_y = (positions[1:,fid_r,1] - positions[:-1,fid_r,1])**2
    feet_r_z = (positions[1:,fid_r,2] - positions[:-1,fid_r,2])**2
    feet_r_h = positions[:-1,fid_r,1]
    feet_r = (((feet_r_x + feet_r_y + feet_r_z) < velfactor) & (feet_r_h < heightfactor)).astype(np.float)
    '''
    """ Get Root Velocity """
    velocity = (positions[1:,0:1] - positions[:-1,0:1]).copy()
    
    """ Remove Translation """
    positions[:,:,0] = positions[:,:,0] - positions[:,0:1,0]
    positions[:,:,2] = positions[:,:,2] - positions[:,0:1,2]
    
    """ Get Forward Direction """
    sdr_l, sdr_r, hip_l, hip_r = 11, 14, 1, 4
    across1 = positions[:,hip_l] - positions[:,hip_r]
    across0 = positions[:,sdr_l] - positions[:,sdr_r]
    across = across0 + across1
    across = across / np.sqrt((across**2).sum(axis=-1))[...,np.newaxis]
    
    direction_filterwidth = 20
    forward = np.cross(across, np.array([[0,1,0]]))
    forward = filters.gaussian_filter1d(forward, direction_filterwidth, axis=0, mode='nearest')
    forward = forward / np.sqrt((forward**2).sum(axis=-1))[...,np.newaxis]

    """ Remove Y Rotation """
    target = np.array([[0,0,1]]).repeat(len(forward), axis=0)
    rotation = Quaternions.between(forward, target)[:,np.newaxis]
    positions = rotation * positions
    
    """ Get Root Rotation """
    velocity = rotation[1:] * velocity
    rvelocity = (rotation[1:] * -rotation[:-1]).to_pivots()
    
    """ Add Velocity, RVelocity to vector """
    positions = positions[:-1]
    positions = positions.reshape(len(positions), -1)
    
    positions = np.concatenate([positions, velocity[:,:,0]], axis=-1)
    positions = np.concatenate([positions, velocity[:,:,2]], axis=-1)
    positions = np.concatenate([positions, rvelocity], axis=-1)
    
    """ Add Foot Contacts """
    '''
    positions = np.concatenate([positions, feet_l, feet_r], axis=-1)
    '''
    """ Slide over windows """
    windows = []
    #print len(positions)
    for j in range(0, len(positions), window_step):
    
        slice = positions[j:j+window]
        if len(slice) < window:
            left  = slice[:1].repeat((window-len(slice))//2 + (window-len(slice))%2, axis=0)
            left[:,-7:-4] = 0.0
            right = slice[-1:].repeat((window-len(slice))//2, axis=0)
            right[:,-7:-4] = 0.0
            slice = np.concatenate([left, slice, right], axis=0)
        
        if len(slice) != window: raise Exception()
        windows.append(slice)
        
    return windows

def process_animation(anim, window=240, window_step=120, export_trajectory=False, H36 = False ):

    #for CMU and Edinburg
    #anim, names, frametime = BVH.load(filename)


    """ Subsample to 60 fps """
    #anim = anim[::1]

    """ Do FK """
    global_positions = Animation.positions_global(anim)
    #print global_positions[0,:,:]
    #print(global_positions.shape)

    window=global_positions.shape[0]
    window_step = window

    """ Remove Uneeded Joints """
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
    positions = global_positions[:,used_joints]
    

    """ Put on Floor """
    positions[:,:,1] -= positions[:,:,1].min()
    '''
    """ Get Foot Contacts """
    velfactor, heightfactor = np.array([0.05,0.05]), np.array([3.0, 2.0])
    
    fid_l, fid_r = np.array([3,4]), np.array([7,8])
    feet_l_x = (positions[1:,fid_l,0] - positions[:-1,fid_l,0])**2
    feet_l_y = (positions[1:,fid_l,1] - positions[:-1,fid_l,1])**2
    feet_l_z = (positions[1:,fid_l,2] - positions[:-1,fid_l,2])**2
    feet_l_h = positions[:-1,fid_l,1]
    feet_l = (((feet_l_x + feet_l_y + feet_l_z) < velfactor) & (feet_l_h < heightfactor)).astype(np.float)
    
    feet_r_x = (positions[1:,fid_r,0] - positions[:-1,fid_r,0])**2
    feet_r_y = (positions[1:,fid_r,1] - positions[:-1,fid_r,1])**2
    feet_r_z = (positions[1:,fid_r,2] - positions[:-1,fid_r,2])**2
    feet_r_h = positions[:-1,fid_r,1]
    feet_r = (((feet_r_x + feet_r_y + feet_r_z) < velfactor) & (feet_r_h < heightfactor)).astype(np.float)
    '''
    """ Get Root Velocity """
    velocity = (positions[1:,0:1] - positions[:-1,0:1]).copy()
    
    """ Remove Translation """
    positions[:,:,0] = positions[:,:,0] - positions[:,0:1,0]
    positions[:,:,2] = positions[:,:,2] - positions[:,0:1,2]
    
    """ Get Forward Direction """
    sdr_l, sdr_r, hip_l, hip_r = 11, 14, 1, 4
    across1 = positions[:,hip_l] - positions[:,hip_r]
    across0 = positions[:,sdr_l] - positions[:,sdr_r]
    across = across0 + across1
    across = across / np.sqrt((across**2).sum(axis=-1))[...,np.newaxis]
    
    direction_filterwidth = 20
    forward = np.cross(across, np.array([[0,1,0]]))
    forward = filters.gaussian_filter1d(forward, direction_filterwidth, axis=0, mode='nearest')
    forward = forward / np.sqrt((forward**2).sum(axis=-1))[...,np.newaxis]

    """ Remove Y Rotation """
    target = np.array([[0,0,1]]).repeat(len(forward), axis=0)
    rotation = Quaternions.between(forward, target)[:,np.newaxis]
    positions = rotation * positions
    
    """ Get Root Rotation """
    velocity = rotation[1:] * velocity
    rvelocity = (rotation[1:] * -rotation[:-1]).to_pivots()
    
    """ Add Velocity, RVelocity to vector """
    positions = positions[:-1]
    positions = positions.reshape(len(positions), -1)
    
    positions = np.concatenate([positions, velocity[:,:,0]], axis=-1)
    positions = np.concatenate([positions, velocity[:,:,2]], axis=-1)
    positions = np.concatenate([positions, rvelocity], axis=-1)
    
    """ Add Foot Contacts """
    '''
    positions = np.concatenate([positions, feet_l, feet_r], axis=-1)
    '''
    """ Slide over windows """
    windows = []
    #print len(positions)
    for j in range(0, len(positions), window_step):
    
        slice = positions[j:j+window]
        if len(slice) < window:
            left  = slice[:1].repeat((window-len(slice))//2 + (window-len(slice))%2, axis=0)
            left[:,-7:-4] = 0.0
            right = slice[-1:].repeat((window-len(slice))//2, axis=0)
            right[:,-7:-4] = 0.0
            slice = np.concatenate([left, slice, right], axis=0)
        
        if len(slice) != window: raise Exception()
        windows.append(slice)
        
    return windows

def plot(animations):
    clips = []
    
    for anim in animations:
        clips += process_animation(anim)

    clips += []

    clips = np.array(clips)
    animation_plot([clips[0]])


def aligned_plot(src, tgts, align, interval=8.33, H36 = False, repeat=False):
    
    clips = []
    L = [src] + tgts 
    for anim in L:
        clips += process_animation(anim)

    clips = np.array(clips)
    animations = [clip for clip in clips]
    footsteps = []
    for ai in range(len(animations)):
        anim = animations[ai].copy()
        
        joints, root_x, root_z, root_r = anim[:,:-3], anim[:,-3], anim[:,-2], anim[:,-1]
        
        joints = joints.reshape((len(joints), -1, 3))
        
        rotation = Quaternions.id(1)
        offsets = []
        translation = np.array([[0,0,0]])
        
        for i in range(len(joints)):
            joints[i,:,:] = rotation * joints[i]
            joints[i,:,0] = joints[i,:,0] + translation[0,0]
            joints[i,:,2] = joints[i,:,2] + translation[0,2]
            rotation = Quaternions.from_angle_axis(-root_r[i], np.array([0,1,0])) * rotation
            offsets.append(rotation * np.array([0,0,1]))
            translation = translation + rotation * np.array([root_x[i], 0, root_z[i]])
        
        animations[ai] = joints
        footsteps.append(anim[:,-4:])
        
    footsteps = np.array(footsteps)
    
    #scale = 1.25*((len(animations))/2)
    scale = 1.25#*((len(animations)))
    #print(animation.writers.list())
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(-scale*30, scale*30)
    ax.set_zlim3d( 0, scale*60)
    ax.set_ylim3d(-scale*30, scale*30)
    ax.set_xticks([], [])
    ax.set_yticks([], [])
    ax.set_zticks([], [])
    ax.set_aspect('equal')
    
    acolors = ['r'] + ['b'] * len(tgts)
    #list(sorted(colors.cnames.keys()))[::-1]
    #print acolors
    lines = []
    if ~H36:
        parents = np.array([-1, 0, 1, 2, 0, 4, 5, 0, 7, 8, 9, 9, 11, 12, 9, 14, 15])
    else:
        parents = np.array([-1, 0, 1, 2, 0, 4, 5, 0, 7, 8, 9, 8, 11, 12, 8, 14, 15])
    for ai, anim in enumerate(animations):
        print anim.shape
        lines.append([plt.plot([0,0], [0,0], [0,0], color=acolors[ai], 
            lw=2, path_effects=[pe.Stroke(linewidth=3, foreground='black'), pe.Normal()])[0] for _ in range(anim.shape[1])])

    
    def animate(i):
        changed = []
        num_frame = animations[0].shape[0]
        ii = i % num_frame
        if align[ii] is not None:
            tnum, fnum = align[ii]
        for ai in range(len(animations)):
            #offset = 25*(ai-((len(animations))/2))
            offset = 0.0 #25*(ai-((len(animations)))/2.0)
            for j in range(len(parents)):
                if parents[j] != -1:
                    if ai == 0:
                        lines[ai][j].set_data(
                            [ animations[ai][ii,j,0]+offset, animations[ai][ii,parents[j],0]+offset],
                            [-animations[ai][ii,j,2],       -animations[ai][ii,parents[j],2]])
                        lines[ai][j].set_3d_properties(
                            [ animations[ai][ii,j,1],        animations[ai][ii,parents[j],1]])
                    elif (ai == tnum+1) and (align[ii] is not None):
                        lines[ai][j].set_data( 
                            [animations[ai][fnum,j,0]+offset, animations[ai][fnum,parents[j],0]+offset],
                            [-animations[ai][fnum,j,2],     -animations[ai][fnum,parents[j],2]])
                        lines[ai][j].set_3d_properties(
                            [ animations[ai][fnum,j,1],     animations[ai][fnum,parents[j],1]])
                    else:
                        lines[ai][j].set_data([], []) 
                        lines[ai][j].set_3d_properties([])

            changed += lines
            
        #if i == 3:
        #    plt.pause(50)
        return changed
    
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    plt.tight_layout()
    num_frames = np.arange(max([len(ai) for ai in animations]))
    
    ani = animation.FuncAnimation(fig, 
        animate, num_frames,
        interval=interval, repeat=repeat)
    #global count
    #ani.save('%s/videos/%d.mp4' % (sys.argv[1], count), writer=writer)
    #count += 1
    plt.show()
