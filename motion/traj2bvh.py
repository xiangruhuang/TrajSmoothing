import sys
import os
import numpy

def traj2bvh(workdir, traj_input, bvh_output):
    with open(workdir+'/list.txt', 'r') as fin:
        lines = fin.readlines()
        filelist = [line.strip() for line in lines]
    
    with open(traj_input, 'r') as fin:
        lines = fin.readlines()
        traces = []
        for line in lines:
            s = numpy.asarray([float(token) for token in line.strip().split(' ')])
            s = numpy.reshape(s, [-1, 93])
            traces.append(s)
    
    for filename, trace in zip(filelist, traces):
        with open(filename, 'r') as fin:
            lines = fin.readlines()
            traj = []
            flag = False
            for line in lines:
                if line.strip().startswith('Frame Time'):
                    flag = True
                    continue
                if not flag:
                    continue
                pose = [float(token) for token in line.strip().split(' ')]
                traj.append(pose)
            traj = numpy.asarray(traj)
            traj = traj[:, :3]
            traj = numpy.concatenate((traj, trace), axis=1)
        notdir = filename.split('/')[-1]
        with open(bvh_output+'/'+notdir, 'w') as fout:
            for line in lines:
                fout.write(line)
                if line.strip().startswith('Frame Time'):
                    break
            for i in range(traj.shape[0]):
                for j in range(traj.shape[1]):
                    if j != 0:
                        fout.write(' ')
                    fout.write('%f' % traj[i, j])
                fout.write('\n')
        

if len(sys.argv) < 3:
    print 'python traj2bvh.py <trajectory input folder> <bvh output folder>'
    exit(1)

traj2bvh('walk_remove_odd', sys.argv[1], sys.argv[2])

#print('converting folder %s' % sys.argv[1])
#filelist = ['/'.join([sys.argv[1], f]) for f in os.listdir(sys.argv[1]) if f.split('.')[-1] == 'bvh']
#with open(sys.argv[1]+'/list.txt', 'w') as fout:
#    for filename in filelist:
#        fout.write(filename +'\n')
#bvh2traj(filelist, sys.argv[1]+'/trajs.txt')
