import sys
import os
import numpy

def bvh2traj(filelist, output_file):
    with open(output_file, 'w') as fout:
        for filename in filelist:
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
                #traj = traj[3:] # remove first three frames
                traj = numpy.asarray(traj)
                traj = traj[:, 3:]
                print filename, traj.shape
                traj = numpy.reshape(traj, [-1])
            for i, num in enumerate(traj):
                if i != 0:
                    fout.write(' ')
                fout.write('%f' % num)
            fout.write('\n')


print('converting folder %s' % sys.argv[1])
filelist = ['/'.join([sys.argv[1], f]) for f in sorted(os.listdir(sys.argv[1])) if f.split('.')[-1] == 'bvh']
with open(sys.argv[1]+'/list.txt', 'w') as fout:
    for filename in filelist:
        fout.write(filename +'\n')
bvh2traj(filelist, sys.argv[1]+'/trajs.txt')
