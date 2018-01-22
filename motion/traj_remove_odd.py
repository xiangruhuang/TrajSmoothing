import sys
import matplotlib.pyplot as plt
import numpy

def remove_odd(trace):

    #diffs = [numpy.linalg.norm(trace[i-1] - trace[i]) for i in range(1, len(trace))]
    #median = numpy.median(diffs)

    #""" remove odd points """
    #i = 0
    #T = median * 10.0
    #final = []
    #last_cutoff = -1
    #keep = [True for pose in trace]
    #while i < len(trace) - 1:
    #    if diffs[i] > T:
    #        if (i - last_cutoff) < 3:
    #            """remove trace[last_cutoff+1:i+1]"""
    #            for t in range(last_cutoff+1, i+1):
    #                keep[t] = False
    #            print('removing frames %d to %d' % (last_cutoff+1, i))
    #        last_cutoff = i
    #    i+=1
    #for i, pose in enumerate(trace):
    #    if keep[i]:
    #        final.append(pose)
    #plt.plot(range(len(diffs)), diffs)
    #plt.plot(range(len(diffs)), [T] * len(diffs))
    #plt.show()
    #return final

    return trace[3:]


if len(sys.argv) < 2:
    print "python traj_remove_odd.py walk/02_01.bvh"
    exit(1)

folder = sys.argv[1].split('/')[0]

with open(sys.argv[1], 'r') as fin, open(sys.argv[1].replace(folder,
    '%s_remove_odd' % folder), 'w') as fout:
    lines = fin.readlines()
    trace = []
    flag = False
    start_point = 0
    for count, line in enumerate(lines):
        if line.startswith('Frame Time'):
            flag = True
            start_point = count
            continue
        if not flag:
            continue
        pose = [float(token) for token in line.strip().split(' ')]
        pose = numpy.asarray(pose)
        trace.append(pose)

    trace = remove_odd(trace)
    trace = numpy.asarray(trace)
    for line in lines[:start_point-1]:
        fout.write(line)
    fout.write('Frames: %d' % len(trace) +'\n')
    fout.write(lines[start_point])
    for pose in trace:
        for t in range(pose.shape[0]):
            if t != 0:
                fout.write(' ')
            fout.write('%f' % pose[t])
        fout.write('\n')


