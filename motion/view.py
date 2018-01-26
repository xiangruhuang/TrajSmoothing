import sys

from util import animation_plot
import numpy as np

count = 0
        
# select the plot

database = np.load('data/%s.npz' % sys.argv[1])['clips']

print type(database), database[3].shape

def compare(idx1, idx2, f1, f2):
    animation_plot([database[idx1][f1:f1+1, :], database[idx2][f2:f2+1, :]],
            repeat=True)

def view_matchings(a, b):
    with open('/home/xiangru/Projects/Qixing/TrajSmoothing/human/matchings', 'r') as fin:
        lines = fin.readlines()
        pairs = lines[0::2]
        matchings = lines[1::2]
        for pair, match in zip(pairs, matchings):
            i, j, T = [int(token) for token in pair.strip().split(' ')]
            if T == 0:
                continue
            tuples = [int(token) for count, token in
                    enumerate(match.strip().split(' ')) if (count % 3 != 2)]
            assert T*2 == len(tuples)
            
            if not (i == a and j == b):
                continue
            for a, b in zip(tuples[0::2], tuples[1::2]):
                compare(i, j, a+3, b+3)
               
#a = int(sys.argv[2])
#b = int(sys.argv[3])
#print "comparing animations a b"
#
#view_matchings(a, b)

f = 


#print type(database[0]), database[0].shape 
#
#step_size = 2
#for index in range(0, len(database), step_size):
#    print 'video', index
#    L = [database[index+i] for i in range(step_size)]
#    animation_plot(L, H36 = False, repeat = True)

#L = [database[a], database[b]]
#animation_plot(L, H36 = False, repeat = True)
