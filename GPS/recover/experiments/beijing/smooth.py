import sys
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import numpy as np

with open(sys.argv[1], 'r') as fin:
    lines = fin.readlines()

    mode = 0
    P = []
    
    lamb = 1e-9

    row = []
    col = []
    data = []
    
    edges = []

    for line in lines:
        if 'Point' in line:
            mode = 1
            continue
        if 'Edge' in line:
            mode = 2
            N = len(P)
            d = [0.0] * N
            continue
        assert mode > 0
        if mode == 1:
            p = [float(token) for token in line.strip().split(' ')]
            P.append(p)
        if mode == 2:
            a, b, w = [float(token) for token in line.strip().split(' ')]
            a = int(a)
            b = int(b)
            edges.append([a, b, w])
            assert a != b
            row.append(a)
            col.append(b)
            data.append(-1.0*lamb*w/100.0)
            d[a] += w/100.0
            d[b] += w/100.0

    for i in range(N):
        row.append(i)
        col.append(i)
        data.append(d[i]*lamb + 1.0)
    row = np.array(row)
    col = np.array(col)
    data = np.array(data)
    P = np.array(P)

    A = csr_matrix((data, (row, col)), shape=(N, N))
    P_smooth = spsolve(A, P)

with open(sys.argv[2], 'w') as fout:
    fout.write('Points:\n')
    for i in range(N):
        for j in range(2):
            if j != 0:
                fout.write(' ')
            fout.write('%.10f' % P_smooth[i, j])
        fout.write('\n')
    fout.write('Edges:\n')
    for i in range(len(edges)):
        a, b, w = edges[i]
        fout.write('%d %d %f\n' % (a, b, w))

