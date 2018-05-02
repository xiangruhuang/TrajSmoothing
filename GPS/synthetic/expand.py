import sys
sys.path.append('../')
from util import repeat

with open(sys.argv[1], 'r') as fin, open('list.txt', 'w') as list_output:
    lines = fin.readlines()
    assert len(lines) % repeat == 0
    num_trace = len(lines) / repeat
    for count, line in enumerate(lines):
        trace_id = count / repeat
        sample_id = count % repeat
        filename = '%d_%d.txt' % (trace_id, sample_id)
        with open(filename, 'w') as fout:
            fout.write(line.strip()+'\n')
        list_output.write('/home/xiangru/Projects/Qixing/TrajSmoothing/GPS/synthetic/'+filename+'\n')

