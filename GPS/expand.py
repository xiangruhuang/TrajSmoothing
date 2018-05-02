import sys
import os

with open(sys.argv[1], 'r') as fin, open('list.txt', 'w') as list_output:
    lines = fin.readlines()
    for count, line in enumerate(lines):
        filename = '%d.txt' % count
        with open(filename, 'w') as fout:
            fout.write(line.strip()+'\n')
        list_output.write(os.getcwd()+'/'+filename+'\n')

