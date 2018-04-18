import sys

with open(sys.argv[1], 'r') as list_in, open(sys.argv[2], 'w') as fout:
    for filename in list_in.readlines():
        with open(filename.strip(), 'r') as fin:
            line = fin.readlines()[0].strip()
            fout.write(line+'\n')

