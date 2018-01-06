import sys
import matplotlib.pyplot as plt

#xmin = 116.0
#xmax = 116.8
#ymin= 39.6
#ymax = 40.2
#xmin = 116.375
#xmax = 116.4
#ymin= 39.875
#ymax = 39.9

def in_range(x, y, xmin, xmax, ymin, ymax):
    if (x <= xmax) and (x >= xmin):
        if (y <= ymax) and (y >= ymin):
            return True
    return False

def truncate_traj(input_file, output_file, xmin, xmax, ymin, ymax):
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        lines = fin.readlines()

        for count, line in enumerate(lines):
            token = line.strip().split(' ')
            
            x = token[2::3]
            y = token[1::3]
            x = [float(x_i) for x_i in x]
            y = [float(y_i) for y_i in y]
            
            workx = []
            worky = []
            for x_i, y_i in zip(x, y):
                if in_range(x_i, y_i, xmin, xmax, ymin, ymax):
                    workx.append(x_i)
                    worky.append(y_i)
                else:
                    if len(workx) > 20:
                        plt.plot(workx, worky, '+')
                        _count = 0
                        for x_i, y_i in zip(workx, worky):
                            if _count != 0:
                                fout.write(' ')
                            fout.write('%f %f' % (x_i, y_i))
                            _count += 1
                        fout.write('\n')
                    workx = []
                    worky = []
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
        plt.show()
