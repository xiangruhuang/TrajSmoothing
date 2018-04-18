from util import *

work_file = 'truncated_traj.txt'

xmin = 116.0
xmax = 116.4
ymin = 39.0
ymax = 40.2
n = 7900

print 'filtering with window [%f, %f] X [%f, %f]' % (xmin, xmax, ymin, ymax)
truncate_traj('real/trajectories.txt', work_file, xmin, xmax, ymin, ymax)
#truncate_traj(work_file, work_file2, xmin, xmax, ymin, ymax)

#print 'computing k means clusters with k = %d' % n
#kmeans = train_kmeans(n, work_file2)
#print 'predicting and ploting'
#predict_kmeans(kmeans, work_file2)
