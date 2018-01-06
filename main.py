from util import *

work_file = 'truncated_traj.txt'

xmin = 116.3
xmax = 116.4
ymin = 39.9
ymax = 40.0
n = 7900

work_file2 = 'truncated_traj.temp'

print 'filtering with window [%f, %f] X [%f, %f]' % (xmin, xmax, ymin, ymax)
truncate_traj('trajectories.txt', work_file, xmin, xmax, ymin, ymax)
#truncate_traj(work_file, work_file2, xmin, xmax, ymin, ymax)

#print 'computing k means clusters with k = %d' % n
#kmeans = train_kmeans(n, work_file2)
#print 'predicting and ploting'
#predict_kmeans(kmeans, work_file2)
