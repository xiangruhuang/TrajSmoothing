import os
import sys
#sys.path.append('D:/Projects/motion/motion')
import numpy as np
import scipy.io as io

from util import process_file

#for database in [
#    'cmu', 'hdm05', 'mhad', 
#    'edin_locomotion', 'edin_kinect', 
#    'edin_xsens', 'edin_misc', 
#    'edin_fight']:

for database in [sys.argv[1]]:
    path = ''
    files = [os.path.join(path+database,f) for f in sorted(list(os.listdir(path+database)))
        if os.path.isfile(os.path.join(path+database,f))
        and f.endswith('.bvh') and f != 'rest.bvh']
    
    clips = []
    
    for i, item in enumerate(files):
        print('Database %s Processing %i of %i (%s)' % (database, i, len(files), item))
        item2 = item.replace(sys.argv[1], sys.argv[1]+'_smoothed')
        clips += process_file(item, H36 = False)
        #clips += process_file(item2, H36 = False)
    
    clips = np.array(clips)
    np.savez_compressed('data/'+database, clips=clips)
    io.savemat('data/'+database+'.mat', dict(clips=clips), do_compression=True)
