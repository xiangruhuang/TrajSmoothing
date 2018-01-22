Edinburgh Deep Learning Motion Database
=======================================

Usage
-----

This is a collection of existing databases and internal captures made at 
Edinburgh University retargetted to a skeleton with common structure 
and joint lengths.

As we do not own all the rights to all the individual parts of this data please 
respect the individual usage terms and licenses of the existing databases if 
you wish to make use of them. These existing databases are as follows:

CMU:   http://mocap.cs.cmu.edu/
HDM05: http://resources.mpi-inf.mpg.de/HDM05/
MHAD:  http://tele-immersion.citris-uc.org/berkeley_mhad

All of these databases are free to use, modify and redistribute, but they do 
ask for citation.

The edinburgh databases are also free to use, modify and redistribute (and the 
license is found under `LICENSE.txt`) but we would also ask that you please 
include the following citations in any published work which uses this data:


    @inproceedings{Holden:2015:LMM,
     author = {Holden, Daniel and Saito, Jun and Komura, Taku and Joyce, Thomas},
     title = {Learning Motion Manifolds with Convolutional Autoencoders},
     booktitle = {SIGGRAPH Asia 2015 Technical Briefs},
     year = {2015},
    } 
    
    @inproceedings{Holden:2016:DLF,
     author = {Holden, Daniel and Saito, Jun and Komura, Taku},
     title = {A Deep Learning Framework for Character Motion Synthesis and Editing},
     booktitle = {SIGGRAPH 2016},
     year = {2016},
    }
    
    

Format
------

All of the data is provided in .bvh format. There are many different ways to 
view and process this data format: http://www.cs.man.ac.uk/~toby/bvh/ 

Data has been retargetted to a common skeleton structure primarily using Full 
Body Inverse Kinematics. For this reason some of the rotations of the end 
effector joints may not be preserved and because the process was fully 
automated there may be occasional issues with some data that could not be found 
without manual inspection. Although all the databases use the same joint 
lengths and skeleton structure there may be a little variety in the posing due 
to the target joint positions coming from a different skeleton structure.

Two scripts are provided written in Numpy/Scipy.

The script `export.py` converts the data from BVH format into the format used 
in the papers "A Deep Learning Framework for Character Motion Synthesis and 
Editing" and "Learning Motion Manifolds with Convolutional Autoencoders".

The script `view.py` converts data from this format back into a viewable format 
and displays it using matplotlib.


Databases 
---------

Here is some more information about the individual databases provided:

`cmu` - This is the full cmu database retargetted to a character with uniform 
joint lengths. It contains a huge variety of motions and has been used in lots 
of different research over the years.

`hdm05` - This is the hdm05 database retargetted to the cmu skeleton structure 
and joint lengths. It contains many small clips of individual motions 
originally intended for action classification.

`mhad` - This is the mhad database retargetted to the cmu skeleton structure 
and joint lengths. This database contains just a few actions repeated many 
times by many different subjects.

`edin_locomotion` - This is a database containing large clips of locomotion 
data including running, walking, jogging, and various sidestepping actions. It 
contains around 20 minutes of raw data and is not segmented into individual 
strides.

`edin_kinect` - This is a database containing a large variety of motions 
captured standing in a small area using the kinect motion capture system. 
Because this was captured with the kinect it contains many errors and 
artefacts so should not be used as example data, but may be useful to 
researchers for other purposes.

`edin_xsens` - This is a database containing exactly the same motions as in the 
`edin_kinect` database, but captured using an xsens inertia based motion 
capture system. Because there is a frame-by-frame coorespondence between the 
motion in this database and `edin_kinect` this database may be of interested to 
researchers trying to improve the output of the kinect.

`edin_misc` - This is a small database of various miscellaneous captures made 
at the university including some different walking styles.

`edin_fight` - This is a database of punching, kicking, and fighting motions 
segmented into many small sections.


Contact
-------

For any questions or queries please contact `tkomura@ed.ac.uk`.

