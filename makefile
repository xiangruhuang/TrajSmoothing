all:
	#g++ -o preprocess read.cpp -O3 -std=c++11
	#g++ -o main main.cpp -O3 -std=c++11 -Wunused-result

.PHONY: align
align:
	g++ -o align align.cpp -O3 -std=c++11 -fopenmp -Wunused-result

.PHONY: motion_smooth
motion_smooth:
	g++ -o motion_smooth motion_smooth.cpp -O3 -std=c++11 -fopenmp -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

.PHONY: knn
knn:
	g++ -o knn knn.cpp -O3 -std=c++11 -fopenmp -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN



run:
	./main truncated_traj.txt

.PHONY:synthetic
synthetic:
	./main synthetic/trajectories.txt 2 1e-3

motion_dir=/home/xiangru/Projects/Qixing/AnimationSequences/animationVideo/preprocess/walk_remove_odd/trajs.txt

lambda=1e-4
threshold=20.0
#list=motion/motion_list.txt
list=motion/lists/small.txt

.PHONY:smooth
smooth:
	./motion_smooth -c equal_dist -l 1e-3 -o 10 -k 10 $(list) motion/smooth 

create_gif:
	rm -f recover/*.png
	for i in `ls recover/*.txt | sort -V`; do \
		make `echo $$i | sed 's/txt/eps/'`; \
	done
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=./converge.pdf `ls recover/*.eps | sort -V`
	#convert -delay 200 -loop 0 `ls recover/*.png | sort -V` converge.gif
	evince converge.pdf

%.eps:
	python plot_traj.py $(basename $@).txt

.PHONY: align
align_walk:
	./align motion/motion_list.txt 14 3000.0 walk

align_fight:
	./align motion/motion_list 2545 3000.0 fight

align_run:
	./align motion 97 3000.0 run
