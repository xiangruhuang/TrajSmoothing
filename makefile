all:
	#g++ -o preprocess read.cpp -O3 -std=c++11
	#g++ -o main main.cpp -O3 -std=c++11 -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

.PHONY: align
align:
	g++ -o align align.cpp -O3 -std=c++11 -fopenmp -Wunused-result

.PHONY: motion_smooth
motion_smooth:
	g++ -g -o motion_smooth motion_smooth.cpp -O3 -std=c++11 -fopenmp -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

.PHONY: knn
knn:
	g++ -o knn knn.cpp -O3 -std=c++11 -fopenmp -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

lambda=1e-2
K=5
o=1
solver=AM

.PHONY:test
test:
	mkdir -p GPS/recover/$@_o$(o)_K$(K)
	./motion_smooth -solver $(solver) -lr $(lambda) -o 0 -k 1 -D 2 -tol 1e-8 -method equal_dist GPS/$@/list.txt recover/$@_o$(o)_K$(K)

delay=1
dds_K=3

synthetic:
	mkdir -p GPS/recover/$@_o$(o)_K$(K)
	./motion_smooth -interval 0.1 -c_align 10.0 -c_smooth 10 -c_reg 1.0 -tol 1e-3 -solver $(solver) -lr $(lambda) -o $(o) -k $(K) -D 2 -method interpolate GPS/$@/list.txt recover/$@_o$(o)_K$(K)
	
dds_synthetic:
	mkdir -p GPS/recover/$@_delay$(delay)_K$(dds_K)
	./dds GPS/synthetic/street.txt $(delay) $(dds_K) GPS/recover/$@_delay$(delay)_K$(dds_K)/dds

dds_real:
	mkdir -p GPS/recover/$@_delay$(delay)_K$(dds_K)
	./dds GPS/real/fake $(delay) $(dds_K) GPS/recover/$@_delay$(delay)_K$(dds_K)/dds

.PHONY:real
real:
	mkdir -p GPS/recover/$@_o$(o)_K$(K)
	./motion_smooth -interval 0.001 -solver $(solver) -lr $(lambda) -o $(o) -k $(K) -c_align 10.0 -c_smooth 10.0 -c_reg 1.0 -D 2 -tol 1e-8 -method interpolate GPS/$@/list.txt recover/$@_o$(o)_K$(K)

threshold=20.0

motion_dir=/home/xiangru/Projects/Qixing/AnimationSequences/animationVideo/preprocess/walk_remove_odd/trajs.txt

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
