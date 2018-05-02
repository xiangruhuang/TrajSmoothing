all:
	g++ -o preprocess read.cpp -O3 -std=c++11
	#g++ -o main main.cpp -O3 -std=c++11 -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

.PHONY: align
align:
	g++ -o align align.cpp -O3 -std=c++11 -fopenmp -Wunused-result

.PHONY: motion_smooth
motion_smooth:
	g++ -g -o motion_smooth motion_smooth.cpp -O3 -std=c++11 -fopenmp -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

.PHONY: recover
recover:
	g++ -g -o recover recover.cpp -O3 -std=c++11 -fopenmp -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN
	cp recover GPS/recover/

.PHONY: knn
knn:
	g++ -o knn knn.cpp -O3 -std=c++11 -fopenmp -Wunused-result -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

lambda=1e-2
K=10
o=3
solver=AM
c_reg=1e-3
c_smooth=1e-10
c_align=1.0
interval=0.1
rad=0.01

.PHONY:test
test:
	mkdir -p GPS/recover/$@_o$(o)_K$(K)
	./motion_smooth -solver $(solver) -lr $(lambda) -o 0 -k 1 -D 2 -tol 1e-8 -method equal_dist GPS/$@/list.txt recover/$@_o$(o)_K$(K)

delay=1
dds_K=50

synthetic/circle:
	mkdir -p GPS/recover/$@/synthetic_o$(o)_radius$(rad)_calign$(c_align)_csmooth$(c_smooth)_creg$(c_reg)_interval$(interval)
	./motion_smooth -radius $(rad) -sigma 0.05 -interval 0.1 -c_align $(c_align) -c_smooth $(c_smooth) -c_reg $(c_reg) -tol 1e-5 -solver $(solver) -lr $(lambda) -o $(o) -k $(K) -D 2 -method shortest_path GPS/synthetic/list.txt recover/$@/synthetic_o$(o)_radius$(rad)_calign$(c_align)_csmooth$(c_smooth)_creg$(c_reg)_interval$(interval)

	
dds_synthetic/circle:
	mkdir -p GPS/recover/synthetic/circle/synthetic_delay$(delay)_K$(dds_K)
	./dds GPS/synthetic/circle.txt $(delay) $(dds_K) GPS/recover/synthetic/circle/synthetic_delay$(delay)_K$(dds_K)/dds_circle

dds_real:
	mkdir -p GPS/recover/real/$@_delay$(delay)_K$(dds_K)
	./dds GPS/real/real.txt $(delay) $(dds_K) GPS/recover/real/$@_delay$(delay)_K$(dds_K)/dds

.PHONY:real
real:
	mkdir -p GPS/recover/$@/$@_o$(o)_K$(K)_calign$(c_align)_csmooth$(c_smooth)_creg$(c_reg)_interval$(interval)
	./motion_smooth -sigma 0.0003 -interval $(interval) -solver $(solver) -lr $(lambda) -o $(o) -k $(K) -c_align $(c_align) -c_smooth $(c_smooth) -c_reg $(c_reg) -D 2 -tol 1e-3 -method shortest_path GPS/$@/list.txt recover/$@/$@_o$(o)_K$(K)_calign$(c_align)_csmooth$(c_smooth)_creg$(c_reg)_interval$(interval)

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
