all:
	g++ main.cpp -o dds -O3 -std=c++11 -I /home/xiangru/Projects/Qixing/DictLearning/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/DictLearning/DDS/ann/ann_1.1.2/lib -lANN

K=50
n=5
batch:
	for i in `seq 1 9`; do \
		./dds ../synthetic/trajectories.txt $$i $K knn_delay$$i_K$K; \
	done

run:
	./dds ../synthetic/trajectories.txt $(n) $(K) knn_delay$(n)_K$(K)

plot:
	for iter in `seq 0 9`; do \
		make knn_delay$(n)_K$(K)_iter$$iter.png; \
	done
	xdg-open knn_delay$(n)_K$(K)_iter0.png

%.png:
	echo $(basename $@)
	python plot_traj.py $(basename $@).out
