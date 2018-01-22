all:
	g++ main.cpp -o dds -O3 -std=c++11 -I /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/include -L /home/xiangru/Projects/Qixing/TrajSmoothing/DDS/ann/ann_1.1.2/lib -lANN

K=50
n=1

grids:
	$(eval folder := $@_knn_delay$(n)_K$(K))
	mkdir -p $(folder)
	./dds ../synthetic/$@.txt $(n) $(K) $(folder)/

circle:
	$(eval folder := $@_knn_delay$(n)_K$(K))
	mkdir -p $(folder)
	./dds ../synthetic/$@.txt $(n) $(K) $(folder)/

plot_grids:
	$(eval folder := grids_knn_delay$(n)_K$(K))
	for iter in `seq 0 9`; do \
		make $(folder)/_iter$$iter.eps; \
	done
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=converge.pdf `ls $(folder)/*.eps | sort -V`
	evince converge.pdf

plot_circle:
	$(eval folder := circle_knn_delay$(n)_K$(K))
	for iter in `seq 0 9`; do \
		make $(folder)/_iter$$iter.eps; \
	done
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=converge.pdf `ls $(folder)/*.eps | sort -V`
	#xdg-open knn_delay$(n)_K$(K)/_iter0.png
	evince converge.pdf

%.eps:
	echo $(basename $@)
	python plot_traj.py $(basename $@).out
