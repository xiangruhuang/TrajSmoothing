all:
	python main.py

folder=o1_K10

plot:
	#python merge.py recover/$(folder)/GPS_300.txt recover
	python plot_traj.py recover/$(folder)/GPS_300.txt figures/$(folder).png

.PHONY: filter
filter:
	g++ filter_with_window.cpp -o filter -O3 -std=c++11 

filter_with_window:
	#./filter real/trajectories.txt trajs.txt 116.0 116.4 39.0 40.2 2
	./filter real/trajectories.txt trajs4_0.txt 116.28 116.30 39.99 40.01 2
	#./filter trajs.txt trajs2.txt 116.21 116.24 39.97 40.0 2
	#./filter trajs.txt trajs4.txt 116.28 116.30 39.99 40.01 2
	#./filter trajs4.txt trajs4_0.txt 116.280 116.30 39.99 40.01 2
