all:
	g++ -o preprocess read.cpp -O3 -std=c++11
	g++ -o main main.cpp -O3 -std=c++11

run:
	./main truncated_traj.txt

.PHONY:synthetic
synthetic:
	./main synthetic/trajectories.txt 1e-3

create_gif:
	rm -f recover/*.png
	for i in `ls recover/*.txt | sort -V`; do \
		python plot_traj.py $$i; \
	done
	convert -delay 200 -loop 0 `ls recover/*.png | sort -V` converge.gif
	eopen converge.gif
