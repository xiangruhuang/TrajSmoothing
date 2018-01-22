all:
	g++ -o preprocess read.cpp -O3 -std=c++11
	g++ -o main main.cpp -O3 -std=c++11

run:
	./main truncated_traj.txt

.PHONY:synthetic
synthetic:
	./main synthetic/trajectories.txt 2 1e-3

motion_dir=/home/xiangru/Projects/Qixing/AnimationSequences/animationVideo/preprocess/walk_remove_odd/trajs.txt

lambda=1e-3
threshold=20.0
human_motion:
	./main $(motion_dir) 93 $(lambda) $(threshold) human

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
