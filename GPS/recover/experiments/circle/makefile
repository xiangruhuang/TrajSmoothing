plots=raw_points raw_graph dds_points dds_graph ours_points ours_graph

all:
	for i in $(plots); do \
		make $$i.epsplot; \
		make $$i.pngplot; \
	done

raw_points:
	cp ../../../synthetic/circle.txt $@

raw_graph:
	cp ../../../synthetic/circle.gt $@

dds_points:
	cp ../../synthetic/circle/synthetic_delay1_K50/dds_circle_iter5.txt $@

recover:
	cp /home/xiangru/Projects/Qixing/TrajSmoothing/recover ./

temp/expand.py:
	mkdir -p temp
	cp /home/xiangru/Projects/Qixing/TrajSmoothing/GPS/expand.py temp/expand.py

dds_graph: dds_points recover temp/expand.py
	mkdir -p temp
	cd temp/ && python expand.py ../$<
	./recover -k 100 -t 50.0 -radius 0.001 -lambda 0.0 -o 100 temp/list.txt $@

ours_points:
	cp ../../synthetic/circle/synthetic_o3_radius0.01_calign1.0_csmooth1e-10_creg1e-3_interval0.1/GPS_20.txt $@

ours_graph: ours_points recover temp/expand.py
	mkdir -p temp
	cd temp/ && python expand.py ../$<
	./recover -k 10 -t 30 -radius 0.001 -lambda 0.0 -o 100 temp/list.txt $@

%.epsplot: %
	mkdir -p figures;
	python plot.py $< figures/$<.eps

%.pngplot: %
	mkdir -p figures;
	python plot.py $< figures/$<.png
