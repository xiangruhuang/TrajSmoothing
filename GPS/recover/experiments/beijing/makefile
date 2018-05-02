plots=raw_points dds_points dds_graph ours_points ours_graph

all:
	for i in $(plots); do \
		make $$i.epsplot; \
		make $$i.pngplot; \
	done

raw_points:
	cp ../../../real/beijing.txt $@

dds_points:
	cp ../../real/dds_real_delay1_K5/dds_iter10.txt $@

recover:
	cp /home/xiangru/Projects/Qixing/TrajSmoothing/recover ./

temp/expand.py:
	mkdir -p temp
	cp /home/xiangru/Projects/Qixing/TrajSmoothing/GPS/expand.py temp/expand.py

dds_graph: dds_points recover temp/expand.py
	mkdir -p temp
	cd temp/ && python expand.py ../$<
	./recover -k 100 -t 100.0 -radius 0.00001 -o 100 -lambda 0.0 -max_iter 0 temp/list.txt $@

ours_points:
	cp ../../real/real_o6_K10_calign1.0_csmooth1e-5_creg0.3_interval0.0005/GPS_20.txt $@

ours_graph: ours_points recover temp/expand.py
	mkdir -p temp
	cd temp/ && python expand.py ../$<
	./recover -k 120 -t 30 -radius 0.00001 -o 100000 -lambda 10.0 -max_iter 15 temp/list.txt $@

%.epsplot: %
	mkdir -p figures;
	python plot.py $< figures/$<.eps

%.pngplot: %
	mkdir -p figures;
	python plot.py $< figures/$<.png
