
name=walk

all:
	make data/$(name).npz
	mkdir -p $(name)/videos
	#python view.py $(name)
	#scp -r $(name)/videos xrhuang@narsil-4.cs.utexas.edu:./public_html/animations/$(name)

%.mp4:
	python view.py $(name)

%.npz:
	python export.py $(notdir $(basename $@))

convert:
	python bvh2traj.py $(name)

remove_odd:
	for i in `ls $(name)/*.bvh`; do \
		echo $$i; \
		python traj_remove_odd.py $$i; \
	done

%/lists.txt:
	$(eval folder:= $(dir $@))
	python bvh2traj.py $(folder)

num=100
traj2bvh:
	mkdir -p walk_smoothed_$(num)
	python traj2bvh.py /home/xiangru/Projects/Qixing/TrajSmoothing/human/$(num).txt walk_smoothed_$(num)

