
name=walk_remove_odd

all:
	make data/$(name).npz
	mkdir -p $(name)/videos
	make $(name)/trajs.txt

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

%/trajs.txt:
	$(eval folder:= $(dir $@))
	python bvh2traj.py $(folder)

num=100
traj2bvh:
	mkdir -p walk_smoothed_$(num)
	python traj2bvh.py /home/xiangru/Projects/Qixing/TrajSmoothing/human/$(num).txt walk_smoothed_$(num)

align-walk:
	python align.py motion_list.txt walk

align-fight:
	python align.py motion_list.txt fight

align-run:
	python align.py motion_list.txt videos/run
