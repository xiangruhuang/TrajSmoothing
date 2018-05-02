grids:
	python create_data.py grids
	python expand.py $@.txt

circle:
	python create_data.py circle
	python expand.py $@.txt
	
street:
	python create_data.py $@
	python expand.py $@.txt

line:
	python create_data.py $@
	python expand.py $@.txt

read:
	python create_data.py ../../picture_generation/traces.txt 
	python expand.py $@.txt
