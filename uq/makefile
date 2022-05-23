make:
	cd ..; bash install.sh -c
	python gen_inputs.py
	spindec
	rm input_*.txt
	python plot.py
clean:
	rm input_*.txt
