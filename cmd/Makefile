all: out/stretch

out/stretch: ../signalsmith-stretch.h main.cpp util/*.h util/*.hxx ../dsp/*.h
	mkdir -p out
	g++ main.cpp -o out/stretch -std=c++11 -O3 -g

examples: out/stretch
	inputs/run-all.sh out/d1- out/stretch --semitones=-1

clean:
	rm -rf out