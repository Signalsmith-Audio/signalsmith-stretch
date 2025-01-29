all: out/stretch

out/stretch: ../signalsmith-stretch.h main.cpp util/*.h util/*.hxx ../dsp/*.h
	mkdir -p out
	g++ -std=c++11 -O3 -g \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		main.cpp -o out/stretch

examples: out/stretch
	inputs/run-all.sh out/d1- out/stretch --semitones=-1

clean:
	rm -rf out