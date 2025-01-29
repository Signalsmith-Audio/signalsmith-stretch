all: out/stretch

out/stretch: ../signalsmith-stretch.h main.cpp util/*.h util/*.hxx ../dsp/*.h
	mkdir -p out
	g++ -std=c++11 -O3 -g \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		main.cpp -o out/stretch

examples: out/stretch
	mkdir -p out/examples
	inputs/run-all.sh out/examples/u2- out/stretch --semitones=2
	inputs/run-all.sh out/examples/d2- out/stretch --semitones=-2
	inputs/run-all.sh out/examples/u4- out/stretch --semitones=4
	inputs/run-all.sh out/examples/d4- out/stretch --semitones=-4
	inputs/run-all.sh out/examples/u8- out/stretch --semitones=8
	inputs/run-all.sh out/examples/d8- out/stretch --semitones=-8
	inputs/run-all.sh out/examples/u16- out/stretch --semitones=16
	inputs/run-all.sh out/examples/d16- out/stretch --semitones=-16
	inputs/run-all.sh out/examples/t_8- out/stretch --time=0.8
	inputs/run-all.sh out/examples/t1_2- out/stretch --time=1.2
	inputs/run-all.sh out/examples/t1_5- out/stretch --time=1.5
	inputs/run-all.sh out/examples/t2- out/stretch --time=2
	inputs/run-all.sh out/examples/t4- out/stretch --time=4

clean:
	rm -rf out

### Example use of CMake

cmake:
	cmake -B build -G Xcode -D CMAKE_RUNTIME_OUTPUT_DIRECTORY=../out
	cmake --build build --config Release