all: out/stretch

DEV_FLAGS := --semitones=4 --time=1

dev: out/stretch
	./out/stretch inputs/dev.wav out/dev-2048.wav --process-chunk=2048 $(DEV_FLAGS)
	./out/stretch inputs/dev.wav out/dev-512.wav --process-chunk=512 $(DEV_FLAGS)
	./out/stretch inputs/dev.wav out/dev-100.wav --process-chunk=100 $(DEV_FLAGS)
	./out/stretch inputs/dev.wav out/dev-2048-sc.wav --process-chunk=2048 --split-computation $(DEV_FLAGS)
	./out/stretch inputs/dev.wav out/dev-512-sc.wav --process-chunk=512 --split-computation $(DEV_FLAGS)
	./out/stretch inputs/dev.wav out/dev-100-sc.wav --process-chunk=100 --split-computation $(DEV_FLAGS)

out/stretch: main.cpp ../signalsmith-stretch.h util/*.h util/*.hxx
	mkdir -p out
	g++ -std=c++11 -O3 -g \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		-framework Accelerate -DSIGNALSMITH_USE_ACCELERATE \
		-I ../include \
		main.cpp -o out/stretch

# Uses input files from: https://signalsmith-audio.co.uk/code/stretch/inputs.zip
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

TEST_WAV ?= "inputs/voice.wav"

clean:
	rm -rf out

### Example use of CMake

cmake:
	# CMAKE_BUILD_TYPE is needed for single-config generators (e.g. Makefiles)
	cmake -B out/build -DCMAKE_BUILD_TYPE=Release
	cmake --build out/build --config Release
