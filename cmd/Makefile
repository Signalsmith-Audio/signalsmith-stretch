all: out/stretch

out/stretch: Makefile ../signalsmith-stretch.h main.cpp util/*.h
	mkdir -p out
	g++ main.cpp -o out/stretch -std=c++11 -O3 -g

clean:
	rm -rf out

TEST_WAV ?= "inputs/dno-solo.wav"

dev: out/stretch
	out/stretch --time=0.8 --semitones=10 $(TEST_WAV) out/shift.wav
	out/stretch --time=0.8 --semitones=10 --formant-comp $(TEST_WAV) out/shift-fc.wav
	out/stretch --time=0.8 --semitones=10 --formant-comp --formant=3 $(TEST_WAV) out/shift-fc-f3.wav
	out/stretch --time=0.8 --semitones=10 --formant-comp --formant=3 --formant-base=500 $(TEST_WAV) out/shift-fc-f3-fb500.wav
