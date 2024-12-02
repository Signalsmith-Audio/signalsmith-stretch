all: out/stretch

out/stretch: ../signalsmith-stretch.h main.cpp util/*.h
	mkdir -p out
	g++ main.cpp -o out/stretch -std=c++11 -Ofast -g

clean:
	rm -rf out