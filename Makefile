all: out/stretch

out/stretch: signalsmith-stretch.h cmd/main.cpp cmd/util/*.h
	mkdir -p out
	g++ cmd/main.cpp -o out/stretch -std=c++11 -Ofast -g

clean:
	rm -rf out