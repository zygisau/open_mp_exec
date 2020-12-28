CC=g++

all: main.cpp
	${CC} -std=c++17 -fopenmp -g main.cpp -o output
	./output

release: main.cpp
	${CC} -std=c++17 -O3 -fopenmp main.cpp -o output
	./output

nomp: main.cpp
	${CC} -std=c++17 main.cpp -o output
	./output
