CC=g++

all: main.cpp
	${CC} -std=c++17 -fopenmp -g main.cpp -o output
	./output

release: main.cpp
	${CC} -std=c++17 -fopenmp main.cpp -o output
	./output