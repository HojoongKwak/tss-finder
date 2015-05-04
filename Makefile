all:
	g++ -std=c++11 -o bin/find_tss src/find_tss.cpp src/smooth.cpp src/bg.cpp src/arguments.cpp src/tabfile.cpp -I./include -pthread

