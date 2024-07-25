# generate easy makefile

main: parameters/parameters.o parameters/readParameters.o main.o
	g++ -o main parameters/parameters.o parameters/readParameters.o main.o

parameters/parameters.o: parameters/parameters.cpp
	g++ -c parameters/parameters.cpp -I../include -I../include/nlohmann -o parameters/parameters.o

parameters/readParameters.o: parameters/readParameters.cpp
	g++ -c parameters/readParameters.cpp -I../include -I../include/nlohmann -o parameters/readParameters.o

main.o: main.cpp
	g++ -c main.cpp -I../include -I./parameters

clean:
	rm -f *.o main