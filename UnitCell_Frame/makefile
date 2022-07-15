CC=/usr/local/Cellar/gcc/10.2.0/bin/g++-10 -std=c++20 -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3

run : test.o atom.o cell.o
	$(CC) -o run test.o cell.o atom.o

atom.o : atom.cpp atom.hpp
	$(CC) -c atom.cpp

cell.o : cell.cpp cell.hpp
	$(CC) -c cell.cpp

test.o : test.cpp
	$(CC) -c test.cpp


clean : 
	rm -rf *.o
