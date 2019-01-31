CCXX := g++
CXXFLAGS := -std=c++11 -O2 -Wall

INCLUDE_EIGEN := -isystem /usr/local/include/eigen3
INCLUDE_GSL := -I/usr/include/gsl/

LINK_GSL := -lgsl -lgslcblas

all: main

main: ./main.o
	$(CCXX) $(CXXFLAGS) $^ -o $@ $(LINK_GSL)

main.o: ./main.cpp
	$(CCXX) $(CXXFLAGS) -c $< -o $@ $(INCLUDE_EIGEN) $(INCLUDE_GSL) 
