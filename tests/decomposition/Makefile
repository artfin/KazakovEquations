CCXX := g++
CXXFLAGS := -std=c++11 -O2

EIGEN := -I /usr/local/include/eigen3/
GSL := -lgsl -lgslcblas

main: ./main.cpp
	$(CCXX) $(CXXFLAGS) $(EIGEN) $< -o $@ $(GSL)
