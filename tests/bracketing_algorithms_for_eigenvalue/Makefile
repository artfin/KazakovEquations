CCXX := g++
CXXFLAGS := -std=c++11 -O2

INCLUDE_EIGEN := -I /usr/local/include/eigen3/
LINK_GSL := -lgsl -lgslcblas

main: ./main.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $< -o $@ $(LINK_GSL)
