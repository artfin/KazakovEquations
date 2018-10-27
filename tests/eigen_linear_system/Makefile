CCXX := g++
CXXFLAGS := -O2 -Wall -std=c++11 

LINK_EIGEN := -I/usr/local/include/eigen3/

main: ./main.o
	$(CCXX) $(CXXFLAGS) $^ -o $@

./main.o: ./main.cpp
	$(CCXX) $(CXXFLAGS) -c $< -o $@ $(LINK_EIGEN)
