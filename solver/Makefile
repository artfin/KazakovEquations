.PHONY: all clean

CCXX := g++
CXXFLAGS := -std=c++11 -O2 -Wall -g -Wextra

EIGEN := -I/usr/local/include/eigen3/
GSL := -lgsl -lgslcblas

BUILDDIR := ./build/

all: main 

main: $(addprefix $(BUILDDIR), $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
	$(CCXX) $(CXXFLAGS) $^ -o $@ $(GSL)

$(BUILDDIR)%.o: %.cpp
	@echo ">> (g++) Compiling $<...";
	@$(CCXX) $(CXXFLAGS) -c -MD $< -o $@ $(EIGEN)

clean:
	@rm -f $(BUILDDIR)*.d
	@rm -f $(BUILDDIR)*.o
	@rm -f main 
