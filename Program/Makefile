UNAME_S := $(shell uname -s)
CCC = g++
TARGETDIR=.
CCFLAGS = -O3 -Wall -std=c++11
ifeq ($(UNAME_S),Linux)
    CCFLAGS += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
    CCFLAGS += -openmp
endif

debug: CCFLAGS += -g
release: CCFLAGS += -O3 -march=native

debug: runTest
release: runTest
 
OBJS2 = \
        $(TARGETDIR)/BRKGA_QL.o \
        $(TARGETDIR)/Problem.o 

$(TARGETDIR)/runTest: $(OBJS2)
	$(CCC) $(CCFLAGS) -o $(TARGETDIR)/runTest $(OBJS2)
    
$(TARGETDIR)/Problem.o: Problem.h Problem.cpp
	$(CCC) $(CCFLAGS) -c Problem.cpp -o $(TARGETDIR)/Problem.o

$(TARGETDIR)/BRKGA_QL.o: BRKGA_QL.h BRKGA_QL.cpp
	$(CCC) $(CCFLAGS) -c BRKGA_QL.cpp -o $(TARGETDIR)/BRKGA_QL.o

test: runTest
	./runTest

clean:
	$(RM) \
	$(TARGETDIR)/BRKGA_QL.o \
	$(TARGETDIR)/Problem.o
    
rebuild:
	$(MAKE) clean
	$(MAKE) debug