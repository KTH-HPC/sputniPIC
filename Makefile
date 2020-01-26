# From https://x.momo86.net/?p=29

CXX=g++
CXXFLAGS=-DMEMCHECK -std=c++11 -I./include -O3 -g -Xcompiler "-fopenmp -Wall -Wno-unknown-pragmas"

NVCC=nvcc
# K420, K80, GTX1050Ti, V100
ARCH= -gencode arch=compute_30,code=sm_30 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70
NVCCFLAGS= -lineinfo -I./include $(ARCH) -std=c++11 -O3 -g -Xcompiler "-fopenmp -Wall -Wno-unknown-pragmas" --compiler-bindir=$(CXX)

SRCDIR=src
SRCS=$(shell find $(SRCDIR) -name '*.cu' -o -name '*.cpp')

OBJDIR=src
OBJS=$(subst $(SRCDIR),$(OBJDIR), $(SRCS))
OBJS:=$(subst .cpp,.o,$(OBJS))
OBJS:=$(subst .cu,.o,$(OBJS))

BIN := ./bin
TARGET=sputniPIC_GPU.out

all: dir $(BIN)/$(TARGET)

dir: ${BIN}
  
${BIN}:
	mkdir -p $(BIN)

$(BIN)/$(TARGET): $(OBJS)
	$(NVCC) $(NVCCFLAGS) $+ -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	[ -d $(OBJDIR) ] || mkdir $(OBJDIR)
	$(NVCC) $(CXXFLAGS) $< -c -o $@

clean:
	rm -rf $(OBJS)
	rm -rf $(TARGET)
