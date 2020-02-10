# From https://x.momo86.net/?p=29

VERSION=CPU

CXX=g++
CXXFLAGS=-DMEMCHECK -std=c++11 -I./include -O3 -g -fopenmp

NVCC=nvcc
ARCH= -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70
NVCCFLAGS= -DMEMCHECK -DUSE_GPU -lineinfo -I./include $(ARCH) -std=c++11 -O3 -g -Xcompiler "-fopenmp -Wall -Wno-unknown-pragmas" --compiler-bindir=$(CXX)

COMPILER=$(CXX)
FLAGS=$(NVCCFLAGS)

SRCDIR=src
ifeq ($(VERSION), CPU)
    FILES=$(shell find $(SRCDIR) -path $(SRCDIR)/gpu -prune -o -name '*.cpp' -print)
#    FILES=$(shell find $(SRCDIR) -maxdepth 1 -name '*.cpp')
    SRCS=$(subst $(SRCDIR)/sputniPIC_GPU.cpp,,${FILES})

    COMPILER=$(CXX)
    COMPILER_FLAGS=$(CXXFLAGS)

    TARGET=sputniPIC_CPU.out
else
    FILES=$(shell find $(SRCDIR) -name '*.cu' -o -name '*.cpp')
    SRCS=$(subst $(SRCDIR)/sputniPIC_CPU.cpp,,${FILES})

    COMPILER=$(NVCC)
    COMPILER_FLAGS=$(NVCCFLAGS)

    TARGET=sputniPIC_GPU.out
endif

OBJDIR=src
OBJS=$(subst $(SRCDIR),$(OBJDIR), $(SRCS))
OBJS:=$(subst .cpp,.o,$(OBJS))
OBJS:=$(subst .cu,.o,$(OBJS))

BIN := ./bin

all: dir $(BIN)/$(TARGET)

dir: ${BIN}
  
${BIN}:
	mkdir -p $(BIN)

$(BIN)/$(TARGET): $(OBJS)
	$(COMPILER) $(COMPILER_FLAGS) $+ -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	[ -d $(OBJDIR) ] || mkdir $(OBJDIR)
	$(COMPILER) $(COMPILER_FLAGS) $< -c -o $@

clean:
	rm -rf $(OBJS)
	rm -rf $(BIN)/*.out
