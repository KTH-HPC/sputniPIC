VERSION=GPU

CXX=mpicxx

HDF5_LIBS=-L/usr/local/lib -L/home/users/project/SHARED/libaoi -L/home/users/project/SHARED/numerical-object-storage/src -laoi -lobjects #-lhdf5
HDF5_CFLAGS=-I/usr/local/include
CXXFLAGS=-std=c++11 \
	$(HDF5_CFLAGS) \
	-I/home/users/project/SHARED/protoc-c/include \
	-I/home/users/project/SHARED/numerical-object-storage/include \
	-I/home/users/project/SHARED/libaoi \
	-I./include -O3 -g -fopenmp -Wall -Wno-unknown-pragmas -DUSE_MERO


NVCC=nvcc --forward-unknown-to-host-compiler

ARCH=#-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70
NVCCFLAGS=-DMEMCHECK -DUSE_GPU -lineinfo $(CXXFLAGS) $(ARCH) -std=c++11 --compiler-bindir=$(CXX)

# Default to use host compiler and flags
COMPILER=$(CXX)
FLAGS=$(CXXFLAGS)
SRCDIR=src

# Check go GPU or CPU path
ifeq ($(VERSION), GPU)
    TARGET=sputniPIC_GPU.out
    FILES=$(shell find $(SRCDIR) -name '*.cu' -o -name '*.cpp')
    SRCS=$(subst $(SRCDIR)/sputniPIC_CPU.cpp,,${FILES})

    COMPILER=$(NVCC)
    COMPILER_FLAGS=$(NVCCFLAGS)
else
    TARGET=sputniPIC_CPU.out
    FILES=$(shell find $(SRCDIR) -path $(SRCDIR)/gpu -prune -o -name '*.cpp' -print)
    SRCS=$(subst $(SRCDIR)/sputniPIC_GPU.cpp,,${FILES})

    COMPILER=$(CXX)
    COMPILER_FLAGS=$(CXXFLAGS)
endif

# Generate list of objects
OBJDIR=src
OBJS=$(subst $(SRCDIR),$(OBJDIR), $(SRCS))
OBJS:=$(subst .cpp,.o,$(OBJS))
OBJS:=$(subst .cu,.o,$(OBJS))

BIN := ./bin

all: dir $(BIN)/$(TARGET)

dir: ${BIN}

# Create bin folder
${BIN}:
	mkdir -p $(BIN)

# Binary linkage
$(BIN)/$(TARGET): $(OBJS)
	$(COMPILER) $(COMPILER_FLAGS) $+ -o $@ $(HDF5_LIBS)

# GPU objects
$(SRCDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

# Host objects
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILER) $(COMPILER_FLAGS) $< -c -o $@

clean:
	rm -rf $(OBJS)
	rm -rf $(BIN)/*.out
