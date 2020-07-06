VERSION=GPU

CXX=mpicxx
CXXFLAGS=-std=c++11 -I./include -O3 -g -fopenmp

NVCC=nvcc
ARCH=-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75
NVCCFLAGS=-DMEMCHECK -DUSE_GPU -lineinfo -I./include $(ARCH) -std=c++11 -O3 -g -Xcompiler "-fopenmp -Wno-unknown-pragmas" --compiler-bindir=$(CXX)

# Default to use host compiler and flags
COMPILER=$(CXX)
FLAGS=$(CXXFLAGS)
TARGET=sputniPIC.out
SRCDIR=src

# Check go GPU or CPU path
ifeq ($(VERSION), GPU)
    FILES=$(shell find $(SRCDIR) -name '*.cu' -o -name '*.cpp')
    SRCS_1=$(subst $(SRCDIR)/sputniPIC_CPU.cpp,,${FILES})
    SRCS=$(subst $(SRCDIR)/sputniPIC_GPU_UVM.cpp,,${SRCS_1})

    COMPILER=$(NVCC)
    COMPILER_FLAGS=$(NVCCFLAGS)
else ifeq ($(VERSION), GPU_UVM)
    FILES=$(shell find $(SRCDIR) -name '*.cu' -o -name '*.cpp')
    SRCS_1=$(subst $(SRCDIR)/sputniPIC_CPU.cpp,,${FILES})
    SRCS=$(subst $(SRCDIR)/sputniPIC_GPU.cpp,,${SRCS_1})

    COMPILER=$(NVCC)
    NVCCFLAGS+=-DCUDA_UVM #-DCUDA_UVM_PREFETCH
    COMPILER_FLAGS=$(NVCCFLAGS)
else
    FILES=$(shell find $(SRCDIR) -path $(SRCDIR)/gpu -prune -o -name '*.cpp' -print)
    SRCS_1=$(subst $(SRCDIR)/sputniPIC_GPU.cpp,,${FILES})
    SRCS=$(subst $(SRCDIR)/sputniPIC_GPU_UVM.cpp,,${SRCS_1})

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
	$(COMPILER) $(COMPILER_FLAGS) $+ -o $@

# GPU objects
$(SRCDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

# Host objects
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILER) $(COMPILER_FLAGS) $< -c -o $@

clean:
	rm -rf */*.o */*/*.o
	rm -rf $(OBJS)
	rm -rf $(BIN)/*.out
