VERSION=GPU

CXX=mpic++
CXXFLAGS=-std=c++11 -I./include -O3 -g -fopenmp -Wall

# Get mpi compiler wrapper arguments (without host compiler arg)
MPIFLAGS=$(shell $(CXX) --show | cut -d' ' -f2-)
# NVCC has a bug with -pthreads, pass to underlying compiler if used
MPI_COMPILE_FLAGS=$(subst -pthread,-Xcompiler="-pthread",$(patsubst -l%,,$(MPIFLAGS)))
MPI_LINK_FLAGS=$(subst -pthread,-Xcompiler="-pthread",$(patsubst -I%,,$(MPIFLAGS)))

NVCC=nvcc
ARCH=-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70
NVCCFLAGS=-DMEMCHECK -DUSE_GPU -lineinfo -I./include $(ARCH) -std=c++11 -O3 -g -Xcompiler "-fopenmp -Wall -Wno-unknown-pragmas" --compiler-bindir=$(CXX) $(MPI_COMPILE_FLAGS)


# Default to use host compiler and flags
COMPILER=$(CXX)
FLAGS=$(CXXFLAGS)
TARGET=sputniPIC.out
SRCDIR=src

# Check go GPU or CPU path

ifeq ($(VERSION), GPU)
    FILES=$(shell find $(SRCDIR) -name '*.cu' -o -name '*.cpp')
    SRCS=$(subst $(SRCDIR)/sputniPIC_CPU.cpp,,${FILES})

    COMPILER=$(NVCC)
    COMPILER_FLAGS=$(NVCCFLAGS)
else
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
	$(COMPILER) $(COMPILER_FLAGS) $(MPI_LINK_FLAGS) $+ -o $@

# GPU objects
$(SRCDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

# Host objects
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILER) $(COMPILER_FLAGS) $< -c -o $@

clean:
	rm -rf $(OBJS)
	rm -rf $(BIN)/*.out
