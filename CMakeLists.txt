cmake_minimum_required(VERSION 3.10)
project(sputniPIC)
enable_language(CXX)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

if(NOT ARCH)
  set(ARCH "-march=native")
else()
  set(ARCH "-march=${ARCH}")
endif()


find_package(OpenMP)
find_package(MPI)
find_package(CUDA)
find_package(ParaView 5.8 REQUIRED COMPONENTS PythonCatalyst)

set(CXX_SRC
    src/mpi_comm.cpp
    src/Solvers.cpp
    src/Mover.cpp
    src/Smoothing.cpp
    src/InterpDensSpecies.cpp
    src/Particles.cpp
    src/EMfield.cpp
    src/IC.cpp
    src/Basic.cpp
    src/BC.cpp
    src/RW_IO.cpp
    src/MaxwellSolver.cpp
    src/ConfigFile.cpp
    src/Alloc.cpp
    src/InterpDensNet.cpp
    src/Grid.cpp
)

set(CUDA_SRC
    src/gpu/BC_gpu.cu
    src/gpu/cuda_helper.cpp
    src/gpu/EMfield_gpu.cpp
    src/gpu/Grid_gpu.cpp
    src/gpu/InterpDensNet_gpu.cu
    src/gpu/InterpDensSpecies_gpu.cpp  
    src/gpu/Particles_gpu.cu
)

if(ParaView_FOUND)
    add_library(catalystAdaptor STATIC src/Adaptor.cpp)
    target_include_directories(catalystAdaptor PRIVATE include)
    target_link_libraries(catalystAdaptor PRIVATE
        ParaView::Catalyst ParaView::PythonCatalyst VTK::CommonDataModel)
endif()

if(OpenMP_CXX_FOUND AND MPI_CXX_FOUND)
    add_library(sputniPIC STATIC ${CXX_SRC})
    target_link_libraries(sputniPIC PUBLIC
        OpenMP::OpenMP_CXX
        MPI::MPI_CXX)
    target_include_directories(sputniPIC PUBLIC 
        include
        ${MPI_INCLUDE_PATH})
    target_compile_options(sputniPIC PUBLIC -O3 ${ARCH} -ffast-math)
    add_executable(sputniPIC_CPU.out src/sputniPIC_CPU.cpp)
    target_link_libraries(sputniPIC_CPU.out sputniPIC)

    if (ParaView_FOUND)
        target_compile_definitions(sputniPIC PUBLIC -DUSE_CATALYST)
        target_link_libraries(sputniPIC PRIVATE catalystAdaptor)
    endif()

    if(CUDA_FOUND)
        enable_language(CUDA)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --forward-unknown-to-host-compiler")
        add_library(sputniPIC_gpu STATIC ${CUDA_SRC})
        target_link_libraries(sputniPIC_gpu PRIVATE
            sputniPIC
            OpenMP::OpenMP_CXX)
        target_include_directories(sputniPIC_gpu PUBLIC 
            ${CUDA_INCLUDE_DIRS}
            include)
        add_executable(sputniPIC_GPU.out src/sputniPIC_GPU.cpp)
        target_compile_definitions(sputniPIC_gpu PRIVATE -DUSE_GPU)
        target_compile_options(sputniPIC_gpu PUBLIC -O3 ${ARCH} -ffast-math)
        target_link_libraries(sputniPIC_GPU.out sputniPIC sputniPIC_gpu)

        target_compile_options(sputniPIC_gpu PRIVATE
            "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75>")
    endif()
endif()
