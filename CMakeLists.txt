cmake_minimum_required(VERSION 3.10)
project(sputniPIC)
enable_language(CXX)

# ParaView adaptor is still under development!
option(Build_ParaView "Build ParaView support" OFF)
option(BUILD_GPU "Build GPU binary" ON)
option(BUILD_CPU "Build CPU binary" ON)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# improve optimization
if(NOT ARCH)
  set(ARCH "")
else()
  set(ARCH "-march=${ARCH}")
endif()

# find mandatory libraries
find_package(OpenMP)
find_package(MPI)
find_package(HDF5 COMPONENTS C)

# source code for sputniPIC
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

# source code for sputniPIC GPU
set(CUDA_SRC
    src/gpu/BC_gpu.cu
    src/gpu/cuda_helper.cpp
    src/gpu/EMfield_gpu.cpp
    src/gpu/Grid_gpu.cpp
    src/gpu/InterpDensNet_gpu.cu
    src/gpu/InterpDensSpecies_gpu.cpp  
    src/gpu/Particles_gpu.cu
)

# build ParaView Catalyst adaptor if paraview found
if(Build_ParaView)
    find_package(ParaView 5.8 REQUIRED COMPONENTS PythonCatalyst)
    if (ParaView_FOUND)
        add_library(catalystAdaptor STATIC src/Adaptor.cpp)
        target_include_directories(catalystAdaptor PRIVATE include ${PARAVIEW_INCLUDE_PATH})
        target_link_libraries(catalystAdaptor PRIVATE
            ParaView::Catalyst ParaView::PythonCatalyst VTK::CommonDataModel)
        target_compile_definitions(catalystAdaptor PUBLIC -DUSE_CATALYST)
        message(STATUS "Building ParaView support")
    endif()
endif()

# if all mandatory libraries are found
if(OpenMP_CXX_FOUND AND MPI_CXX_FOUND AND HDF5_FOUND)
    # create sputniPIC core
    add_library(sputniPIC STATIC ${CXX_SRC})

    # link sputniPIC with libraries
    target_link_libraries(sputniPIC PRIVATE
        ${HDF5_C_LIBRARIES}
        /home/users/project/SHARED/libaoi/libaoi.so
        /home/users/project/SHARED/numerical-object-storage/src/libobjects.so
        OpenMP::OpenMP_CXX
        MPI::MPI_CXX)

    # include header paths for sputniPIC
    target_include_directories(sputniPIC PRIVATE
        ${HDF5_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
        /home/users/project/SHARED/protoc-c/include
        /home/users/project/SHARED/libaoi
        /home/users/project/SHARED/numerical-object-storage/include)
    target_include_directories(sputniPIC PUBLIC
        include)

    # compile with highest optimization
    target_compile_options(sputniPIC PUBLIC -O3 ${ARCH} -ffast-math -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    target_compile_definitions(sputniPIC PUBLIC -DUSE_MERO -DVERIFY_MERO)

    # build CPU executable
    if (BUILD_CPU)
        add_executable(sputniPIC_CPU.out src/sputniPIC_CPU.cpp)
        target_link_libraries(sputniPIC_CPU.out PRIVATE sputniPIC)
        # link to Catalyst adaptor if adaptor is built
        if (ParaView_FOUND AND Build_ParaView)
            target_link_libraries(sputniPIC_CPU.out PRIVATE catalystAdaptor)
            target_compile_definitions(sputniPIC_CPU.out PUBLIC -DUSE_CATALYST)
        endif()
        message(STATUS "Building CPU binary")
    endif()

    # build GPU executable
    if(BUILD_GPU)
        # build only if CUDA is found
        find_package(CUDA)
        if(CUDA_FOUND)
            # enable CUDA and forward any unknown flags to host compiler
            enable_language(CUDA)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --forward-unknown-to-host-compiler")

            # create sputniPIC core for GPU
            add_library(sputniPIC_gpu STATIC ${CUDA_SRC})

            # link to sputniPIC core
            target_link_libraries(sputniPIC_gpu PRIVATE
                sputniPIC
                OpenMP::OpenMP_CXX)

            # include header paths to sputniPIC GPU and CUDA
            target_include_directories(sputniPIC_gpu PUBLIC 
                ${CUDA_INCLUDE_DIRS}
                include)

            # create GPU executable
            add_executable(sputniPIC_GPU.out src/sputniPIC_GPU.cpp)

            # enable all GPU parts of code and optimize
            target_compile_definitions(sputniPIC_gpu PUBLIC -DUSE_GPU)
            target_compile_options(sputniPIC_gpu PUBLIC -O3 ${ARCH} -ffast-math -Wall -Wextra -pedantic -Wno-unknown-pragmas)

            # support compute 3.5 onwards
            set_property(TARGET sputniPIC_gpu PROPERTY CUDA_ARCHITECTURES 35 37 50 60 61 70 75 80)

            # link to Catalyst adaptor if adaptor is built
            if (ParaView_FOUND AND Build_ParaView)
                target_link_libraries(sputniPIC_gpu PUBLIC catalystAdaptor)
                target_compile_definitions(sputniPIC_gpu PUBLIC -DUSE_CATALYST)
            endif()

            # link sputniPIC GPU executable to core GPU library
            target_link_libraries(sputniPIC_GPU.out sputniPIC_gpu)

            message(STATUS "Building CPU binary")
        endif()
    endif()
endif()
