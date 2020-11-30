VERSION=CPU

CXX=mpicxx
CXXFLAGS=-std=c++11 -I./include -O3 -g -fopenmp

NVCC=nvcc
ARCH=-gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70
NVCCFLAGS=-USE_CATALYST -DMEMCHECK -DUSE_GPU -lineinfo -I./include $(ARCH) -std=c++11 -O3 -g -Xcompiler "-fopenmp -Wno-unknown-pragmas -isystem /usr/local/include/paraview-5.8 -isystem /home/steven/anaconda3/include/python3.7m -fPIC" --compiler-bindir=$(CXX)

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
	$(COMPILER) -DUSE_CATALYST --compiler-bindir=$(CXX) $+ -o $@ -Xcompiler "-fopenmp -std=c++11" --linker-options="--no-as-needed -W,--unresolved-symbols=ignore-in-shared-libs -W,-rpath /usr/local/openmpi/lib --enable-new-dtags -pthread -W,-rpath,/usr/local/lib:/home/steven/anaconda3/lib:/usr/local/openmpi/lib /usr/local/lib/libvtkWrappingTools-pv5.8.so.5.8 /usr/local/lib/libvtkWebCore-pv5.8.so.5.8 /usr/local/lib/libvtkWebGLExporter-pv5.8.so.5.8 /usr/local/lib/libvtkViewsContext2D-pv5.8.so.5.8 /usr/local/lib/libvtkTestingRendering-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingRayTracing-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingSceneGraph-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingVolumeAMR-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingVolumeOpenGL2-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingParallel-pv5.8.so.5.8 /usr/local/lib/libvtkImagingMath-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingMatplotlib-pv5.8.so.5.8 /usr/local/lib/libvtkPythonInterpreter-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingLabel-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingLICOpenGL2-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingContextOpenGL2-pv5.8.so.5.8 /usr/local/lib/libvtkParallelMPI4Py-pv5.8.so.5.8 /usr/local/lib/libvtkIOXdmf2-pv5.8.so.5.8 /usr/local/lib/libvtkxdmf2-pv5.8.so.5.8 /usr/local/lib/libvtkIOVeraOut-pv5.8.so.5.8 /usr/local/lib/libvtkIOVPIC-pv5.8.so.5.8 /usr/local/lib/libvtkvpic-pv5.8.so.5.8 /usr/local/lib/libvtkIOTecplotTable-pv5.8.so.5.8 /usr/local/lib/libvtkIOTRUCHAS-pv5.8.so.5.8 /usr/local/lib/libvtkIOSegY-pv5.8.so.5.8 /usr/local/lib/libvtkIOParallelXML-pv5.8.so.5.8 /usr/local/lib/libvtkIOParallelNetCDF-pv5.8.so.5.8 /usr/local/lib/libvtkIOParallelLSDyna-pv5.8.so.5.8 /usr/local/lib/libvtkIOLSDyna-pv5.8.so.5.8 /usr/local/lib/libvtkIOParallelExodus-pv5.8.so.5.8 /usr/local/lib/libvtkIOPLY-pv5.8.so.5.8 /usr/local/lib/libvtkIOPIO-pv5.8.so.5.8 /usr/local/lib/libvtkIOOggTheora-pv5.8.so.5.8 /usr/local/lib/libvtktheora-pv5.8.so.5.8 /usr/local/lib/libvtkogg-pv5.8.so.5.8 /usr/local/lib/libvtkIONetCDF-pv5.8.so.5.8 /usr/local/lib/libvtkIOParallel-pv5.8.so.5.8 /usr/local/lib/libvtkIOMPIImage-pv5.8.so.5.8 /usr/local/lib/libvtkIOInfovis-pv5.8.so.5.8 /usr/local/lib/libvtkIOImport-pv5.8.so.5.8 /usr/local/lib/libvtkIOH5part-pv5.8.so.5.8 /usr/local/lib/libvtkh5part-pv5.8.so.5.8 /usr/local/lib/libvtkIOGeometry-pv5.8.so.5.8 /usr/local/lib/libvtkIOFFMPEG-pv5.8.so.5.8 /usr/local/lib/libvtkIOVideo-pv5.8.so.5.8 /usr/local/lib/libvtkIOMovie-pv5.8.so.5.8 /usr/local/lib/libvtkIOExportGL2PS-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingGL2PSOpenGL2-pv5.8.so.5.8 /usr/local/lib/libvtkgl2ps-pv5.8.so.5.8 /usr/local/lib/libvtkIOExodus-pv5.8.so.5.8 /usr/local/lib/libvtkIOEnSight-pv5.8.so.5.8 /usr/local/lib/libvtkIOCityGML-pv5.8.so.5.8 /usr/local/lib/libvtkIOAsynchronous-pv5.8.so.5.8 /usr/local/lib/libvtkIOAMR-pv5.8.so.5.8 /usr/local/lib/libvtkViewsCore-pv5.8.so.5.8 /usr/local/lib/libvtkGUISupportQt-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersPython-pv5.8.so.5.8 /usr/local/lib/libvtkWrappingPythonCore-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersProgrammable-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersPoints-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersParallelVerdict-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersVerdict-pv5.8.so.5.8 /usr/local/lib/libvtkverdict-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersParallelStatistics-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersParallelGeometry-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersParallelFlowPaths-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersParallelDIY2-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersHyperTree-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersGeneric-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersFlowPaths-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersAMR-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersParallelMPI-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersParallel-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersTexture-pv5.8.so.5.8 /usr/local/lib/libvtkDomainsChemistryOpenGL2-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingOpenGL2-pv5.8.so.5.8 /usr/local/lib/libvtkglew-pv5.8.so.5.8 /usr/local/lib/libvtkDomainsChemistry-pv5.8.so.5.8 /usr/local/lib/libvtkChartsCore-pv5.8.so.5.8 /usr/local/lib/libvtkAcceleratorsVTKm-pv5.8.so.5.8 /usr/local/lib/libospray.so.1.8.5 /usr/local/lib/libospray_common.so.1.8.5 /usr/lib/x86_64-linux-gnu/libtbb.so.2 /usr/lib/x86_64-linux-gnu/libtbbmalloc.so.2 /usr/local/lib/libvtkexodusII-pv5.8.so.5.8 /usr/local/lib/libvtknetcdf-pv5.8.so.5.8 /usr/local/lib/libvtklibxml2-pv5.8.so.5.8 /usr/local/lib/libvtkhdf5_hl-pv5.8.so.5.8 /usr/local/lib/libvtkhdf5-pv5.8.so.5.8 /usr/local/lib/libvtkIOExport-pv5.8.so.5.8 /usr/local/lib/libvtklibharu-pv5.8.so.5.8 /usr/local/lib/libvtkjsoncpp-pv5.8.so.5.8 /usr/local/lib/libvtkInteractionWidgets-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingVolume-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingAnnotation-pv5.8.so.5.8 /usr/local/lib/libvtkImagingColor-pv5.8.so.5.8 /usr/local/lib/libvtkImagingHybrid-pv5.8.so.5.8 /usr/local/lib/libvtkImagingGeneral-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersHybrid-pv5.8.so.5.8 /usr/local/lib/libvtkInteractionStyle-pv5.8.so.5.8 /home/steven/anaconda3/lib/libQt5Widgets.so.5.12.9 /home/steven/anaconda3/lib/libQt5Gui.so.5.12.9 /home/steven/anaconda3/lib/libQt5Core.so.5.12.9 /home/steven/anaconda3/lib/libpython3.7m.so /usr/local/lib/libvtkFiltersModeling-pv5.8.so.5.8 /usr/lib/x86_64-linux-gnu/libXt.so /usr/lib/x86_64-linux-gnu/libX11.so /usr/lib/x86_64-linux-gnu/libICE.so /usr/lib/x86_64-linux-gnu/libSM.so /usr/lib/x86_64-linux-gnu/libGLX.so /usr/lib/x86_64-linux-gnu/libOpenGL.so /usr/local/lib/libvtkInfovisCore-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersExtraction-pv5.8.so.5.8 /usr/local/lib/libvtkParallelDIY-pv5.8.so.5.8 /usr/local/lib/libvtkIOXML-pv5.8.so.5.8 /usr/local/lib/libvtkIOXMLParser-pv5.8.so.5.8 /usr/local/lib/libvtkexpat-pv5.8.so.5.8 /usr/local/lib/libvtkParallelMPI-pv5.8.so.5.8 /usr/local/openmpi/lib/libmpi.so /usr/local/lib/libvtkParallelCore-pv5.8.so.5.8 /usr/local/lib/libvtkIOLegacy-pv5.8.so.5.8 /usr/local/lib/libvtkIOCore-pv5.8.so.5.8 /usr/local/lib/libvtklzma-pv5.8.so.5.8 /usr/local/lib/libvtklz4-pv5.8.so.5.8 /usr/local/lib/libvtkdoubleconversion-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersStatistics-pv5.8.so.5.8 /usr/local/lib/libvtkImagingFourier-pv5.8.so.5.8 /usr/local/lib/libvtkImagingSources-pv5.8.so.5.8 /usr/local/lib/libvtkIOImage-pv5.8.so.5.8 /usr/local/lib/libvtkpng-pv5.8.so.5.8 /usr/local/lib/libvtkDICOMParser-pv5.8.so.5.8 /usr/local/lib/libvtkmetaio-pv5.8.so.5.8 /usr/local/lib/libvtktiff-pv5.8.so.5.8 /usr/local/lib/libvtkjpeg-pv5.8.so.5.8 /usr/local/lib/libvtkpugixml-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingContext2D-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingFreeType-pv5.8.so.5.8 /usr/local/lib/libvtkfreetype-pv5.8.so.5.8 /usr/local/lib/libvtkzlib-pv5.8.so.5.8 /usr/local/lib/libvtkRenderingCore-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersSources-pv5.8.so.5.8 /usr/local/lib/libvtkCommonColor-pv5.8.so.5.8 /usr/local/lib/libvtkm_filter-pv5.8.a /usr/local/lib/libvtkm_worklet-pv5.8.a /usr/local/lib/libvtkm_cont-pv5.8.a /usr/local/lib/libvtkImagingCore-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersGeometry-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersGeneral-pv5.8.so.5.8 /usr/local/lib/libvtkCommonComputationalGeometry-pv5.8.so.5.8 /usr/local/lib/libvtkFiltersCore-pv5.8.so.5.8 /usr/local/lib/libvtkCommonExecutionModel-pv5.8.so.5.8 /usr/local/lib/libvtkCommonDataModel-pv5.8.so.5.8 /usr/local/lib/libvtkCommonSystem-pv5.8.so.5.8 /usr/local/lib/libvtkCommonMisc-pv5.8.so.5.8 /usr/local/lib/libvtkCommonTransforms-pv5.8.so.5.8 /usr/local/lib/libvtkCommonMath-pv5.8.so.5.8 /usr/local/lib/libvtkCommonCore-pv5.8.so.5.8 /usr/local/lib/libvtkloguru-pv5.8.so.5.8 /usr/local/lib/libvtksys-pv5.8.so.5.8 /usr/local/lib/libvtkPVCatalyst-pv5.8.so /usr/local/lib/libvtkPVPythonCatalyst-pv5.8.so -ldl -pthread -W,-rpath-link,/home/steven/anaconda3/lib" 

# GPU objects
$(SRCDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@

# Host objects
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILER) $(COMPILER_FLAGS) $< -c -o $@

clean:
	rm -rf $(OBJS)
	rm -rf $(BIN)/*.out
