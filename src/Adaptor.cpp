#ifdef USE_CATALYST
#include "Adaptor.h"
#include "PrecisionTypes.h"

#include <iostream>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkImageData.h>

namespace
{
vtkCPProcessor* Processor = nullptr;
vtkImageData* VTKGrid = nullptr;
const char* InputName = "particles";

int _start_x;
int _nx;
int _dx;
int _start_y;
int _ny;
int _dy;
int _start_z;
int _nz;
int _dz;


//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
void UpdateVTKAttributes(vtkCPInputDataDescription* idd, FPfield ***Bx, FPfield ***By, FPfield ***Bz, FPinterp ***rhons0_array, FPinterp ***rhons1_array)
{
  // I am not sure whether we need to do this check 
  if (idd->IsFieldNeeded("B", vtkDataObject::POINT) == true)
  {
    // Create a VTK object representing magnetic field array

    // Get a reference to the grid's point data object.
    vtkPointData* vtk_point_data = VTKGrid->GetPointData();

    // We need to create a new VTK array object and attach it to the point data,
    // if it hasn't been done yet.
    if (vtk_point_data->GetNumberOfArrays() == 0)
    {
      vtkNew<vtkDoubleArray> field_array;
      field_array->SetName("B");
      field_array->SetNumberOfComponents(3);
      field_array->SetNumberOfTuples(static_cast<vtkIdType>(_nx * _ny * _nz));
      vtk_point_data->AddArray(field_array);
    }
    vtkDoubleArray* field_array = vtkDoubleArray::SafeDownCast(vtk_point_data->GetArray("B"));

    // Feed the data into VTK array. Since we don't know the memory layout of our B field data, 
    // we feed it point-by-point, in a very slow way

    // Array of grid's dimensions
    int* dims = VTKGrid->GetDimensions();

    // electric charge density for species 0 and 1
    vtkNew<vtkDoubleArray> rhons0{}, rhons1{};
    rhons0->SetName("RhoN_s0");
    rhons0->SetNumberOfComponents(1);
    rhons0->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
    vtk_point_data->AddArray(rhons0);

    rhons1->SetName("RhoN_s1");
    rhons1->SetNumberOfComponents(1);
    rhons1->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
    vtk_point_data->AddArray(rhons1);

    // Cycle over all VTK grid's points, get their indices and copy the data.
    // We want to have only one cycle over point's ID to efficiently use multi-threading.
    for (long p = 0; p < VTKGrid->GetNumberOfPoints(); ++p)
    {
      // Get cells's indices i, j , k
      unsigned long k = p / (dims[0] * dims[1]);
      unsigned long j = (p - k * dims[0] * dims[1]) / dims[0];
      unsigned long i = p - k * dims[0] * dims[1] - j * dims[0];

      // CAUTION!!! K should be always zero in the 2D case?
      field_array->SetComponent(p, 0, Bx[i][j][k]);
      field_array->SetComponent(p, 1, By[i][j][k]);
      field_array->SetComponent(p, 2, Bz[i][j][k]);

      rhons0->SetValue(p, rhons0_array[i][j][k]);
      rhons1->SetValue(p, rhons1_array[i][j][k]);
    }

    /// Fast way, if memry layout is correct.
    //velocityData->SetArray(const_cast<double*>(velocity.data()), static_cast<vtkIdType>(velocity.size()), 1);
  }
//  if (idd->IsFieldNeeded("collision", vtkDataObject::POINT) == true)
//  {
//    if (VTKGrid->GetPointData()->GetArray("collision") == nullptr)
//    {
//      // velocity array
//      vtkNew<vtkIntArray> collisionData;
//      collisionData->SetName("collision");
//      collisionData->SetNumberOfComponents(1);
//      collisionData->SetNumberOfTuples(static_cast<vtkIdType>(collisions.size()));
//      VTKGrid->GetPointData()->AddArray(collisionData);
//    }
//    vtkIntArray* collisionData =
//      vtkIntArray::SafeDownCast(VTKGrid->GetPointData()->GetArray("collision"));
//
//    collisionData->SetArray(const_cast<int*>(collisions.data()), static_cast<vtkIdType>(collisions.size()), 1);
//  }
}

//----------------------------------------------------------------------------
void BuildVTKDataStructures(vtkCPInputDataDescription *idd, FPfield ***Bx, FPfield ***By, FPfield ***Bz, FPinterp ***rhons0_array, FPinterp ***rhons1_array)
{
  // feed data to grid
  UpdateVTKAttributes(idd, Bx, By, Bz, rhons0_array, rhons1_array);
}
}

namespace Adaptor
{

//----------------------------------------------------------------------------
void Initialize(const char* script, const int start_x, const int start_y, const int start_z, \
                          const int nx, const int ny, const int nz, \
                          const double dx, const double dy, const double dz)
{
  if (Processor == NULL)
  {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
  }
  else
  {
    Processor->RemoveAllPipelines();
  }
  vtkNew<vtkCPPythonScriptPipeline> pipeline;
  pipeline->Initialize(script);
  Processor->AddPipeline(pipeline);

  _start_x = start_x;
  _nx = nx;
  _dx = dx;
  _start_y = start_y;
  _ny = ny;
  _dy = dy;
  _start_z = start_z;
  _nz = nz;
  _dz = dz;

  if (VTKGrid == NULL)
  {
    // The grid structure isn't changing so we only build it
    // the first time it's needed. If we needed the memory
    // we could delete it and rebuild as necessary.
    VTKGrid = vtkImageData::New();
    printf("%d %d %d %d\n", start_x, start_z, nx, nz);
    VTKGrid->SetExtent(start_x, start_x+nx-1, start_y, start_y+ny-1, start_z, start_z+nz-1);
    VTKGrid->SetSpacing(dx, dy, dz);
  }


}

//----------------------------------------------------------------------------
void Finalize()
{
  if (Processor)
  {
    Processor->Delete();
    Processor = NULL;
  }
  if (VTKGrid)
  {
    VTKGrid->Delete();
    VTKGrid = NULL;
  }
}

//----------------------------------------------------------------------------
void CoProcess(double time, unsigned int timeStep, FPfield ***Bx, FPfield ***By, FPfield ***Bz, FPinterp ***rhons0_array, FPinterp ***rhons1_array)
{
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput(InputName);
  dataDescription->SetTimeData(time, timeStep);

  if (Processor->RequestDataDescription(dataDescription) != 0)
  {
    vtkCPInputDataDescription* idd = dataDescription->GetInputDescriptionByName(InputName);
    BuildVTKDataStructures(idd, Bx, By, Bz, rhons0_array, rhons1_array);
    idd->SetGrid(VTKGrid);
    Processor->CoProcess(dataDescription);
  }
}
} // end of Catalyst namespace

#endif
