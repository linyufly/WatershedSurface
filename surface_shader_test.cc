// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_shader.h"

#include <vtkCellArray.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

const char *kOriginalSurfaceFile = "smoothed_surfaces.vtk";
const char *kShadedSurfaceFile = "shaded_surfaces.vtk";

void shade_by_region_test() {
  printf("shaded_by_region_test {\n");

  vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();

  reader->SetFileName(kOriginalSurfaceFile);
  reader->Update();

  vtkSmartPointer<vtkPolyData> surfaces =
      vtkSmartPointer<vtkPolyData>::New();
  surfaces->DeepCopy(reader->GetOutput());

  vtkPolyData *positive_region = NULL;
  vtkPolyData *negative_region = NULL;

  SurfaceShader shader;
  shader.shade_by_region(surfaces, &positive_region, &negative_region);

  vtkPolyData *colored_surfaces = shader.get_colored_surfaces(
      positive_region, negative_region);
  
  vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(kShadedSurfaceFile);
  writer->SetInputData(colored_surfaces);
  writer->Write();

  if (positive_region) {
    positive_region->Delete();
  }

  if (negative_region) {
    negative_region->Delete();
  }

  if (colored_surfaces) {
    colored_surfaces->Delete();
  }

  printf("} shade_by_region_test\n\n");
}

int main() {
  shade_by_region_test();

  return 0;
}
