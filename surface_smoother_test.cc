// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_smoother.h"

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

const char *kOriginalSurfaceFile = "poly_mesh.vtk";
const char *kSmoothedSurfaceFile = "smoothed_surfaces.vtk";

void smooth_surfaces_test() {
  printf("smooth_surfaces_test {\n");

  vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();

  reader->SetFileName(kOriginalSurfaceFile);
  reader->Update();

  vtkSmartPointer<vtkPolyData> surfaces =
      vtkSmartPointer<vtkPolyData>::New();
  surfaces->ShallowCopy(reader->GetOutput());

  SurfaceSmoother smoother;
  vtkPolyData *smoothed_surfaces = smoother.smooth_surfaces(surfaces, 0.5);

  vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();

  writer->SetFileName(kSmoothedSurfaceFile);
  writer->SetInputData(smoothed_surfaces);
  writer->Write();

  if (smoothed_surfaces) {
    smoothed_surfaces->Delete();
  }

  printf("} smooth_surfaces_test\n\n");
}

int main() {
  smooth_surfaces_test();

  return 0;
}
