// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_extractor.h"

#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkIndent.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

#include <cstdio>
#include <cstdlib>

#include <iostream>

// const char *kBasinFile = "data/smoothed_basins.vtk";
// const char *kBasinFile = "../Watershed/basin_index.vtk";
const char *kBasinFile = "../Watershed/filtered_basin.vtk";
const char *kFTLEFile = "data/smoothed_ftle.vtk";
// const char *kBasinFile = "data/basins.vtk";
// const char *kFTLEFile = "data/ftle.vtk";
// const char *kBasinFile = "data/sphere_basins.vtk";
// const char *kFTLEFile = "data/sphere_ftle.vtk";
const char *kPolyDataFile = "poly_mesh.vtk";

void extract_surfaces_test() {
  printf("extract_surfaces_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kBasinFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> basins =
      vtkSmartPointer<vtkStructuredPoints>::New();
  basins->ShallowCopy(reader->GetOutput());

  reader->SetFileName(kFTLEFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> ftle =
      vtkSmartPointer<vtkStructuredPoints>::New();
  ftle->ShallowCopy(reader->GetOutput());

  SurfaceExtractor extractor;
  vtkPolyData *surfaces = extractor.extract_surfaces(ftle, basins);

  int num_points = surfaces->GetNumberOfPoints();
  int num_cells = surfaces->GetNumberOfCells();

  printf("num_points = %d\n", num_points);
  printf("num_cells = %d\n", num_cells);

  vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();

  writer->SetFileName(kPolyDataFile);
  writer->SetInputData(surfaces);
  writer->Write();

  if (surfaces) {
    surfaces->Delete();
  }

  printf("} extract_surfaces_test\n\n");
}

int main() {
  extract_surfaces_test();

  return 0;
}
