// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_extractor.h"

#include "surface_shader.h"
#include "surface_smoother.h"

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
// const char *kBasinFile = "data/basins_cell_8.vtk";
const char *kFTLEFile = "data/smoothed_ftle.vtk";
// const char *kBasinFile = "data/basins.vtk";
// const char *kFTLEFile = "data/ftle.vtk";
// const char *kBasinFile = "data/sphere_basins.vtk";
// const char *kFTLEFile = "data/sphere_ftle.vtk";
const char *kPolyDataFile = "poly_mesh.vtk";
const char *kColoredSurfaceFile = "colored_surfaces.vtk";

const int kNumberOfSmoothing = 50;  // 200 for brute-force mesh
const double kLambda = 1.0;

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

void extract_surfaces_with_regions_test() {
  printf("extract_surfaces_with_regions_test {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kBasinFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> basins =
      vtkSmartPointer<vtkStructuredPoints>::New();
  basins->ShallowCopy(reader->GetOutput());

  SurfaceExtractor extractor;
  vtkPolyData *positive_region = NULL;
  vtkPolyData *negative_region = NULL;
  extractor.extract_surfaces_with_regions(basins,
                                          &positive_region,
                                          &negative_region);

  int num_points = positive_region->GetNumberOfPoints();
  int num_cells = positive_region->GetNumberOfCells();

  printf("num_points = %d\n", num_points);
  printf("num_cells = %d\n", num_cells);

  SurfaceSmoother smoother;
  vtkPolyData *smoothed_surfaces = smoother.smooth_surfaces(
      positive_region, kLambda);

  for (int count = 1; count < kNumberOfSmoothing; count++) {
    printf("count = %d\n", count);

    vtkPolyData *new_smoothed_surfaces =
        smoother.smooth_surfaces(smoothed_surfaces, kLambda);
    smoothed_surfaces->Delete();

    smoothed_surfaces = new_smoothed_surfaces;
  }

  printf("Finished smoothing.\n");

  SurfaceShader shader;
  vtkPolyData *colored_surfaces = shader.get_colored_surfaces(
      smoothed_surfaces, negative_region);

  vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();

  writer->SetFileName(kColoredSurfaceFile);
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

  if (smoothed_surfaces) {
    smoothed_surfaces->Delete();
  }

  printf("} extract_surfaces_with_regions_test\n\n");
}

int main() {
  // extract_surfaces_test();
  extract_surfaces_with_regions_test();

  return 0;
}
