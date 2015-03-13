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
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkIndent.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

#include <cstdio>
#include <cstdlib>

#include <iostream>

// const char *kBasinFile = "data/smoothed_basins.vtk";
// const char *kBasinFile = "../Watershed/basin_index.vtk";
// const char *kBasinFile = "../Watershed/filtered_basin.vtk";
// const char *kBasinFile = "data/basins_cell_8.vtk";

// const char *kBasinFile = "two_spheres_basins.vtk";
// const char *kFTLEFile = "two_spheres_ftle.vtk";

const char *kBasinFile = "one_sphere_basins.vtk";
const char *kFTLEFile = "one_sphere_ftle.vtk";

// const char *kFTLEFile = "data/smoothed_ftle.vtk";
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
/*
void extract_surfaces_with_regions_test_2() {
  printf("extract_surfaces_with_regions_test_2 {\n");

  vtkSmartPointer<vtkStructuredPointsReader> reader =
      vtkSmartPointer<vtkStructuredPointsReader>::New();

  reader->SetFileName(kBasinFile);
  reader->Update();

  vtkSmartPointer<vtkStructuredPoints> basins =
      vtkSmartPointer<vtkStructuredPoints>::New();
  basins->ShallowCopy(reader->GetOutput());

  SurfaceExtractor extractor;
  vtkPolyData *surfaces = extractor.extract_surfaces_with_regions(
      NULL, basins);

  int num_points = surfaces->GetNumberOfPoints();
  int num_cells = surfaces->GetNumberOfCells();

  printf("num_points = %d\n", num_points);
  printf("num_cells = %d\n", num_cells);

  // vtkPolyData *smoothed_surfaces = vtkPolyData::New();
  // smoothed_surfaces->DeepCopy(surfaces);

  SurfaceSmoother smoother;

  vtkPolyData *smoothed_surfaces = smoother.smooth_surfaces(
      surfaces, kLambda);

  for (int count = 1; count < kNumberOfSmoothing; count++) {
    printf("count = %d\n", count);

    vtkPolyData *new_smoothed_surfaces =
        smoother.smooth_surfaces(smoothed_surfaces, kLambda);
    smoothed_surfaces->Delete();

    smoothed_surfaces = new_smoothed_surfaces;
  }

  printf("Finished smoothing.\n");

  vtkPolyData *positive_region = vtkPolyData::New();
  vtkPolyData *negative_region = vtkPolyData::New();

  positive_region->CopyStructure(smoothed_surfaces);
  negative_region->CopyStructure(smoothed_surfaces);

  vtkSmartPointer<vtkDoubleArray> positive_array =
      vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> negative_array =
      vtkSmartPointer<vtkDoubleArray>::New();

  positive_array->SetNumberOfComponents(1);
  negative_array->SetNumberOfComponents(1);

  positive_array->SetNumberOfTuples(surfaces->GetNumberOfCells());
  negative_array->SetNumberOfTuples(surfaces->GetNumberOfCells());

  for (int cell_id = 0;
       cell_id < smoothed_surfaces->GetNumberOfCells();
       cell_id++) {
    smoothed_surfaces->GetCellData()->SetActiveScalars("first_face");
    positive_array->SetTuple1(cell_id, smoothed_surfaces->GetCellData()
                                                        ->GetScalars()
                                                        ->GetTuple1(cell_id));

    smoothed_surfaces->GetCellData()->SetActiveScalars("second_face");
    negative_array->SetTuple1(cell_id, smoothed_surfaces->GetCellData()
                                                        ->GetScalars()
                                                        ->GetTuple1(cell_id));
  }

  positive_region->GetCellData()->SetScalars(positive_array);
  negative_region->GetCellData()->SetScalars(negative_array);

  SurfaceShader shader;
  vtkPolyData *colored_surfaces = shader.get_colored_surfaces(
      positive_region, negative_region);

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

  if (surfaces) {
    surfaces->Delete();
  }

  printf("} extract_surfaces_with_regions_test_2\n\n");
}
*/

int main() {
  extract_surfaces_test();
  // extract_surfaces_with_regions_test();
  // extract_surfaces_with_regions_test_2();

  return 0;
}
