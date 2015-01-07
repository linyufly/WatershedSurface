// Author: Mingcheng Chen (linyufly@gmail.com)

#include "sphere_generator.h"

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

const char *kBasinFile1 = "one_sphere_basins.vtk";
const char *kFTLEFile1 = "one_sphere_ftle.vtk";

const char *kBasinFile2 = "two_spheres_basins.vtk";
const char *kFTLEFile2 = "two_spheres_ftle.vtk";

void generate_one_sphere_test() {
  printf("generate_one_sphere_test {\n");

  vtkStructuredPoints *basins = NULL;
  vtkStructuredPoints *scalar = NULL;

  SphereGenerator generator;
  generator.generate_one_sphere(100, 100, 100, 1.0, &scalar, &basins);

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();

  writer->SetFileName(kBasinFile1);
  writer->SetInputData(basins);
  writer->Write();

  writer->SetFileName(kFTLEFile1);
  writer->SetInputData(scalar);
  writer->Write();

  if (basins) {
    basins->Delete();
  }

  if (scalar) {
    scalar->Delete();
  }

  printf("} generate_one_sphere_test\n\n");
}

void generate_two_spheres_test() {
  printf("generate_two_spheres_test {\n");

  vtkStructuredPoints *basins = NULL;
  vtkStructuredPoints *ftle = NULL;

  SphereGenerator generator;
  generator.generate_two_spheres(200, 160, 160, 1.0, &ftle, &basins);

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
      vtkSmartPointer<vtkStructuredPointsWriter>::New();

  writer->SetFileName(kBasinFile2);
  writer->SetInputData(basins);
  writer->Write();

  writer->SetFileName(kFTLEFile2);
  writer->SetInputData(ftle);
  writer->Write();

  if (basins) {
    basins->Delete();
  }

  if (ftle) {
    ftle->Delete();
  }

  printf("} generate_two_spheres_test\n\n");
}

int main() {
  generate_one_sphere_test();
  generate_two_spheres_test();

  return 0;
}
