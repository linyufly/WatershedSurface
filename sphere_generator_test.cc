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

const char *kBasinFile = "sphere_basins.vtk";
const char *kFTLEFile = "sphere_ftle.vtk";

void generate_two_spheres_test() {
  printf("generate_two_spheres_test {\n");

  vtkStructuredPoints *basins = NULL;
  vtkStructuredPoints *ftle = NULL;

  SphereGenerator generator;
  generator.generate_two_spheres(250, 200, 200, 1.0, &ftle, &basins);

  vtkSmartPointer<vtkStructuredPointsWriter> writer =
    vtkSmartPointer<vtkStructuredPointsWriter>::New();

  writer->SetFileName(kBasinFile);
  writer->SetInputData(basins);
  writer->Write();

  writer->SetFileName(kFTLEFile);
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
  generate_two_spheres_test();

  return 0;
}
