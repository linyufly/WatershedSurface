// Author: Mingcheng Chen (linyufly@gmail.com)

#include "sphere_generator.h"

#include <cmath>

#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

namespace {

double sqr(double a) {
  return a * a;
}

double distance(double x1, double y1, double z1,
                double x2, double y2, double z2) {
  return sqrt(sqr(x1 - x2) + sqr(y1 - y2) + sqr(z1 - z2));
}

};

void SphereGenerator::generate_two_spheres(
    int nx, int ny, int nz, double radius,
    vtkStructuredPoints **scalar_field,
    vtkStructuredPoints **basins) {
  double spacing[3] = {5 * radius / (nx - 1),
                       4 * radius / (ny - 1),
                       4 * radius / (nz - 1)};
  int dimensions[3] = {nx, ny, nz};
  double origin[3] = {0.0, 0.0, 0.0};
  double c1_x = 1.75 * radius;
  double c2_x = 3.25 * radius;
  double c_y = 2 * radius;
  double c_z = 2 * radius;

  *scalar_field = vtkStructuredPoints::New();
  (*scalar_field)->SetDimensions(dimensions);
  (*scalar_field)->SetOrigin(origin);
  (*scalar_field)->SetSpacing(spacing);

  *basins = vtkStructuredPoints::New();
  (*basins)->CopyStructure(*scalar_field);

  vtkSmartPointer<vtkDoubleArray> scalars =
      vtkSmartPointer<vtkDoubleArray>::New();
  scalars->SetName("ftle");
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(nx * ny * nz);

  vtkSmartPointer<vtkIntArray> codes =
      vtkSmartPointer<vtkIntArray>::New();
  codes->SetName("region");
  codes->SetNumberOfComponents(1);
  codes->SetNumberOfTuples(nx * ny * nz);

  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        int index = (z * ny + y) * nx + x;
        double coord_x = x * spacing[0];
        double coord_y = y * spacing[1];
        double coord_z = z * spacing[2];

        double d1 = distance(coord_x, coord_y, coord_z, c1_x, c_y, c_z);
        double d2 = distance(coord_x, coord_y, coord_z, c2_x, c_y, c_z);

        double d = (fabs(d1 - radius) < fabs(d2 - radius)) ?
            fabs(d1 - radius) : fabs(d2 - radius);

        scalars->SetTuple1(index, 1 / (d + 1));

        int region = -1;
        if (d1 >= radius && d2 >= radius) {
          region = 0;
        } else if (d1 >= radius && d2 < radius) {
          region = 1;
        } else if (d1 < radius && d2 >= radius) {
          region = 2;
        } else {
          region = 3;
        }

        codes->SetTuple1(index, region);
      }
    }
  }

  (*scalar_field)->GetPointData()->SetScalars(scalars);
  (*basins)->GetPointData()->SetScalars(codes);
}
