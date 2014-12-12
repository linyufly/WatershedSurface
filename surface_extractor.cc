// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_extractor.h"

#include <set>

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

vtkPolyData *SurfaceExtractor::extract_surfaces(
    vtkStructuredPoints *ftle, vtkStructuredPoints *basins,
    int num_regions) {
  int dimensions[3];
  double spacing[3], origin[3];
  ftle->GetDimensions(dimensions);
  ftle->GetSpacing(spacing);
  ftle->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int ****edge_mark = create_4d_array<int>(nx, ny, nz, 3);

  vtkSmartPointer<vtkPoints> mesh_points =
      vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> mesh_cells =
      vtkSmartPointer<vtkCellArray>::New();

  // Get mid-edge points
  for (int x = 0; x + 1 < nx; x++) {
    for (int y = 0; y + 1 < ny; y++) {
      for (int z = 0; z + 1 < nz; z++) {
        int code[2][2][2];
        std::set<int> code_sets;
        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            for (int dz = 0; dz < 2; dz++) {
              int curr_x = x + dx;
              int curr_y = y + dy;
              int curr_z = z + dz;
              int point_id = (curr_z * ny + curr_y) * nx + curr_x;
              code[dx][dy][dz] = basins->GetPointData()
                                       ->GetScalars()
                                       ->GetTuple1(point_id);
              code_sets.insert(code[dx][dy][dz]);
            }
          }
        }
        if (code_sets.size() == 1) {
          continue;
        }
        if (code_sets.size() == 2) {  // Similar to Marching Cubes

        } else {  // Need in-cell point
        }
      }
    }
  }

  delete_matrix(edge_mark);
}
