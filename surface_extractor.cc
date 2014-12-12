// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_extractor.h"

#include "marchingCubesTable.h"
#include "util.h"

#include <set>
#include <algorithm>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

namespace {

const int kVertexList[8][3] = {
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
};

const int kEdgeList[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

}

vtkPolyData *SurfaceExtractor::extract_surfaces(
    vtkStructuredPoints *ftle, vtkStructuredPoints *basins) {
  int dimensions[3];
  double spacing[3], origin[3];
  ftle->GetDimensions(dimensions);
  ftle->GetSpacing(spacing);
  ftle->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int ****edge_mark = create_4d_array<int>(nx, ny, nz, 3);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        for (int k = 0; k < 3; k++) {
          edge_mark[x][y][z][k] = -1;
        }
      }
    }
  }

  vtkSmartPointer<vtkPoints> mesh_points =
      vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> mesh_cells =
      vtkSmartPointer<vtkCellArray>::New();

  /// DEBUG ///
  int cnt = 0;
  int max_region = 0;

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
          int cube_code = 0;
          for (int i = 0; i < 8; i++) {
            int dx = kVertexList[i][0];
            int dy = kVertexList[i][1];
            int dz = kVertexList[i][2];
            if (code[dx][dy][dz] == code[0][0][0]) {
              cube_code |= 1 << i;
            }
          }
          for (int i = 0; i < numVertsTable[cube_code]; i += 3) {
            mesh_cells->InsertNextCell(3);

            for (int j = 0; j < 3; j++) {
              int edge_idx = triTable[cube_code][i + j];
              int vtx_1 = kEdgeList[edge_idx][0];
              int vtx_2 = kEdgeList[edge_idx][1];

              int dim;
              for (dim = 0; dim < 3; dim++) {
                if (kVertexList[vtx_1][dim] != kVertexList[vtx_2][dim]) {
                  break;
                }
              }

              if (kVertexList[vtx_1][dim] > kVertexList[vtx_2][dim]) {
                std::swap(vtx_1, vtx_2);
              }

              int start_x = kVertexList[vtx_1][0];
              int start_y = kVertexList[vtx_1][1];
              int start_z = kVertexList[vtx_1][2];

              int finish_x = kVertexList[vtx_2][0];
              int finish_y = kVertexList[vtx_2][1];
              int finish_z = kVertexList[vtx_2][2];

              // Insert a new point to the mesh if necessary
              if (edge_mark[x + start_x][y + start_y][z + start_z][dim]
                  == -1) {
                edge_mark[x + start_x][y + start_y][z + start_z][dim] =
                    mesh_points->GetNumberOfPoints();

                double aug_x = (finish_x - start_x) * spacing[0] / 2;
                double aug_y = (finish_y - start_y) * spacing[1] / 2;
                double aug_z = (finish_z - start_z) * spacing[2] / 2;

                double point_x = origin[0] + spacing[0] * (x + start_x)
                                 + aug_x;
                double point_y = origin[1] + spacing[1] * (y + start_y)
                                 + aug_y;
                double point_z = origin[2] + spacing[2] * (z + start_z)
                                 + aug_z;

                mesh_points->InsertNextPoint(point_x, point_y, point_z);
              }

              mesh_cells->InsertCellPoint(
                  edge_mark[x + start_x][y + start_y][z + start_z][dim]);
            }
          }
        } else {  // Need in-cell point and face point
          /// DEBUG ///
          cnt++;
          if (code_sets.size() > max_region) {
            max_region = code_sets.size();
          }
        }
      }
    }
  }

  /// DEBUG ///
  printf("cnt = %d\n", cnt);
  printf("max_region = %d\n", max_region);

  vtkPolyData *mesh = vtkPolyData::New();
  mesh->SetPoints(mesh_points);
  mesh->SetPolys(mesh_cells);

  delete_matrix(edge_mark);

  return mesh;
}
