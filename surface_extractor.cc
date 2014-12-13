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

// Every face is given in clockwise order.
const int kFaceList[6][4] = {
    {0, 1, 2, 3}, {0, 4, 5, 1}, {0, 3, 7, 4},
    {1, 5, 6, 2}, {2, 6, 7, 3}, {4, 7, 6, 5}
};

// vtx_1 and vtx_2 must share the same edge.
bool reorder_edge_points(int *vtx_1, int *vtx_2, int *dim) {
  for (*dim = 0; *dim < 3; (*dim)++) {
    if (kVertexList[*vtx_1][*dim] != kVertexList[*vtx_2][*dim]) {
      break;
    }
  }

  if (*dim >= 3) {
    report_error("In reorder_edge_points, *dim >= 3.\n");
  }

  if (kVertexList[*vtx_1][*dim] > kVertexList[*vtx_2][*dim]) {
    std::swap(*vtx_1, *vtx_2);
    return true;
  }

  return false;
}

int insert_face_point(int x, int y, int z, int face_id, int ****face_mark,
                      double point_x, double point_y, double point_z,
                      vtkPoints *mesh_points) {
  switch (face_id) {
    case 3: {
      x++;
      face_id = 2;
    } break;
    case 4: {
      y++;
      face_id = 1;
    } break;
    case 5: {
      z++;
      face_id = 0;
    } break;
  }

  if (face_mark[x][y][z][face_id] == -1) {
    face_mark[x][y][z][face_id] = mesh_points->GetNumberOfPoints();
    mesh_points->InsertNextPoint(point_x, point_y, point_z);
  }

  return face_mark[x][y][z][face_id];
}

// alpha is the ratio from (x, y, z) to the other end.
int insert_edge_point_basic_1(
    int x, int y, int z, int dim, int ****edge_mark,
    double *origin, double *spacing, double alpha,
    vtkPoints *mesh_points) {
  if (edge_mark[x][y][z][dim] == -1) {
    edge_mark[x][y][z][dim] = mesh_points->GetNumberOfPoints();

    int coordinates[3] = {x, y, z};
    double point[3];
    for (int i = 0; i < 3; i++) {
      point[i] = origin[i] + spacing[i] * coordinates[i];
    }

    point[dim] += spacing[dim] * alpha;

    mesh_points->InsertNextPoint(point[0], point[1], point[2]);
  }

  return edge_mark[x][y][z][dim];
}

int insert_edge_point_basic_2(
    int x, int y, int z, int point_id, int dim,
    int ****edge_mark, double *origin, double *spacing,
    double alpha, vtkPoints *mesh_points) {
  return insert_edge_point_basic_1(
      x + kVertexList[point_id][0],
      y + kVertexList[point_id][1],
      z + kVertexList[point_id][2],
      dim, edge_mark, origin, spacing, alpha,
      mesh_points);
}

// alpha is the ratio from vtx_1 to vtx_2.
int insert_edge_point(int x, int y, int z, int vtx_1, int vtx_2,
                       int ****edge_mark, double *origin, double *spacing,
                       double alpha, vtkPoints *mesh_points) {
  int dim = -1;

  if (reorder_edge_points(&vtx_1, &vtx_2, &dim)) {
    alpha = 1.0 - alpha;
  }
 
  return insert_edge_point_basic_2(x, y, z, vtx_1, dim,
                                   edge_mark, origin, spacing, 0.5,
                                   mesh_points);
}

// This method also adds the id of the edge point into the current cell of
// mesh_cells.
int insert_edge_point(int x, int y, int z, int vtx_1, int vtx_2,
                       int ****edge_mark, double *origin, double *spacing,
                       double alpha, vtkPoints *mesh_points,
                       vtkCellArray *mesh_cells) {
  int point_id = insert_edge_point(x, y, z, vtx_1, vtx_2,
                                   edge_mark, origin, spacing, 0.5,
                                   mesh_points);

  mesh_cells->InsertCellPoint(point_id);

  return point_id;
}

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
  int ****face_mark = create_4d_array<int>(nx, ny, nz, 3);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        for (int k = 0; k < 3; k++) {
          edge_mark[x][y][z][k] = -1;
          face_mark[x][y][z][k] = -1;
        }
      }
    }
  }

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

        /// DEBUG ///
        bool too_small = false;

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

              /// DEBUG ///
              if (ftle->GetPointData()->GetScalars()->GetTuple1(point_id) < 0.06) {
                too_small = true;
              }
            }
          }
        }

        /// DEBUG ///
        //if (too_small) {
        //  continue;
        //}

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

              insert_edge_point(x, y, z, vtx_1, vtx_2,
                                edge_mark, origin, spacing, 0.5,
                                mesh_points, mesh_cells);
            }
          }
        } else {  // Need in-cell point and face point
            double center_x = origin[0] + spacing[0] * (x + 0.5);
            double center_y = origin[1] + spacing[1] * (y + 0.5);
            double center_z = origin[2] + spacing[2] * (z + 0.5);
            int center_id = mesh_points->GetNumberOfPoints();
            mesh_points->InsertNextPoint(center_x, center_y, center_z);

            for (int face_id = 0; face_id < 6; face_id++) {
              std::set<int> face_code_sets;
              int first_code = -1;
              bool local_code[4];
              int local_non_binary_code[4];

              for (int i = 0; i < 4; i++) {
                int point_id = kFaceList[face_id][i];
                int dx = kVertexList[point_id][0];
                int dy = kVertexList[point_id][1];
                int dz = kVertexList[point_id][2];
                face_code_sets.insert(code[dx][dy][dz]);
                if (i == 0) {
                  first_code = code[dx][dy][dz];
                }
                local_code[i] = code[dx][dy][dz] == first_code;
                local_non_binary_code[i] = code[dx][dy][dz];
              }

              if (face_code_sets.size() == 1) {
                continue;
              }

              if (face_code_sets.size() == 2) {
                int num_true = local_code[0] + local_code[1]
                               + local_code[2] + local_code[3];
                if (num_true == 3 || num_true == 1) {  // 1 + 3
                  int single = -1;
                  if (num_true == 1) {
                    single = 0;
                  } else {
                    for (int i = 0; i < 4; i++) {
                      if (!local_code[i]) {
                        single = i;
                        break;
                      }
                    }
                  }

                  int prev_id = kFaceList[face_id][(single + 3) % 4];
                  int succ_id = kFaceList[face_id][(single + 1) % 4];
                  int curr_id = kFaceList[face_id][single];

                  mesh_cells->InsertNextCell(3);
                  mesh_cells->InsertCellPoint(center_id);
                  insert_edge_point(x, y, z, prev_id, curr_id,
                      edge_mark, origin, spacing, 0.5,
                      mesh_points, mesh_cells);
                  insert_edge_point(x, y, z, succ_id, curr_id,
                      edge_mark, origin, spacing, 0.5,
                      mesh_points, mesh_cells);

                  continue;
                } else if (num_true == 2
                           && (local_code[0] == local_code[1]
                               || local_code[1] == local_code[2])) {
                  // AA    BB
                  // BB or AA
                  mesh_cells->InsertNextCell(3);
                  mesh_cells->InsertCellPoint(center_id);
                  for (int i = 0; i < 4; i++) {
                    if (local_code[i] != local_code[(i + 1) % 4]) {
                      int curr_id = kFaceList[face_id][i];
                      int succ_id = kFaceList[face_id][(i + 1) % 4];
                      insert_edge_point(x, y, z, curr_id, succ_id,
                                        edge_mark, origin, spacing, 0.5,
                                        mesh_points, mesh_cells);
                    }
                  }
                  continue;
                }
              }

              double face_point_x = 0.0;
              double face_point_y = 0.0;
              double face_point_z = 0.0;
              for (int i = 0; i < 4; i++) {
                int point_id = kFaceList[face_id][i];
                face_point_x += spacing[0] * kVertexList[point_id][0];
                face_point_y += spacing[1] * kVertexList[point_id][1];
                face_point_z += spacing[2] * kVertexList[point_id][2];
              }

              face_point_x = origin[0] + spacing[0] * x + face_point_x / 4;
              face_point_y = origin[1] + spacing[1] * y + face_point_y / 4;
              face_point_z = origin[2] + spacing[2] * z + face_point_z / 4;

              int face_point_id = insert_face_point(
                  x, y, z, face_id, face_mark,
                  face_point_x, face_point_y, face_point_z,
                  mesh_points);

              for (int i = 0; i < 4; i++) {
                if (local_non_binary_code[i] != local_non_binary_code[(i + 1) % 4]) {
                  int curr_id = kFaceList[face_id][i];
                  int succ_id = kFaceList[face_id][(i + 1) % 4];
                  mesh_cells->InsertNextCell(3);
                  mesh_cells->InsertCellPoint(center_id);
                  mesh_cells->InsertCellPoint(face_point_id);
                  insert_edge_point(x, y, z, curr_id, succ_id,
                                    edge_mark, origin, spacing, 0.5,
                                    mesh_points, mesh_cells);
                }
              }
            }
        }
      }
    }
  }

  vtkPolyData *mesh = vtkPolyData::New();
  mesh->SetPoints(mesh_points);
  mesh->SetPolys(mesh_cells);

  delete_matrix(edge_mark);
  delete_matrix(face_mark);

  return mesh;
}
