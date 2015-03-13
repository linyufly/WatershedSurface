// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_extractor.h"

#include "marchingCubesTable.h"
#include "util.h"

#include <algorithm>
#include <set>
#include <vector>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

namespace {

const int kNumberOfBisections = 8;

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

double triple_product(double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3) {
  return x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 * z2
         - x1 * y3 * z2 - x2 * y1 * z3 - x3 * y2 * z1;
}

double sign_of_face(double viewer_x, double viewer_y, double viewer_z,
                    double x1, double y1, double z1,
                    double x2, double y2, double z2,
                    double x3, double y3, double z3) {
  return triple_product(x1 - viewer_x, y1 - viewer_y, z1 - viewer_z,
                        x2 - viewer_x, y2 - viewer_y, z2 - viewer_z,
                        x3 - viewer_x, y3 - viewer_y, z3 - viewer_z);
}

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

void get_grid_point(int x, int y, int z,
                    double *origin, double *spacing,
                    double *coordinates) {
  int index[3] = {x, y, z};
  for (int dim = 0; dim < 3; dim++) {
    coordinates[dim] = origin[dim] + spacing[dim] * index[dim];
  }
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

void get_coefficients_for_3d_equation(const double inv_a[3][3],
                                      const double values[4],
                                      double coefficients[3],
                                      double *right_value) {
  for (int index_unknown = 0; index_unknown < 3; index_unknown++) {
    coefficients[index_unknown] = 0.0;
  }

  for (int index_value = 1; index_value <= 3; index_value++) {
    for (int index_unknown = 0; index_unknown < 3; index_unknown++) {
      coefficients[index_unknown] +=
          inv_a[index_value - 1][index_unknown]
          * (values[index_value] - values[0]);
    }
  }

  *right_value = -values[0];
}

void clamp(int high, int *val) {
  if (*val < 0) {
    *val += 2;
  }

  if (*val >= high) {
    *val -= 2;
  }
}

double get_scalar(vtkStructuredPoints *scalar_field,
                  int nx, int ny, int nz, int x, int y, int z) {
  return scalar_field->GetPointData()
                     ->GetScalars()
                     ->GetTuple1((z * ny + y) * nx + x);
}

void get_gradients(int x, int y, int z,
                   vtkStructuredPoints *scalar_field,
                   double gradients[2][2][2][3]) {
  int dimensions[3];
  double spacing[3];

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  scalar_field->GetDimensions(dimensions);
  scalar_field->GetSpacing(spacing);

  for (int dx = 0; dx < 2; dx++) {
    for (int dy = 0; dy < 2; dy++) {
      for (int dz = 0; dz < 2; dz++) {
        int global_x = x + dx;
        int global_y = y + dy;
        int global_z = z + dz;

        int other_x = dx == 0 ? global_x - 1 : global_x + 1;
        int other_y = dy == 0 ? global_y - 1 : global_y + 1;
        int other_z = dz == 0 ? global_z - 1 : global_z + 1;

        clamp(nx, &other_x);
        clamp(ny, &other_y);
        clamp(nz, &other_z);

        double curr_value = get_scalar(
            scalar_field, nx, ny, nz, global_x, global_y, global_z);

        gradients[dx][dy][dz][0] = (get_scalar(
            scalar_field, nx, ny, nz, other_x, global_y, global_z)
            - curr_value) / ((other_x - global_x) * spacing[0]);

        gradients[dx][dy][dz][1] = (get_scalar(
            scalar_field, nx, ny, nz, global_x, other_y, global_z)
            - curr_value) / ((other_y - global_y) * spacing[1]);

        gradients[dx][dy][dz][2] = (get_scalar(
            scalar_field, nx, ny, nz, global_x, global_y, other_z)
            - curr_value) / ((other_z - global_z) * spacing[2]);
      }
    }
  }
}

// Mesh: (1) No edge is shared by more than two faces;
//       (2) All the faces are connected by edges.

// We do not need to set face point on minimal point. They
// should in the segment connected by two center points.

double get_gradient_norm_3d(double gradients[2][2][2][3],
                            double dx, double dy, double dz) {
  double gradient[3];
  for (int dimension = 0; dimension < 3; dimension++) {
    double g00 = gradients[0][0][0][dimension] * (1.0 - dz)
                 + gradients[0][0][1][dimension] * dz;
    double g01 = gradients[0][1][0][dimension] * (1.0 - dz)
                 + gradients[0][1][1][dimension] * dz;
    double g10 = gradients[1][0][0][dimension] * (1.0 - dz)
                 + gradients[1][0][1][dimension] * dz;
    double g11 = gradients[1][1][0][dimension] * (1.0 - dz)
                 + gradients[1][1][1][dimension] * dz;

    double g0 = g00 * (1.0 - dy) + g01 * dy;
    double g1 = g10 * (1.0 - dy) + g11 * dy;

    gradient[dimension] = g0 * (1.0 - dx) + g1 * dx;
  }

  double norm = 0.0;
  for (int dimension = 0; dimension < 3; dimension++) {
    norm += gradient[dimension] * gradient[dimension];
  }

  return norm;
}

void get_center_point(double gradients[2][2][2][3],
                      double *dx, double *dy, double *dz) {
  double lower_x = 0.0, upper_x = 1.0;
  double lower_y = 0.0, upper_y = 1.0;
  double lower_z = 0.0, upper_z = 1.0;

  for (int iteration_id = 0; iteration_id < kNumberOfBisections;
       iteration_id++) {
    double middle_x = (lower_x + upper_x) / 2.0;
    double middle_y = (lower_y + upper_y) / 2.0;
    double middle_z = (lower_z + upper_z) / 2.0;

    double new_lower_x = -1.0, new_upper_x = -1.0;
    double new_lower_y = -1.0, new_upper_y = -1.0;
    double new_lower_z = -1.0, new_upper_z = -1.0;

    double minimal = -1.0;

    for (int block_x = 0; block_x < 2; block_x++) {
      for (int block_y = 0; block_y < 2; block_y++) {
        for (int block_z = 0; block_z < 2; block_z++) {
          double curr_lower_x = block_x == 0 ? lower_x : middle_x;
          double curr_upper_x = block_x == 0 ? middle_x : upper_x;
          double curr_lower_y = block_y == 0 ? lower_y : middle_y;
          double curr_upper_y = block_y == 0 ? middle_y : upper_y;
          double curr_lower_z = block_z == 0 ? lower_z : middle_z;
          double curr_upper_z = block_z == 0 ? middle_z : upper_z;

          double curr_x = (curr_lower_x + curr_upper_x) / 2.0;
          double curr_y = (curr_lower_y + curr_upper_y) / 2.0;
          double curr_z = (curr_lower_z + curr_upper_z) / 2.0;

          double magnitude = get_gradient_norm_3d(gradients, curr_x, curr_y, curr_z);

          if (minimal < 0.0 || minimal > magnitude) {
            minimal = magnitude;

            new_lower_x = curr_lower_x;
            new_upper_x = curr_upper_x;

            new_lower_y = curr_lower_y;
            new_upper_y = curr_upper_y;

            new_lower_z = curr_lower_z;
            new_upper_z = curr_upper_z;
          }
        }
      }
    }

    lower_x = new_lower_x;
    upper_x = new_upper_x;

    lower_y = new_lower_y;
    upper_y = new_upper_y;

    lower_z = new_lower_z;
    upper_z = new_upper_z;
  }

  *dx = (lower_x + upper_x) / 2.0;
  *dy = (lower_y + upper_y) / 2.0;
  *dz = (lower_z + upper_z) / 2.0;
}

void add_a_single_triangle(int single, int face_id, int center_id,
                           int x, int y, int z,
                           double origin[3], double spacing[3],
                           int local_non_binary_code[4],
                           std::vector<int> *positive_face,
                           std::vector<int> *negative_face,
                           vtkCellArray *mesh_cells,
                           vtkPoints *mesh_points,
                           int ****edge_mark) {
  int negative_color = local_non_binary_code[single];
  int positive_color = local_non_binary_code[(single + 1) % 4];
  positive_face->push_back(positive_color);
  negative_face->push_back(negative_color);

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
}

}

vtkPolyData *SurfaceExtractor::extract_surfaces(
    vtkStructuredPoints *ftle, vtkStructuredPoints *basins) {
  int dimensions[3];
  double spacing[3], origin[3];
  basins->GetDimensions(dimensions);
  basins->GetSpacing(spacing);
  basins->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  /// DEBUG ///
  printf("nx, ny, nz: %d, %d, %d\n", nx, ny, nz);

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
        // bool too_small = false;

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
              // if (ftle->GetPointData()->GetScalars()->GetTuple1(point_id) < 0.06) {
              //   too_small = true;
              // }
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

void SurfaceExtractor::extract_surfaces_with_regions(
    vtkStructuredPoints *basins,
    vtkPolyData **positive_region,
    vtkPolyData **negative_region) {
  int dimensions[3];
  double spacing[3], origin[3];
  basins->GetDimensions(dimensions);
  basins->GetSpacing(spacing);
  basins->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  /// DEBUG ///
  printf("nx, ny, nz: %d, %d, %d\n", nx, ny, nz);

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

  std::vector<int> positive_markers;
  std::vector<int> negative_markers;

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
        /* if (code_sets.size() == 2) {  // Similar to Marching Cubes
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
        } else */ {  // Need in-cell point and face point
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

                  int positive_color = local_non_binary_code[single];
                  int negative_color = local_non_binary_code[(single + 1) % 4];
                  positive_markers.push_back(positive_color);
                  negative_markers.push_back(negative_color);

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

                  bool first = true;

                  for (int i = 0; i < 4; i++) {
                    if (local_code[i] != local_code[(i + 1) % 4]) {
                      if (first) {
                        first = false;
                        int positive_color =
                            local_non_binary_code[(i + 1) % 4];
                        int negative_color =
                            local_non_binary_code[i];
                        positive_markers.push_back(positive_color);
                        negative_markers.push_back(negative_color);
                      }

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

                  // The negative_color and positive_color have not been tested.
                  int negative_color = local_non_binary_code[(i + 1) % 4];
                  int positive_color = local_non_binary_code[i];
                  positive_markers.push_back(positive_color);
                  negative_markers.push_back(negative_color);

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

  /// DEBUG ///
  printf("positive_markers.size() = %d\n",
         static_cast<int>(positive_markers.size()));
  printf("negative_markers.size() = %d\n",
         static_cast<int>(negative_markers.size()));

  vtkSmartPointer<vtkIntArray> positive_array =
      vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> negative_array =
      vtkSmartPointer<vtkIntArray>::New();

  positive_array->SetNumberOfComponents(1);
  negative_array->SetNumberOfComponents(1);

  positive_array->SetNumberOfTuples(positive_markers.size());
  negative_array->SetNumberOfTuples(negative_markers.size());

  for (int cell_index = 0;
       cell_index < static_cast<int>(positive_markers.size());
       cell_index++) {
    positive_array->SetTuple1(cell_index, positive_markers[cell_index]);
    negative_array->SetTuple1(cell_index, negative_markers[cell_index]);
  }

  *positive_region = vtkPolyData::New();
  *negative_region = vtkPolyData::New();

  (*positive_region)->SetPoints(mesh_points);
  (*negative_region)->SetPoints(mesh_points);

  (*positive_region)->SetPolys(mesh_cells);
  (*negative_region)->SetPolys(mesh_cells);

  (*positive_region)->GetCellData()->SetScalars(positive_array);
  (*negative_region)->GetCellData()->SetScalars(negative_array);

  delete_matrix(edge_mark);
  delete_matrix(face_mark);
}

/*
vtkPolyData *SurfaceExtractor::extract_surfaces_with_regions(
      vtkStructuredPoints *scalar_field,
      vtkStructuredPoints *basin_field) {
  int dimensions[3];
  double spacing[3], origin[3];
  basin_field->GetDimensions(dimensions);
  basin_field->GetSpacing(spacing);
  basin_field->GetOrigin(origin);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  /// DEBUG ///
  printf("nx, ny, nz: %d, %d, %d\n", nx, ny, nz);

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

  std::vector<int> positive_face;
  std::vector<int> negative_face;

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
              code[dx][dy][dz] = basin_field->GetPointData()
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

            double triangle[3][3];

            for (int j = 0; j < 3; j++) {
              int edge_idx = triTable[cube_code][i + j];
              int vtx_1 = kEdgeList[edge_idx][0];
              int vtx_2 = kEdgeList[edge_idx][1];

              int point_id = insert_edge_point(x, y, z, vtx_1, vtx_2,
                                               edge_mark, origin, spacing, 0.5,
                                               mesh_points, mesh_cells);

              mesh_points->GetPoint(point_id,
                                    triangle[j]);
            }

            // Examine the first vertex of the triangle
            int edge_idx = triTable[cube_code][i];
            int vtx_1 = kEdgeList[edge_idx][0];
            int vtx_2 = kEdgeList[edge_idx][1];
            int dx_1 = kVertexList[vtx_1][0];
            int dy_1 = kVertexList[vtx_1][1];
            int dz_1 = kVertexList[vtx_1][2];
            int dx_2 = kVertexList[vtx_2][0];
            int dy_2 = kVertexList[vtx_2][1];
            int dz_2 = kVertexList[vtx_2][2];
            double grid_point_1[3];
            get_grid_point(x + dx_1, y + dy_1, z + dz_1,
                           origin, spacing, grid_point_1);
            double sign = sign_of_face(
                grid_point_1[0], grid_point_1[1], grid_point_1[2],
                triangle[0][0], triangle[0][1], triangle[0][2],
                triangle[1][0], triangle[1][1], triangle[1][2],
                triangle[2][0], triangle[2][1], triangle[2][2]);

            // Assign positive and negative faces
            if (sign < 0.0) {
              negative_face.push_back(code[dx_1][dy_1][dz_1]);
              positive_face.push_back(code[dx_2][dy_2][dz_2]);
            } else {
              negative_face.push_back(code[dx_2][dy_2][dz_2]);
              positive_face.push_back(code[dx_1][dy_1][dz_1]);
            }

            if (code[dx_1][dy_1][dz_1] == code[dx_2][dy_2][dz_2]) {
              report_error("A triangle separates the same region.");
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

                  // int negative_color = local_non_binary_code[single];
                  // int positive_color = local_non_binary_code[(single + 1) % 4];
                  positive_face.push_back(positive_color);
                  negative_face.push_back(negative_color);

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

                  bool first = true;
                  for (int i = 0; i < 4; i++) {
                    if (local_code[i] != local_code[(i + 1) % 4]) {
                      if (first) {
                        first = false;
                        int negative_color =
                            local_non_binary_code[(i + 1) % 4];
                        int positive_color =
                            local_non_binary_code[i];
                        positive_face.push_back(positive_color);
                        negative_face.push_back(negative_color);
                      }

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
                if (local_non_binary_code[i]
                    != local_non_binary_code[(i + 1) % 4]) {
                  int positive_color = local_non_binary_code[(i + 1) % 4];
                  int negative_color = local_non_binary_code[i];
                  positive_face.push_back(positive_color);
                  negative_face.push_back(negative_color);

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

  // Calculate locations of mesh points

  vtkSmartPointer<vtkDoubleArray> positive_face_array =
      vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> negative_face_array =
      vtkSmartPointer<vtkDoubleArray>::New();

  positive_face_array->SetNumberOfComponents(1);
  negative_face_array->SetNumberOfComponents(1);

  positive_face_array->SetNumberOfTuples(positive_face.size());
  negative_face_array->SetNumberOfTuples(negative_face.size());
  
  positive_face_array->SetName("second_face");
  negative_face_array->SetName("first_face");

  for (int cell_index = 0;
       cell_index < static_cast<int>(positive_face.size());
       cell_index++) {
    positive_face_array->SetTuple1(cell_index, positive_face[cell_index]);
    negative_face_array->SetTuple1(cell_index, negative_face[cell_index]);
  }

  vtkPolyData *mesh = vtkPolyData::New();
  mesh->SetPoints(mesh_points);
  mesh->SetPolys(mesh_cells);
  mesh->GetCellData()->AddArray(positive_face_array);
  mesh->GetCellData()->AddArray(negative_face_array);

  delete_matrix(edge_mark);
  delete_matrix(face_mark);

  return mesh;
}
*/
