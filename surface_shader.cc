// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_shader.h"

#include "lcsGeometry.h"
#include "util.h"

#include <cstdio>

#include <algorithm>
#include <map>
#include <vector>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

const double kGap = 0.00001;

namespace {

int get_positive(int cell_id) {
  return cell_id * 2;
}

int get_negative(int cell_id) {
  return cell_id * 2 + 1;
}

void get_point_ids(vtkCellArray *cell_array, int cell_id, int *point_ids) {
  vtkIdList *id_list = vtkIdList::New();

  cell_array->GetCell(cell_id * 4, id_list);

  if (id_list->GetNumberOfIds() != 3) {
    report_error("Wrong cell size at cell %d: %d\n",
                 cell_id, id_list->GetNumberOfIds());
  }

  for (int loc = 0; loc < 3; loc++) {
    point_ids[loc] = id_list->GetId(loc);
  }

  id_list->Delete();
}

}

void SurfaceShader::shade_by_region(
    vtkPolyData *surfaces,
    vtkPolyData **positive_region,
    vtkPolyData **negative_region) {
  vtkCellArray *cells = surfaces->GetPolys();
  vtkPoints *points = surfaces->GetPoints();

  int num_cells = cells->GetNumberOfCells();
  int num_points = points->GetNumberOfPoints();

  /// DEBUG ///
  printf("num_cells = %d\n", num_cells);
  printf("num_points = %d\n", num_points);

  std::vector<std::vector<int> > point_connectivity(num_points);

  /// DEBUG ///
  int pp[3];
  get_point_ids(cells, 20, pp);
  printf("%d %d %d\n", pp[0], pp[1], pp[2]);

  for (int cell_index = 0; cell_index < num_cells; cell_index++) {
    int point_ids[3];
    get_point_ids(cells, cell_index, point_ids);

    for (int loc = 0; loc < 3; loc++) {
      point_connectivity[point_ids[loc]].push_back(cell_index);
    }
  }

  /// DEBUG ///
  printf("after building point_connectivity\n");

  int *colors = new int[num_cells * 2];
  int *queue = new int[num_cells * 2];

  for (int face_index = 0; face_index < num_cells * 2; face_index++) {
    colors[face_index] = -1;
  }

  int num_colors = 0;

  for (int cell_index = 0; cell_index < num_cells; cell_index++) {
    int head = 0, tail = -1;
    if (colors[get_positive(cell_index)] == -1) {
      colors[get_positive(cell_index)] = num_colors++;
      queue[++tail] = get_positive(cell_index);
    }
    if (colors[get_negative(cell_index)] == -1) {
      colors[get_negative(cell_index)] = num_colors++;
      queue[++tail] = get_negative(cell_index);
    }
    while (head <= tail) {
      int curr_face_id = queue[head++];
      int curr_color = colors[curr_face_id];
      int curr_cell_id = curr_face_id / 2;
      int curr_point_ids[3];
      get_point_ids(cells, curr_cell_id, curr_point_ids);
      if (curr_face_id & 1) {
        std::swap(curr_point_ids[1], curr_point_ids[2]);
      }

      for (int loc = 0; loc < 3; loc++) {
        int curr_point_id = curr_point_ids[loc];
        int prev_point_id = curr_point_ids[(loc + 2) % 3];
        int succ_point_id = curr_point_ids[(loc + 1) % 3];

        for (int conn_itr = 0;
             conn_itr < static_cast<int>(
                 point_connectivity[curr_point_id].size());
             conn_itr++) {
          int next_cell_id = point_connectivity[curr_point_id][conn_itr];
          if (next_cell_id == curr_cell_id) {
            continue;
          }
          int next_point_ids[3];
          get_point_ids(cells, next_cell_id, next_point_ids);
          int pos = 0;
          for (; next_point_ids[pos] != curr_point_id; pos++) {}
          int new_prev_point_id = next_point_ids[(pos + 2) % 3];
          int new_succ_point_id = next_point_ids[(pos + 1) % 3];

          if ((new_prev_point_id == succ_point_id
              || new_succ_point_id == prev_point_id)
              && colors[get_positive(next_cell_id)] == -1) {
            colors[get_positive(next_cell_id)] = curr_color;
            queue[++tail] = get_positive(next_cell_id);
          }

          if ((new_prev_point_id == prev_point_id
              || new_succ_point_id == succ_point_id)
              && colors[get_negative(next_cell_id)] == -1) {
            colors[get_negative(next_cell_id)] = curr_color;
            queue[++tail] = get_negative(next_cell_id);
          }
        }
      }
    }
  }

  *positive_region = vtkPolyData::New();
  *negative_region = vtkPolyData::New();

  (*positive_region)->DeepCopy(surfaces);
  (*negative_region)->DeepCopy(surfaces);

  vtkSmartPointer<vtkIntArray> pos_regions =
      vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> neg_regions =
      vtkSmartPointer<vtkIntArray>::New();

  pos_regions->SetNumberOfComponents(1);
  neg_regions->SetNumberOfComponents(1);

  pos_regions->SetNumberOfTuples(num_cells);
  neg_regions->SetNumberOfTuples(num_cells);

  for (int cell_index = 0; cell_index < num_cells; cell_index++) {
    pos_regions->SetTuple1(cell_index, colors[get_positive(cell_index)]);
    neg_regions->SetTuple1(cell_index, colors[get_negative(cell_index)]);
  }

  (*positive_region)->GetCellData()->SetScalars(pos_regions);
  (*negative_region)->GetCellData()->SetScalars(neg_regions);

  delete [] queue;
  delete [] colors;
}

vtkPolyData *SurfaceShader::get_colored_surfaces(
      vtkPolyData *positive_region,
      vtkPolyData *negative_region) {
  vtkPolyData *colored_surfaces =
      vtkPolyData::New();

  vtkCellArray *cells = positive_region->GetPolys();
  vtkPoints *points = positive_region->GetPoints();

  int num_cells = cells->GetNumberOfCells();
  int num_points = points->GetNumberOfPoints();

  vtkSmartPointer<vtkPoints> colored_points =
      vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> colored_cells =
      vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkIntArray> colors =
      vtkSmartPointer<vtkIntArray>::New();
  colors->SetNumberOfComponents(1);
  colors->SetNumberOfTuples(positive_region->GetPolys()->GetNumberOfCells() * 2);

  for (int cell_index = 0; cell_index < num_cells; cell_index++) {
    int point_ids[3];
    get_point_ids(cells, cell_index, point_ids);

    double point_array[3];
    lcs::Vector curr_points[3];
    for (int loc = 0; loc < 3; loc++) {
      int point_id = point_ids[loc];
      points->GetPoint(point_id, point_array);
      curr_points[loc] = lcs::Vector(point_array);
    }

    lcs::Vector normal =
        lcs::Cross(curr_points[1] - curr_points[0],
                   curr_points[2] - curr_points[0]);

    if (normal.Length() < 0.000001) {
      // report_error("Small length of a vector: %lf", normal.Length());
      normal = lcs::Vector(0.0, 0.0, 0.0);
    } else {
      normal = normal / normal.Length();
      normal = normal * kGap;
    }

    // Add positive triangle
    for (int loc = 0; loc < 3; loc++) {
      colored_points->InsertNextPoint(curr_points[loc].GetX() + normal.GetX(),
                                      curr_points[loc].GetY() + normal.GetY(),
                                      curr_points[loc].GetZ() + normal.GetZ());
    }

    colored_cells->InsertNextCell(3);
    for (int loc = 0; loc < 3; loc++) {
      colored_cells->InsertCellPoint(cell_index * 6 + loc);
    }

    colors->SetTuple1(cell_index * 2, positive_region->GetCellData()
                                                     ->GetScalars()
                                                     ->GetTuple1(cell_index));

    // Add negative triangle
    for (int loc = 0; loc < 3; loc++) {
      colored_points->InsertNextPoint(curr_points[loc].GetX() - normal.GetX(),
                                      curr_points[loc].GetY() - normal.GetY(),
                                      curr_points[loc].GetZ() - normal.GetZ());
    }

    colored_cells->InsertNextCell(3);
    for (int loc = 0; loc < 3; loc++) {
      colored_cells->InsertCellPoint(cell_index * 6 + 3 + loc);
    }

    colors->SetTuple1(cell_index * 2 + 1,
                      negative_region->GetCellData()
                                     ->GetScalars()
                                     ->GetTuple1(cell_index));
  }

  colored_surfaces->SetPolys(colored_cells);
  colored_surfaces->SetPoints(colored_points);
  colored_surfaces->GetCellData()->SetScalars(colors);

  return colored_surfaces;
}

