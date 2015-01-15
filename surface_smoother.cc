// Author: Mingcheng Chen (linyufly@gmail.com)

#include "surface_smoother.h"

#include "util.h"

#include <cstdio>

#include <map>
#include <vector>

#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>

vtkPolyData *SurfaceSmoother::smooth_surfaces(
    vtkPolyData *surfaces, double lambda) {
  vtkCellArray *cells = surfaces->GetPolys();
  vtkPoints *points = surfaces->GetPoints();

  int num_points = points->GetNumberOfPoints();
  std::vector<std::map<int, int> > point_link(num_points);

  int num_cells = cells->GetNumberOfCells();

  /// DEBUG ///
  // printf("num_points = %d\n", num_points);
  // printf("num_cells = %d\n", num_cells);

  cells->InitTraversal();
  for (int cell_index = 0; cell_index < num_cells; cell_index++) {
    vtkIdList *id_list = vtkIdList::New();
    cells->GetNextCell(id_list);

    if (id_list->GetNumberOfIds() != 3) {
      report_error("Wrong cell size at cell %d: %d\n",
                   cell_index, id_list->GetNumberOfIds());
    }

    int cell_ids[3];
    for (int id_loc = 0; id_loc < id_list->GetNumberOfIds(); id_loc++) {
      cell_ids[id_loc] = id_list->GetId(id_loc);
    }

    for (int id_loc = 0; id_loc < 3; id_loc++) {
      int curr_id = cell_ids[id_loc];
      int next_id = cell_ids[(id_loc + 1) % 3];
      point_link[curr_id][next_id]++;
      point_link[next_id][curr_id]++;
    }

    id_list->Delete();
  }

  double **new_position = create_matrix<double>(num_points, 3);

  for (int point_id = 0; point_id < num_points; point_id++) {
    int max_repetitions = 0;
    for (std::map<int, int>::iterator visitor = point_link[point_id].begin();
         visitor != point_link[point_id].end(); visitor++) {
      if (visitor->second > max_repetitions) {
        max_repetitions = visitor->second;
      }
    }

    double curr_point[3];
    points->GetPoint(point_id, curr_point);

    double accumulation[3] = {0.0, 0.0, 0.0};
    int count = 0;

    for (std::map<int, int>::iterator visitor = point_link[point_id].begin();
         visitor != point_link[point_id].end(); visitor++) {
      if (visitor->second == max_repetitions) {
        int neighbor_id = visitor->first;
        double neighbor_point[3];
        points->GetPoint(neighbor_id, neighbor_point);
        for (int dimension = 0; dimension < 3; dimension++) {
          accumulation[dimension] +=
              neighbor_point[dimension] - curr_point[dimension];
          count++;
        }
      }
    }

    for (int dimension = 0; dimension < 3; dimension++) {
      new_position[point_id][dimension] =
          curr_point[dimension] + lambda * accumulation[dimension] / count;
    }
  }

  vtkPolyData *smoothed_surfaces = vtkPolyData::New();
  smoothed_surfaces->DeepCopy(surfaces);

  /// DEBUG ///
  // printf("smoothed_surfaces->GetNumberOfPoints() = %d\n",
  //        static_cast<int>(smoothed_surfaces->GetNumberOfPoints()));
  // printf("smoothed_surfaces->GetNumberOfCells() = %d\n",
  //        static_cast<int>(smoothed_surfaces->GetNumberOfCells()));

  for (int point_id = 0; point_id < num_points; point_id++) {
    smoothed_surfaces->GetPoints()->SetPoint(point_id, new_position[point_id]);
  }

  delete_matrix(new_position);

  return smoothed_surfaces;
}
