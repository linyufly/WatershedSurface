// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef SURFACE_SMOOTHER_H_
#define SURFACE_SMOOTHER_H_

class vtkPolyData;

class SurfaceSmoother {
 public:
  vtkPolyData *smooth_surfaces(vtkPolyData *surfaces, double lambda);
};

#endif  // SURFACE_SMOOTHER_H_
