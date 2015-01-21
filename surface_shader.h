// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef SURFACE_SHADER_H_
#define SURFACE_SHADER__H_

class vtkPolyData;

class SurfaceShader {
 public:
  void shade_by_region(
      vtkPolyData *surfaces,
      vtkPolyData **positive_region,
      vtkPolyData **negative_region);

  vtkPolyData *get_colored_surfaces(
      vtkPolyData *positive_region,
      vtkPolyData *negative_region);
};

#endif  // SURFACE_SHADER_H_
