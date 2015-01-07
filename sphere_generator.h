// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef SPHERE_GENERATOR_H_
#define SPHERE_GENERATOR_H_

class vtkStructuredPoints;

class SphereGenerator {
 public:
  // The resulting field is 2R x 2R x 2R.
  void generate_one_sphere(int nx, int ny, int nz, double radius,
                           vtkStructuredPoints **scalar_field,
                           vtkStructuredPoints **basins);

  // The resulting field is 5R x 4R x 4R.
  void generate_two_spheres(int nx, int ny, int nz, double radius,
                            vtkStructuredPoints **scalar_field,
                            vtkStructuredPoints **basins);
};

#endif  // SPHERE_GENERATOR_H_
