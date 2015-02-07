// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef SURFACE_EXTRACTOR_H_
#define SURFACE_EXTRACTOR_H_

class vtkStructuredPoints;
class vtkPolyData;

class SurfaceExtractor {
 public:
  // Region numbers are consecutive 0-based integers.
  vtkPolyData *extract_surfaces(vtkStructuredPoints *ftle,
                                vtkStructuredPoints *basins);

  // The positive face of a triangle is in counter clockwise order.
  void extract_surfaces_with_regions(vtkStructuredPoints *basins,
                                     vtkPolyData **positive_region,
                                     vtkPolyData **negative_region);

  // The first face of a triangle is in counter clockwise order.
  vtkPolyData *extract_surfaces_with_regions(
      vtkStructuredPoints *scalar_field,
      vtkStructuredPoints *basin_field);
};

#endif  // SURFACE_EXTRACTOR_H_
