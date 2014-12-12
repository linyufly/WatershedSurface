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
};

#endif  // SURFACE_EXTRACTOR_H_
