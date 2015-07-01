#pragma once

#include <ctime>

#include "geometry.hh"

#include "IgesExport.hh"

class IGES::IGESImpl {
public:
  IGESImpl(std::string filename) {
    std::time_t curr_time = std::time(nullptr);
    RCPtr<IgesGlobal> iges_global =
      IgesGlobal::create("Sketches", "Sketches", "1.0.2", "Product", "Product",
                         "Author", "ShapEx", "protocol", curr_time, 0.1, filename.c_str(),
                         IgesGlobal::IGU_MILLIMETERS, IgesGlobal::IGDSTD_NONE, 1.0, 5.0);
    filter_ = IgesExport(iges_global, IgesExport::IGVERSION_5_2);
  }
  void writeSurface(const Surface &surfaces) const { filter_ << surface; }
  void writeTrimmedSurface(const TrimmedSurface &surface) const { filter_ << surface; }
  void close() { filter_.close(); }
private:
  IgesExport filter_;
};
