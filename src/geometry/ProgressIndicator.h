// -*- mode: c++ -*-

#pragma once

#include <string>

// Dummmy class

class ProgressIndicator {
public:
  void processEvents() {}
  bool wasCancelled() const { return false; }

  void reset() {}
  void setLabelText(std::string) {}
  void setProgress(size_t) {}
  void setTotalSteps(size_t) {}

  void hide() {}
  void show() {}
};
