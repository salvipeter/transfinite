#include "with-handler.hh"

namespace Transfinite {

WithHandler::WithHandler() : multiplier_(1.0), handler_initialized_(false) {
}

WithHandler::~WithHandler() {
}

void
WithHandler::setMultiplier(double m) {
  multiplier_ = m;
}

void
WithHandler::setHandler(const Vector3D &h) {
  handler_ = h;
  handler_initialized_ = true;
}

void
WithHandler::resetHandler() {
  multiplier_ = 1.0;
  handler_initialized_ = true;
}

} // namespace Transfinite
