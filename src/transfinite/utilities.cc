double
blendHermite(double x) {
  double x2 = x * x;
  return 2.0 * x * x2 - 3.0 * x2 + 1.0;
}
