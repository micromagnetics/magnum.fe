#ifndef _METRIC_MATRIX_CUBE_2_
#define _METRIC_MATRIX_CUBE_2_

#include <dolfin.h>

namespace magnumfe {

  class MetricMatrixCube2: public dolfin::Expression
  {
    protected:
      const dolfin::Array<double> size;
      const int  coord;

    public:
      MetricMatrixCube2(const dolfin::Array<double>& size, const int coord): dolfin::Expression(3, 3), size(size), coord(coord) {}
      ~MetricMatrixCube2(){}

      // returns the biggest of three double values
      static double max3(double a, double b, double c);

      // returns the biggest of three double values
      static double min3(double a, double b, double c);

      // alias Power for mathematica C export
      static double Power(double x, int y);

      void eval(dolfin::Array<double>& values, const dolfin::Array<double>& pos) const;
  };

};

#endif
